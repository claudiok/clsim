/**
 * Copyright (c) 2013
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file I3MuonSliceRemoverAndPulseRelabeler.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3ParticleID.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/I3Constants.h"

#include "simclasses/I3MCPE.h"
#include "simclasses/I3Photon.h"
#include "simclasses/I3CompressedPhoton.h"
#include "simclasses/I3ParticleIDMap.hpp"

/**
 * @brief Removes muon slices in I3MCTree objects
 * as written by I3MuonSlicer and re-labels I3MCPEs generated from
 * these so that they refer back to the original muon. 
 *
 * Applying this after having generated I3MCPEs using clsim
 * will generate an MCTree fully compatible with what ppc
 * would generate (i.e. the light will look light it had been
 * created by the original muon)
 *
 * TODO: THIS IS NOT READY FOR USE YET!
 *
 */
class I3MuonSliceRemoverAndPulseRelabeler : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3MuonSliceRemoverAndPulseRelabeler(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    virtual ~I3MuonSliceRemoverAndPulseRelabeler();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
    /**
     * The module needs to process Physics frames
     */
    virtual void DAQ(I3FramePtr frame);

    
private:
    // parameters
    
    /// Parameter: Name of the input I3MCTree frame object as written by I3MuonSlicer. 
    std::string inputMCTreeName_;

    /// Parameter: Name of the original MCTree. Will be used for sanity checks
    std::string oldMCTreeName_;

    /// Parameter: Name of the input I3MCPESeriesMap frame object.
    std::string inputMCPESeriesMapName_;

    /// Parameter: Name of the output I3MCPESeriesMap frame object (can be identical to the input).
    std::string outputMCPESeriesMapName_;

    /// Parameter: Name of the input I3PhotonSeriesMap frame object.
    std::string inputPhotonSeriesMapName_;

    /// Parameter: Name of the output I3PhotonSeriesMap frame object (can be identical to the input).
    std::string outputPhotonSeriesMapName_;
    
private:
    // default, assignment, and copy constructor declared private
    I3MuonSliceRemoverAndPulseRelabeler();
    I3MuonSliceRemoverAndPulseRelabeler(const I3MuonSliceRemoverAndPulseRelabeler&);
    I3MuonSliceRemoverAndPulseRelabeler& operator=(const I3MuonSliceRemoverAndPulseRelabeler&);

    
    SET_LOGGER("I3MuonSliceRemoverAndPulseRelabeler");
};

// The module
I3_MODULE(I3MuonSliceRemoverAndPulseRelabeler);

I3MuonSliceRemoverAndPulseRelabeler::I3MuonSliceRemoverAndPulseRelabeler(const I3Context& context) 
: I3ConditionalModule(context)
{
    
    inputMCTreeName_="I3MCTree_sliced";
    AddParameter("InputMCTreeName",
                 "Name of the input I3MCTree frame object.",
                 inputMCTreeName_);

    oldMCTreeName_="I3MCTree";
    AddParameter("OldMCTreeName",
                 "Name of the original MCTree. Will be used for sanity checks",
                 oldMCTreeName_);

    inputMCPESeriesMapName_="MCPESeriesMap";
    AddParameter("InputMCPESeriesMapName",
                 "Name of the input I3MCPESeriesMap frame object.",
                 inputMCPESeriesMapName_);

    outputMCPESeriesMapName_="MCPESeriesMap";
    AddParameter("OutputMCPESeriesMapName",
                 "Name of the output I3MCPESeriesMap frame object.",
                 outputMCPESeriesMapName_);

    inputPhotonSeriesMapName_="";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object.",
                 inputPhotonSeriesMapName_);

    outputPhotonSeriesMapName_="";
    AddParameter("OutputPhotonSeriesMapName",
                 "Name of the output I3PhotonSeriesMap frame object.",
                 outputPhotonSeriesMapName_);


    // add an outbox
    AddOutBox("OutBox");

}

I3MuonSliceRemoverAndPulseRelabeler::~I3MuonSliceRemoverAndPulseRelabeler()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3MuonSliceRemoverAndPulseRelabeler::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("InputMCTreeName", inputMCTreeName_);
    GetParameter("OldMCTreeName", oldMCTreeName_);
    GetParameter("InputMCPESeriesMapName", inputMCPESeriesMapName_);
    GetParameter("OutputMCPESeriesMapName", outputMCPESeriesMapName_);
    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputPhotonSeriesMapName", outputPhotonSeriesMapName_);

    if (inputMCTreeName_=="")
        log_fatal("The \"InputMCTreeName\" parameter must not be empty.");

}

namespace {
    void RemoveMuonSlicesIfDarkMuon(const I3MCTree &mcTree,
                                    const I3Particle &particle,
                                    std::map<I3ParticleID, I3ParticleID> &re_label_map
                                    )
    {
        // is this a dark muon? (no, not *that* kind of dark muon, Kurt!)
        if ((particle.GetShape() == I3Particle::Dark) &&
            ((particle.GetType()==I3Particle::MuMinus) ||
             (particle.GetType()==I3Particle::MuPlus)))
        {
            // yes, it's a dark muon.

            for(const I3Particle& daughter : mcTree.children(particle))
            {
                if ((daughter.GetType() != I3Particle::MuMinus) &&
                    (daughter.GetType() != I3Particle::MuPlus))
                    continue; // it's not a muon

                if (daughter.GetShape() == I3Particle::Dark)
                    log_fatal("Dark muon has dark muon daughter particle. Your MCTree is messed up.");

                if ((bool) mcTree.first_child(daughter))
                    log_fatal("Non-dark daughter muon of dark muon has children of its own. Your MCTree is messed up.");

                // insert into map so we know how to re-label hits later on
                re_label_map.insert( std::pair<I3ParticleID, I3ParticleID>( daughter, particle ) );
            }

        }

        for(const I3Particle& daughter : mcTree.children(particle))
        {
            RemoveMuonSlicesIfDarkMuon(mcTree, daughter, re_label_map);
        }

    }

    ///\return whether the photons were found and relabeled
    template<typename PhotonType>
    bool RelabelPhotons(I3FramePtr frame,
                        const std::string& inputPhotonSeriesMapName,
                        const std::string& outputPhotonSeriesMapName,
                        const std::map<I3ParticleID, I3ParticleID>& re_label_map)
    {
        using PhotonSeriesMap=I3Map<ModuleKey, I3Vector<PhotonType> >;
        auto inputPhotonSeriesMap = frame->Get<boost::shared_ptr<const PhotonSeriesMap>>(inputPhotonSeriesMapName);
        if (!inputPhotonSeriesMap)
            return false;
        
        // allocate the output map (start with copies of the input objects)
        boost::shared_ptr<PhotonSeriesMap> outputPhotonSeriesMap(new PhotonSeriesMap(*inputPhotonSeriesMap));
        
        for (auto& dom_entry : *outputPhotonSeriesMap)
        {
            I3Vector<PhotonType>& photonSeries = dom_entry.second;
            
            for(PhotonType& photon : photonSeries)
            {
                I3ParticleID oldID(photon.GetParticleMajorID(), photon.GetParticleMinorID());
                
                auto re_label_it = re_label_map.find(oldID);
                if (re_label_it==re_label_map.end())
                    // no need to re-label
                    continue;
                
                const I3ParticleID newID = re_label_it->second;
                photon.SetParticleMajorID(newID.majorID);
                photon.SetParticleMinorID(newID.minorID);
            }
        }
        
        // store the output I3MCPhotonSeriesMap
        if (inputPhotonSeriesMapName==outputPhotonSeriesMapName) {
            frame->Delete(inputPhotonSeriesMapName);
        }
        frame->Put(outputPhotonSeriesMapName, outputPhotonSeriesMap);
        
        return true;
    }
}

void I3MuonSliceRemoverAndPulseRelabeler::DAQ(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3MCTreeConstPtr inputMCTree = frame->Get<I3MCTreeConstPtr>(inputMCTreeName_);
    if (!inputMCTree) {
        log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                  inputMCTreeName_.c_str());
        PushFrame(frame);
        return;
    }

    I3MCTreeConstPtr oldMCTree = frame->Get<I3MCTreeConstPtr>(oldMCTreeName_);
    if (!oldMCTree) {
        log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                  oldMCTreeName_.c_str());
        PushFrame(frame);
        return;
    }
    std::set<I3ParticleID> oldTree_IDs;
    for (I3MCTree::const_iterator j=oldMCTree->begin(); j!=oldMCTree->end(); ++j)
        oldTree_IDs.insert(*j);

    I3MCPESeriesMapPtr outputMCPESeriesMap;
    boost::shared_ptr<const I3ParticleIDMap> inputPEParentTable;

    if (inputMCPESeriesMapName_ != "")
    {
        I3MCPESeriesMapConstPtr inputMCPESeriesMap;
        inputMCPESeriesMap = frame->Get<I3MCPESeriesMapConstPtr>(inputMCPESeriesMapName_);
        if (!inputMCPESeriesMap) {
            log_fatal("Frame does not contain an I3MCPESeriesMap named \"%s\".",
                      inputMCPESeriesMapName_.c_str());
            PushFrame(frame);
            return;
        }

        // allocate the output map (start with copies of the input objects)
        outputMCPESeriesMap = I3MCPESeriesMapPtr(new I3MCPESeriesMap(*inputMCPESeriesMap));
        
        // check for a side table which may also exist
        inputPEParentTable=frame->Get<boost::shared_ptr<const I3ParticleIDMap>>(inputMCPESeriesMapName_+"ParticleIDMap");
    }

    // now just work on the output objects

    // oldID->newID
    std::map<I3ParticleID, I3ParticleID> re_label_map;

    // get a list of primaries
    const std::vector<const I3Particle*> primaries = I3MCTreeUtils::GetPrimariesPtr(inputMCTree);
    
    // add each one to the output tree and check their children
    for(const I3Particle* primary : primaries)
    {
        RemoveMuonSlicesIfDarkMuon(*inputMCTree, *primary, re_label_map);
    }
    
    // sanity check
    for (std::map<I3ParticleID, I3ParticleID>::const_iterator it = re_label_map.begin();
        it != re_label_map.end(); ++it)
    {
        const I3ParticleID &fromID = it->first;
        const I3ParticleID &toID = it->second;

        if (oldTree_IDs.count(toID)==0) 
            log_fatal("Original muon ID is not in \"old\" tree.");

        if (oldTree_IDs.count(fromID)!=0) {
            log_fatal("Muon slice ID *is* in \"old\" tree.");
        }
    }

    // re-label hits if requested
    if (!inputMCPESeriesMapName_.empty())
    {
        for (auto& it : *outputMCPESeriesMap)
        {
            I3MCPESeries &peSeries = it.second;

            for (I3MCPE& pe : peSeries)
            {
                auto re_label_it = re_label_map.find(pe.ID);
                if (re_label_it!=re_label_map.end())
                    pe.ID = re_label_it->second;
                // otherwise no need to relabel
            }
        }

        // store the output I3MCPESeriesMap
        if (inputMCPESeriesMapName_==outputMCPESeriesMapName_) {
            frame->Delete(inputMCPESeriesMapName_);
        }
        frame->Put(outputMCPESeriesMapName_, outputMCPESeriesMap);
        
        // if parent information is stored externally, fix it up too
        if (inputPEParentTable)
        {
            boost::shared_ptr<I3ParticleIDMap> outputPEParentTable(new I3ParticleIDMap);
            for (auto& omdata : *inputPEParentTable)
            {
                OMKey om = omdata.first;
                ParticlePulseIndexMap& dest = (*outputPEParentTable)[om];
                //relabel and probably merge the index lists
                for (auto& p : omdata.second)
                {
                    I3ParticleID key = p.first;
                    auto re_label_it = re_label_map.find(key);
                    if (re_label_it==re_label_map.end())
                        key = re_label_it->second;
                    dest[key].insert(dest[key].end(),p.second.begin(),p.second.end());
                }
                //sort index lists and remove duplicates merging may have introduced
                for (auto& p : dest)
                {
                    std::sort(p.second.begin(),p.second.end());
                    std::vector<uint32_t> temp;
                    std::unique_copy(p.second.begin(),p.second.end(),std::back_inserter(temp));
                    p.second.swap(temp);
                }
            }
            if (inputMCPESeriesMapName_==outputMCPESeriesMapName_) {
                frame->Delete(inputMCPESeriesMapName_+"ParticleIDMap");
            }
            frame->Put(outputMCPESeriesMapName_+"ParticleIDMap", outputPEParentTable);
        }
    }
    
    // re-label photons if requested
    if (!inputPhotonSeriesMapName_.empty())
    {
        //Try regular photons first
        bool done = RelabelPhotons<I3Photon>(frame, inputPhotonSeriesMapName_,
                                             outputPhotonSeriesMapName_, re_label_map);
        //Try compressed photons
        if (!done)
            done = RelabelPhotons<I3CompressedPhoton>(frame, inputPhotonSeriesMapName_,
                                            outputPhotonSeriesMapName_, re_label_map);
        //If we still didn't do it, something is wrong
        if (!done)
            log_fatal_stream("Frame does not contain photons named \""
                             << inputPhotonSeriesMapName_ << '"');
    }
    
    // that's it!
    PushFrame(frame);
}
