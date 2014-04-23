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

#include <string>

#include <boost/foreach.hpp>

#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3ParticleID.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/I3TreeUtils.h"
#include "dataclasses/I3Constants.h"

#include "simclasses/I3MCPE.h"
#include "clsim/I3Photon.h"

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
    inline
    I3MCTree::iterator
    GetMCTreeIterator(const I3MCTree &t, const I3Particle& p){
        I3MCTree::iterator i;
        for(i=t.begin() ; i!= t.end(); i++)
            if((i->GetMinorID() == p.GetMinorID()) &&
               (i->GetMajorID() == p.GetMajorID()))
                return i;
        return t.end();
    }

    void RemoveMuonSlicesIfDarkMuon(const I3MCTree &mcTree,
                                    I3MCTree::iterator particle_it,
                                    std::map<I3ParticleID, I3ParticleID> &re_label_map
                                    )
    {
        const I3Particle &particle = *particle_it;

        // is this a dark muon? (no, not *that* kind of dark muon, Kurt!)
        if ((particle.GetShape() == I3Particle::Dark) &&
            ((particle.GetType()==I3Particle::MuMinus) ||
             (particle.GetType()==I3Particle::MuPlus)))
        {
            // yes, it's a dark muon.

            std::vector<I3MCTree::iterator> daughterList;
            I3MCTree::sibling_iterator j(particle_it);
            for (j=mcTree.begin(particle_it); j!=mcTree.end(particle_it); ++j)
                daughterList.push_back(j);

            BOOST_FOREACH(I3MCTree::iterator &daughter_it, daughterList)
            {
                const I3Particle &daughter = *daughter_it;

                if ((daughter.GetType() != I3Particle::MuMinus) &&
                    (daughter.GetType() != I3Particle::MuPlus))
                    continue; // it's not a muon

                if (daughter.GetShape() == I3Particle::Dark)
                    log_fatal("Dark muon has dark muon daughter particle. Your MCTree is messed up.");

                if (mcTree.number_of_children(daughter_it) > 0)
                    log_fatal("Non-dark daughter muon of dark muon has children of its own. Your MCTree is messed up.");

                // this needs to work on older versions of dataclasses without
                // "I3ParticleID I3Particle::operator()" - just make IDs and set
                // their properties by hand for now.
                I3ParticleID daughter_id;
                daughter_id.majorID = daughter.GetMajorID();
                daughter_id.minorID = daughter.GetMinorID();

                I3ParticleID particle_id;
                particle_id.majorID = particle.GetMajorID();
                particle_id.minorID = particle.GetMinorID();

                // insert into map so we know how to re-label hits later on
                re_label_map.insert( std::pair<I3ParticleID, I3ParticleID>( daughter_id, particle_id ) );
            }

        }

        // now decend the tree and recursively call this function again
        std::vector<I3MCTree::iterator> daughterList;
        I3MCTree::sibling_iterator j(particle_it);
        for (j=mcTree.begin(particle_it); j!=mcTree.end(particle_it); ++j)
            daughterList.push_back(j);

        BOOST_FOREACH(I3MCTree::iterator &daughter_it, daughterList)
        {
            RemoveMuonSlicesIfDarkMuon(mcTree, daughter_it, re_label_map);
        }

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
    for (I3MCTree::iterator j=oldMCTree->begin(); j!=oldMCTree->end(); ++j)
        oldTree_IDs.insert(*j);

    I3MCPESeriesMapConstPtr inputMCPESeriesMap;
    I3MCPESeriesMapPtr outputMCPESeriesMap;

    if (inputMCPESeriesMapName_ != "")
    {
        inputMCPESeriesMap = frame->Get<I3MCPESeriesMapConstPtr>(inputMCPESeriesMapName_);
        if (!inputMCPESeriesMap) {
            log_fatal("Frame does not contain an I3MCPESeriesMap named \"%s\".",
                      inputMCPESeriesMapName_.c_str());
            PushFrame(frame);
            return;
        }

        // allocate the output map (start with copies of the input objects)
        outputMCPESeriesMap = I3MCPESeriesMapPtr(new I3MCPESeriesMap(*inputMCPESeriesMap));
    }

    I3PhotonSeriesMapConstPtr inputPhotonSeriesMap;
    I3PhotonSeriesMapPtr outputPhotonSeriesMap;

    if (inputPhotonSeriesMapName_ != "")
    {
        inputPhotonSeriesMap = frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
        if (!inputPhotonSeriesMap) {
            log_fatal("Frame does not contain an I3PhotonSeriesMap named \"%s\".",
                      inputPhotonSeriesMapName_.c_str());
            PushFrame(frame);
            return;
        }

        // allocate the output map (start with copies of the input objects)
        outputPhotonSeriesMap = I3PhotonSeriesMapPtr(new I3PhotonSeriesMap(*inputPhotonSeriesMap));
    }

    // now just work on the output objects

    // oldID->newID
    std::map<I3ParticleID, I3ParticleID> re_label_map;

    // get a list of primaries
    const std::vector<I3Particle> primaries = I3MCTreeUtils::GetPrimaries(inputMCTree);
    
    // add each one to the output tree and check their children
    BOOST_FOREACH(const I3Particle &primary, primaries)
    {
        I3MCTree::iterator primary_in_output_tree = GetMCTreeIterator(*inputMCTree, primary);

        RemoveMuonSlicesIfDarkMuon(*inputMCTree, primary_in_output_tree, re_label_map);
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
    if (inputMCPESeriesMapName_ != "")
    {
        for (I3MCPESeriesMap::iterator it = outputMCPESeriesMap->begin();
             it != outputMCPESeriesMap->end(); ++it)
        {
            I3MCPESeries &peSeries = it->second;

            BOOST_FOREACH(I3MCPE &pe, peSeries)
            {
                I3ParticleID oldID;
                oldID.majorID = pe.major_ID;
                oldID.minorID = pe.minor_ID;

                std::map<I3ParticleID, I3ParticleID>::const_iterator re_label_it = re_label_map.find(oldID);
                if (re_label_it==re_label_map.end())
                    // no need to re-label
                    continue;

                const I3ParticleID newID = re_label_it->second;
                pe.major_ID = newID.majorID;
                pe.minor_ID = newID.minorID;
            }
        }

        // store the output I3MCPESeriesMap
        if (inputMCPESeriesMapName_==outputMCPESeriesMapName_) {
            frame->Delete(inputMCPESeriesMapName_);
        }
        frame->Put(outputMCPESeriesMapName_, outputMCPESeriesMap);
    }
    
    // re-label photons if requested
    if (inputPhotonSeriesMapName_ != "")
    {
        for (I3PhotonSeriesMap::iterator it = outputPhotonSeriesMap->begin();
             it != outputPhotonSeriesMap->end(); ++it)
        {
            I3PhotonSeries &photonSeries = it->second;

            BOOST_FOREACH(I3Photon &photon, photonSeries)
            {
                I3ParticleID oldID;
                oldID.majorID = photon.GetParticleMajorID();
                oldID.minorID = photon.GetParticleMinorID();

                std::map<I3ParticleID, I3ParticleID>::const_iterator re_label_it = re_label_map.find(oldID);
                if (re_label_it==re_label_map.end())
                    // no need to re-label
                    continue;

                const I3ParticleID newID = re_label_it->second;
                photon.SetParticleMajorID(newID.majorID);
                photon.SetParticleMinorID(newID.minorID);
            }
        }

        // store the output I3MCPhotonSeriesMap
        if (inputPhotonSeriesMapName_==outputPhotonSeriesMapName_) {
            frame->Delete(inputPhotonSeriesMapName_);
        }
        frame->Put(outputPhotonSeriesMapName_, outputPhotonSeriesMap);
    }
    
    // that's it!
    PushFrame(frame);
}
