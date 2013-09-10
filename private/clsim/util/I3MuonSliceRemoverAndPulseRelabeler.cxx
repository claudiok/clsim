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
 * $Id: I3MuonSliceRemoverAndPulseRelabeler.cxx 108199 2013-07-12 21:33:08Z nwhitehorn $
 *
 * @file I3MuonSliceRemoverAndPulseRelabeler.cxx
 * @version $Revision: 108199 $
 * @date $Date: 2013-07-12 16:33:08 -0500 (Fri, 12 Jul 2013) $
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
#ifdef IS_Q_FRAME_ENABLED
    virtual void DAQ(I3FramePtr frame);
#else
    virtual void Physics(I3FramePtr frame);
#endif

    
private:
    // parameters
    
    /// Parameter: Name of the input I3MCTree frame object as written by I3MuonSlicer. 
    std::string inputMCTreeName_;

    /// Parameter: Name of the output I3MCTree frame object to generate (can be identical to the input). 
    std::string outputMCTreeName_;

    /// Parameter: Name of the original MCTree if available. Will be used for sanity checks
    std::string oldMCTreeForSanityCheckName_;

    /// Parameter: Name of the input I3MCPESeriesMap frame object.
    std::string inputMCPESeriesMapName_;

    /// Parameter: Name of the output I3MCPESeriesMap frame object (can be identical to the input).
    std::string outputMCPESeriesMapName_;
    
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

    outputMCTreeName_="I3MCTree_unsliced";
    AddParameter("OutputMCTreeName",
                 "Name of the output I3MCTree frame object.",
                 outputMCTreeName_);

    oldMCTreeForSanityCheckName_="I3MCTree";
    AddParameter("OldMCTreeForSanityCheckName",
                 "Name of the original MCTree if available. Will be used for sanity checks",
                 oldMCTreeForSanityCheckName_);

    inputMCPESeriesMapName_="MCPESeriesMap";
    AddParameter("InputMCPESeriesMapName",
                 "Name of the input I3MCPESeriesMap frame object.",
                 inputMCPESeriesMapName_);

    outputMCPESeriesMapName_="MCPESeriesMap";
    AddParameter("OutputMCPESeriesMapName",
                 "Name of the output I3MCPESeriesMap frame object.",
                 outputMCPESeriesMapName_);


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
    GetParameter("OutputMCTreeName", outputMCTreeName_);
    GetParameter("OldMCTreeForSanityCheckName", oldMCTreeForSanityCheckName_);
    GetParameter("InputMCPESeriesMapName", inputMCPESeriesMapName_);
    GetParameter("OutputMCPESeriesMapName", outputMCPESeriesMapName_);

    if (inputMCTreeName_=="")
        log_fatal("The \"InputMCTreeName\" parameter must not be empty.");

    if (inputMCPESeriesMapName_=="")
        log_fatal("The \"InputMCPESeriesMapName\" parameter must not be empty.");
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

    void RemoveMuonSlicesIfDarkMuon(I3MCTree &mcTree,
                                    I3MCTree::iterator particle_it,
                                    std::map<I3ParticleID, I3ParticleID> &re_label_map
                                    )
    {
        I3Particle &particle = *particle_it;

        // is this a dark muon? (no, not *that* kind of dark muon, Kurt!)
        if ((particle.GetShape() == I3Particle::Dark) &&
            ((particle.GetType()==I3Particle::MuMinus) ||
             (particle.GetType()==I3Particle::MuPlus)))
        {
            // yes, it's a dark muon. delete all of its children
            // that are a) muons, b) not dark and c) have no children
            // themselves

            std::vector<I3MCTree::iterator> daughterList;
            I3MCTree::sibling_iterator j(particle_it);
            for (j=mcTree.begin(particle_it); j!=mcTree.end(particle_it); ++j)
                daughterList.push_back(j);

            bool had_at_least_one_muon_daughter = false;

            I3Particle::ParticleShape resetToShape = I3Particle::InfiniteTrack;

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

                resetToShape = daughter.GetShape();

                // insert into map so we know how to re-label hits later on
                re_label_map.insert( std::pair<I3ParticleID, I3ParticleID>( daughter, particle ) );

                // it's all okay. remove it from the tree.
                mcTree.erase(daughter_it);

                had_at_least_one_muon_daughter=true;
            }

            if (!had_at_least_one_muon_daughter)
            {
                // log_error("Would have expected at least one non-dark muon daugther particle. Your MCTree is messed up.");
            }
            else
            {
                // now set the particle to "in-ice"
                particle.SetShape(resetToShape);
            }

        }
        else
        {
            // it's something else.
            // do nothing for now.
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

#ifdef IS_Q_FRAME_ENABLED
void I3MuonSliceRemoverAndPulseRelabeler::DAQ(I3FramePtr frame)
#else
void I3MuonSliceRemoverAndPulseRelabeler::Physics(I3FramePtr frame)
#endif
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3MCTreeConstPtr inputMCTree = frame->Get<I3MCTreeConstPtr>(inputMCTreeName_);
    if (!inputMCTree) {
        log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                  inputMCTreeName_.c_str());
        PushFrame(frame);
        return;
    }

    I3MCTreeConstPtr oldMCTreeForSanityCheck;
    if (oldMCTreeForSanityCheckName_ != "")
    {
        oldMCTreeForSanityCheck = frame->Get<I3MCTreeConstPtr>(oldMCTreeForSanityCheckName_);
        if (!oldMCTreeForSanityCheck) {
            log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                      oldMCTreeForSanityCheckName_.c_str());
            PushFrame(frame);
            return;
        }
    }

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

    // allocate the output objects (start with copies of the input objects)
    I3MCTreePtr outputMCTree(new I3MCTree(*inputMCTree));
    
    // now just work on the output objects

    // oldID->newID
    std::map<I3ParticleID, I3ParticleID> re_label_map;
    std::set<I3ParticleID> sanity_check_IDs;

    if (bool(oldMCTreeForSanityCheck)) 
    {
        for (I3MCTree::iterator j=oldMCTreeForSanityCheck->begin(); j!=oldMCTreeForSanityCheck->end(); ++j)
            sanity_check_IDs.insert(*j);
    }

    // get a list of primaries
    const std::vector<I3Particle> primaries = I3MCTreeUtils::GetPrimaries(outputMCTree);
    
    // add each one to the output tree and check their children
    BOOST_FOREACH(const I3Particle &primary, primaries)
    {
        I3MCTree::iterator primary_in_output_tree = GetMCTreeIterator(*outputMCTree, primary);

        RemoveMuonSlicesIfDarkMuon(*outputMCTree, primary_in_output_tree, re_label_map);
    }
    
    // store the output I3MCTree
    if (outputMCTreeName_==inputMCTreeName_) {
        frame->Delete(outputMCTreeName_);
        frame->Put(outputMCTreeName_, outputMCTree);
    } else if (outputMCTreeName_!="") {
        frame->Put(outputMCTreeName_, outputMCTree);
    }


    // TODO: This should als re-label I3Photons!


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

                if (sanity_check_IDs.size()>0)
                {
                    if (sanity_check_IDs.count(oldID) > 0)
                    {
                        log_fatal("particle which should not have been in the old tree has been found anyway");
                    }
                }

                const I3ParticleID newID = re_label_it->second;
                pe.major_ID = newID.majorID;
                pe.minor_ID = newID.minorID;

                // log_warn("changed ID");
            }
        }

        if (sanity_check_IDs.size()>0)
        {
            log_warn("checking..");

            // sanity check
            for (I3MCPESeriesMap::const_iterator it = outputMCPESeriesMap->begin();
                 it != outputMCPESeriesMap->end(); ++it)
            {
                const I3MCPESeries &peSeries = it->second;
                BOOST_FOREACH(const I3MCPE &pe, peSeries)
                {
                    I3ParticleID particleID;
                    particleID.majorID = pe.major_ID;
                    particleID.minorID = pe.minor_ID;

                    if (sanity_check_IDs.count(particleID)==0)
                        log_fatal("new tree and sanity check tree are incompatible!");
                }
            }
        }

        // store the output I3MCPESeriesMap
        if (inputMCPESeriesMapName_==outputMCPESeriesMapName_) {
            frame->Delete(inputMCPESeriesMapName_);
        }
        frame->Put(outputMCPESeriesMapName_, outputMCPESeriesMap);
    }
    
    
    // that's it!
    PushFrame(frame);
}
