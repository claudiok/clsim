/**
 * Copyright (c) 2011, 2012
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
 * @file I3MuonSlicer.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "clsim/util/I3MuonSlicer.h"

#include <boost/foreach.hpp>

#include <gsl/gsl_sys.h>

#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/I3TreeUtils.h"
#include "dataclasses/I3Constants.h"

#include "simclasses/I3MMCTrack.h"

// The module
I3_MODULE(I3MuonSlicer);

I3MuonSlicer::I3MuonSlicer(const I3Context& context) 
: I3ConditionalModule(context)
{
    
    inputMCTreeName_="I3MCTree";
    AddParameter("InputMCTreeName",
                 "Name of the I3MCTree frame object.",
                 inputMCTreeName_);

    inputMCTreeName_="MMCTrackList";
    AddParameter("MMCTrackListName",
                 "Name of the I3MMCTrackList frame object.",
                 MMCTrackListName_);

    inputMCTreeName_="I3MCTree_sliced";
    AddParameter("OutputMCTreeName",
                 "Name of the output I3MCTree frame object. If identical to the\n"
                 "input or empty, the input object will be replaced.",
                 outputMCTreeName_);


    // add an outbox
    AddOutBox("OutBox");

}

I3MuonSlicer::~I3MuonSlicer()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3MuonSlicer::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("InputMCTreeName", inputMCTreeName_);
    GetParameter("MMCTrackListName", MMCTrackListName_);
    GetParameter("OutputMCTreeName", outputMCTreeName_);

    if (inputMCTreeName_=="")
        log_fatal("The \"InputMCTreeName\" parameter must not be empty.");
}

// helper functions
namespace {
    inline bool IsAnyOfType(const std::vector<I3MCTree::iterator> &particles, I3Particle::ParticleType type)
    {
        BOOST_FOREACH(const I3MCTree::iterator &particle_it, particles)
        {
            const I3Particle &particle = *particle_it;
            if (particle.GetType()==type) return true;
        }
        return false;
    }
    
    inline bool DoesItHaveChildren(const I3MCTree &mcTree, const I3Particle &particle)
    {
        const std::vector<I3Particle> daughters =
        I3MCTreeUtils::GetDaughters(mcTree, particle);
        return (daughters.size() > 0);
    }
    
    inline bool AreParticlesSortedInTime(const std::vector<I3MCTree::iterator> &particles)
    {
        bool firstIt=true;
        double previousTime=NAN;
        
        BOOST_FOREACH(const I3MCTree::iterator &particle_it, particles)
        {
            const I3Particle &particle = *particle_it;

            if (firstIt) {
                if (isnan(particle.GetTime())) return false; // not sorted..
                previousTime=particle.GetTime();
                continue;
            }

            if (particle.GetTime() < previousTime) return false; // not sorted..
            
            previousTime = particle.GetTime();
        }
        return true;
    }
    
    inline double GetTotalEnergyOfParticles(const std::vector<I3MCTree::iterator> &particles,
                                            double tMin, double tMax)
    {
        double totalEnergy=0.;
        
        BOOST_FOREACH(const I3MCTree::iterator &particle_it, particles)
        {
            const I3Particle &particle = *particle_it;

            if ((particle.GetTime() < tMin) || (particle.GetTime() > tMax))
                continue;

            totalEnergy+=particle.GetEnergy();
        }
        
        return totalEnergy;
    }

    inline double DistanceFromInfiniteTrack(const I3Position &position, const I3Position &muonPos, const I3Direction &muonDir, double &distanceOnTrack)
    {
        double pos0_x = muonPos.GetX();
        double pos0_y = muonPos.GetY();
        double pos0_z = muonPos.GetZ();
        
        double e_x = muonDir.GetX();
        double e_y = muonDir.GetY();
        double e_z = muonDir.GetZ();
        
        double h_x = position.GetX() - pos0_x;  //vector between particle position and "OM"/"cascade"
        double h_y = position.GetY() - pos0_y;
        double h_z = position.GetZ() - pos0_z;
        
        double s = e_x*h_x + e_y*h_y + e_z*h_z; //distance between particle position and closest approach position
        
        double pos2_x = pos0_x + s*e_x; //closest approach position
        double pos2_y = pos0_y + s*e_y;
        double pos2_z = pos0_z + s*e_z;
        
        distanceOnTrack=s; //record the distance on the track
        
        return gsl_hypot3(pos2_x-position.GetX(), pos2_y-position.GetY(), pos2_z-position.GetZ());
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

        inline 
        const std::vector<I3MCTree::iterator> 
        GetDaughterIteratorsFromParentIterator(const I3MCTree& t, I3MCTree::iterator parent_it)
        {
            std::vector<I3MCTree::iterator> daughterList;
            I3MCTree::sibling_iterator j(parent_it);
            for (j=t.begin(parent_it); j!=t.end(parent_it); ++j)
                daughterList.push_back(j);
            return daughterList;
        }
    }
    
    void SliceMuonOrCopySubtree(const I3MCTree &inputTree,
                                const I3MMCTrackList &mmcTrackList,
                                const std::map<std::pair<uint64_t, int>,
                                const I3MMCTrack *> &mmcTrackListIndex,
                                I3MCTree &outputTree,
                                const I3MCTree::iterator &particle_it_inputTree,
                                I3MCTree::iterator particle_it_outputTree
                                )
    {
        const I3Particle &particle = *particle_it_inputTree;

        // cache the iterator to speed up adding things to the tree
        if (particle_it_inputTree == inputTree.end()) log_fatal("internal error. output particle not in output tree.");
        if (particle_it_outputTree == outputTree.end()) log_fatal("internal error. output particle not in output tree.");

        const std::vector<I3MCTree::iterator> daughterIterators =
        GetDaughterIteratorsFromParentIterator(inputTree, particle_it_inputTree);

        // special treatment for muons with a length only
        if (((!isnan(particle.GetLength())) && (particle.GetLength() > 0.)) &&
            ((particle.GetType()==I3Particle::MuMinus) ||
            (particle.GetType()==I3Particle::MuPlus)))
        {
            // is any of the daughters a muon?
            if (IsAnyOfType(daughterIterators, I3Particle::MuMinus) ||
                IsAnyOfType(daughterIterators, I3Particle::MuPlus) ||
                IsAnyOfType(daughterIterators, I3Particle::unknown))
            {
                log_fatal("It seems you either ran MMC with the \"-recc\" option or I3MuonSlicer has already been applied.");
            }

            if (std::abs(particle.GetSpeed()-I3Constants::c) > 1e-5)
                log_fatal("Found a muon that does not travel with the speed of light. v=%gm/ns",
                          particle.GetSpeed()/(I3Units::m/I3Units::ns));
            
            if (isnan(particle.GetEnergy())) log_fatal("Muon must have an energy");
            if (isnan(particle.GetTime())) log_fatal("Muon must have a time");
            
            // get MMC track times and energies
            double ti=NAN;
            double Ei=NAN;
            double tf=NAN;
            double Ef=NAN;

            // find it in the I3MMCTrackList
            std::map<std::pair<uint64_t, int>, const I3MMCTrack *>::const_iterator it =
            mmcTrackListIndex.find(std::make_pair(particle.GetMajorID(), particle.GetMinorID()));
            if (it==mmcTrackListIndex.end())
            {
                log_debug("Muon is not in I3MMCTrackList. (minorID=%i, majorID=%" PRIu64 ") length=%fm.",
                          particle.GetMinorID(), particle.GetMajorID(), particle.GetLength()/I3Units::m);

                ti=NAN;
                Ei=NAN;
                tf=NAN;
                Ef=NAN;
            }
            else
            {
                const I3MMCTrack &mmcTrack = *(it->second);
                
                // get MMC track times and energies
                ti = mmcTrack.GetTi();
                Ei = mmcTrack.GetEi();
                tf = mmcTrack.GetTf();
                Ef = mmcTrack.GetEf();
            }
            
            // daughters need to be sorted in time (ascending)
            if (!AreParticlesSortedInTime(daughterIterators))
            {
                log_fatal("Muon daughters are not sorted in time (ascending).");
            }

            bool hadInvalidEi=false;
            bool hadInvalidEf=false;
            
            // correct values with information from I3Particle
            if ((Ei<=0.) || (isnan(Ei)))
            {
                Ei = particle.GetEnergy();
                ti = particle.GetTime();
                hadInvalidEi=true;
            }
            if ((Ef<0.) || (isnan(Ef)))
            {
                Ef = 0.;
                tf = particle.GetTime() + particle.GetLength()/I3Constants::c;
                hadInvalidEf=true;
            }

            if (isnan(ti)) log_fatal("t_initial is NaN");
            if (isnan(tf)) log_fatal("t_final is NaN");

            if (Ei<=0) 
            {
                if (daughterIterators.size() > 0)
                    log_fatal("Muon with Energy==0 has children.");
            }
            else if (tf<ti)
            {
                log_warn("Muon stops before it starts.. [Setting to shape \"Dark\"] (t_final==%fns < t_initial==%fns) (GetTime()=%fns, minorID=%i, majorID=%" PRIu64 ") ignoring.",
                         tf/I3Units::ns, ti/I3Units::ns, particle.GetTime()/I3Units::ns, particle.GetMinorID(), particle.GetMajorID());

                // set the current muon shape to "Dark"
                particle_it_outputTree->SetShape(I3Particle::Dark);
            }
            else if (tf==ti)
            {
                if (particle.GetLength()>0.)
                    log_warn("Particle has length but t_final==t_initial==%fns. [Setting to shape \"Dark\"] (GetTime()=%fns, minorID=%i, majorID=%" PRIu64 ") ignoring.",
                              tf/I3Units::ns, particle.GetTime()/I3Units::ns, particle.GetMinorID(), particle.GetMajorID());

                // set the current muon shape to "Dark"
                particle_it_outputTree->SetShape(I3Particle::Dark);
            }
            else // (ti<tf) && (Ei>0.)
            {
                // set the current muon shape to "Dark"
                particle_it_outputTree->SetShape(I3Particle::Dark);
                
                const double totalEnergyInCascades = 
                GetTotalEnergyOfParticles(daughterIterators, ti, tf);
                const double dEdt_calc = ((Ef-Ei+totalEnergyInCascades)/(ti-tf));
                const double dEdt_max = (0.21+8.8e-3*log(Ei/I3Units::GeV)/log(10.))*(I3Units::GeV/I3Units::m)*I3Constants::c;  // for ice only at Ecut=500 MeV (stolen from PPC)
                const double dEdt = std::min(dEdt_calc,dEdt_max);
                
                // add all daughters to the muon (in the output tree)
                // while inserting the muon slices between them
                
                double currentEnergy = Ei;
                double currentTime = ti;
                
                unsigned int iterationNum=0;
                
                BOOST_FOREACH(const I3MCTree::iterator daughter_it, daughterIterators)
                {
                    // make a copy of the particle here (we might need to change it later)
                    I3Particle daughter = *daughter_it;

                    if (currentEnergy<0.) {
                        log_error("Muon loses more energy than it has. Ecurrent=%gGeV, Ei=%fGeV, now reset to E=0", currentEnergy/I3Units::GeV, Ei/I3Units::GeV);
                        currentEnergy=0.;
                    }

                    double distanceOnTrack=NAN;
                    const double distanceFromMuonTrack = DistanceFromInfiniteTrack(daughter.GetPos(), particle.GetPos(), particle.GetDir(), distanceOnTrack);
                    if (distanceFromMuonTrack > 1.*I3Units::mm) {
                        log_error("cascade is not on muon track! (distance from (infinite) track=%gm). Not splitting the muon.",
                                  distanceFromMuonTrack/I3Units::m);

                        // append it to the output tree anyway
                        outputTree.append_child(particle_it_outputTree,daughter);
                        continue;
                    }
                    
                    if (isnan(daughter.GetTime())) continue;
                    
                    // calculate an expected time for the cascade (from its position on the track)
                    const double expectedTime = particle.GetTime() + distanceOnTrack/I3Constants::c;
                    
                    if (std::abs(distanceOnTrack-particle.GetLength()) < 5.*I3Units::mm) {
                        // do NOT correct the cascade time, it might be a delayed muon deacy (which should not be corrected)
                        log_debug("decaying muon detected, no timing correction for cascade at the track end.");
// for now, do NOT try to correct what we are given by MMC.
#ifdef TRY_TO_CORRECT_MMC
                    } else {
                        if (std::abs(expectedTime-daughter.GetTime()) > 2.*I3Units::ns) {
                            log_warn("Expected a cascade at time %fns (from its position on the track), but found it at t=%fns. Correcting.",
                                      expectedTime/I3Units::ns, daughter.GetTime()/I3Units::ns);

                            // correct the particle time
                            daughter.SetTime(expectedTime);
                        }

                        {
                            // if the cascade is at the very beginning or end of the track
                            // make sure it gets the exact same time as the track
                            double newTime=daughter.GetTime();
                            if ((newTime < ti) && (newTime > ti-0.1*I3Units::ns)) newTime=ti;
                            if ((newTime > tf) && (newTime < tf+0.1*I3Units::ns)) newTime=tf;
                            daughter.SetTime(newTime);
                        }

                        if ((daughter.GetTime() < ti) || (daughter.GetTime() > tf)) {
                            log_error("skipped a cascade that is not within the muon track time bounds!, muon_time_range=[%f,%f]ns cascade_time=%fns. LIGHT IS LOST!",
                                      ti/I3Units::ns, tf/I3Units::ns, daughter.GetTime());
                            continue;
                        }
#endif
                    }
                    
                    double sliceDuration = expectedTime-currentTime;
                    if (sliceDuration<0.) sliceDuration=0.;
                    
                    const double sliceLength = sliceDuration*I3Constants::c;
                    const double lengthFromVertexToSliceStart = (currentTime-particle.GetTime())*I3Constants::c;
                    
                    I3Particle muonSlice;
                    muonSlice.SetDir(particle.GetDir());
                    muonSlice.SetPos(particle.GetPos().GetX() + particle.GetDir().GetX() * lengthFromVertexToSliceStart,
                                     particle.GetPos().GetY() + particle.GetDir().GetY() * lengthFromVertexToSliceStart,
                                     particle.GetPos().GetZ() + particle.GetDir().GetZ() * lengthFromVertexToSliceStart);
                    muonSlice.SetTime(currentTime);
                    muonSlice.SetLength(sliceLength);
                    muonSlice.SetEnergy(currentEnergy);
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
                    muonSlice.SetPdgEncoding(particle.GetPdgEncoding());
#else
                    muonSlice.SetType(particle.GetType());
#endif
                    muonSlice.SetShape(particle.GetShape());
                    muonSlice.SetFitStatus(I3Particle::NotSet);
                    muonSlice.SetLocationType(particle.GetLocationType());
                
                    if (sliceLength>1.*I3Units::km)
                        log_warn("Extremely long slice detected in iteration %u len=%fkm (hadInvalidEi=%s, hadInvalidEf=%s)",
                                 iterationNum,
                                 sliceLength/I3Units::km,
                                 hadInvalidEi?"YES":"NO",
                                 hadInvalidEf?"YES":"NO");
                    
                    if (sliceLength>=0.1*I3Units::mm) {
                        outputTree.append_child(particle_it_outputTree, muonSlice);
                    }

                    I3MCTree::iterator daughter_it_outputTree =
                    outputTree.append_child(particle_it_outputTree, daughter);
                    
                    SliceMuonOrCopySubtree(inputTree,
                       mmcTrackList,
                       mmcTrackListIndex,
                       outputTree,
                       daughter_it,
                       daughter_it_outputTree);

                    currentTime+=sliceDuration;
                    currentEnergy-=daughter.GetEnergy()+dEdt*sliceDuration;
                    
                    ++iterationNum;
                }

                if (iterationNum==0)
                {
                    // if it had no daughters and Ei and Ef are invalid, it
                    // seems to be an outgoing muon behind the can. Ignore it.
                    if ((hadInvalidEi) && (hadInvalidEf)) {
                        log_debug("Ignored an outgoing muon that starts behind the detector.");
                        return;
                    }
                }
                
                if (currentEnergy < 0.) {
                    log_debug("muon decayed. (currentEnergy=%gGeV)", currentEnergy/I3Units::GeV);
                    return;
                }
                
                // slice after last daughter
                double sliceDuration = tf-currentTime;
                const double sliceLength = sliceDuration*I3Constants::c;
                const double lengthFromVertexToSliceStart = (currentTime-particle.GetTime())*I3Constants::c;

                if (sliceLength <= 0.1*I3Units::mm) {
                    log_debug("ignoring muon slice with negative length. L=%fm", sliceLength/I3Units::m);
                    return; // do not write tracks with no (or negative) length
                }

                I3Particle muonSlice;
                muonSlice.SetDir(particle.GetDir());
                muonSlice.SetPos(particle.GetPos().GetX() + particle.GetDir().GetX() * lengthFromVertexToSliceStart,
                                 particle.GetPos().GetY() + particle.GetDir().GetY() * lengthFromVertexToSliceStart,
                                 particle.GetPos().GetZ() + particle.GetDir().GetZ() * lengthFromVertexToSliceStart);
                muonSlice.SetTime(currentTime);
                muonSlice.SetLength(sliceLength);
                muonSlice.SetEnergy(currentEnergy);
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
                muonSlice.SetPdgEncoding(particle.GetPdgEncoding());
#else
                muonSlice.SetType(particle.GetType());
#endif
                muonSlice.SetShape(particle.GetShape());
                muonSlice.SetFitStatus(I3Particle::NotSet);
                muonSlice.SetLocationType(particle.GetLocationType());
                
                if (sliceLength>1.*I3Units::km)
                    log_warn("Extremely long slice detected in iteration %u len=%fkm (hadInvalidEi=%s, hadInvalidEf=%s)",
                             iterationNum,
                             sliceLength/I3Units::km,
                             hadInvalidEi?"YES":"NO",
                             hadInvalidEf?"YES":"NO");

                outputTree.append_child(particle_it_outputTree, muonSlice);

                return; // return here in order _not_ to do the default loop below
            }
            
        }

        // It's either something else or a muon without a length.
        // Add all the daughters and recurse.
        BOOST_FOREACH(const I3MCTree::iterator &daughter_it_inputTree, daughterIterators)
        {
            // add the particle to the output tree and get an iterator
            I3MCTree::iterator daughter_it_outputTree =
            outputTree.append_child(particle_it_outputTree, *daughter_it_inputTree);
            
            SliceMuonOrCopySubtree(inputTree,
                                   mmcTrackList,
                                   mmcTrackListIndex,
                                   outputTree,
                                   daughter_it_inputTree,
                                   daughter_it_outputTree);
        }
        
    }
    
    
}

void I3MuonSlicer::DAQ(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3MCTreeConstPtr inputMCTree = frame->Get<I3MCTreeConstPtr>(inputMCTreeName_);
    if (!inputMCTree) {
        log_debug("Frame does not contain an I3MCTree named \"%s\".",
                  inputMCTreeName_.c_str());
        PushFrame(frame);
        return;
    }

    I3MMCTrackListConstPtr MMCTrackList = frame->Get<I3MMCTrackListConstPtr>(MMCTrackListName_);
    if (!MMCTrackList) {
        log_debug("Frame does not contain an I3MMCTrackList named \"%s\". Using an empty list for this frame.",
                                 MMCTrackListName_.c_str());

        // make our own empty MMCTrackList
        MMCTrackList = I3MMCTrackListConstPtr(new I3MMCTrackList);
    }

    // build an index into the MMCTrackList (by particle ID)
    std::map<std::pair<uint64_t, int>, const I3MMCTrack *> MMCTrackListIndex;
    BOOST_FOREACH(const I3MMCTrack &mmcTrack, *MMCTrackList)
    {
        const std::pair<uint64_t, int> identifier(mmcTrack.GetI3Particle().GetMajorID(),
                                                  mmcTrack.GetI3Particle().GetMinorID());
        
        if (!MMCTrackListIndex.insert(std::make_pair(identifier, &mmcTrack)).second) {
            log_fatal("Particle ID exists more than once in I3MMCTrackList.");
        }
    }
    
    // allocate the output I3MCTree
    I3MCTreePtr outputMCTree(new I3MCTree());
    
    // get a list of primaries
    const std::vector<I3Particle> primaries = I3MCTreeUtils::GetPrimaries(inputMCTree);
    
    // add each one to the output tree and check their children
    BOOST_FOREACH(const I3Particle &primary, primaries)
    {
        if ((primary.GetShape() != I3Particle::Primary) && (primary.GetShape() != I3Particle::Null) && (primary.GetShape() != I3Particle::Dark))
            log_warn("Input tree contains a particle with shape!=(Primary or Null or Dark) at its root. (shape=%s, type=%s)",
                      primary.GetShapeString().c_str(), primary.GetTypeString().c_str());
        
        // have to use I3TreeUtils here, I3MCTreeUtils::AddPrimary expects a non-const I3Particle..
        //I3MCTreeUtils::AddPrimary(outputMCTree, primary);
        I3TreeUtils::AddTopLevel<I3Particle>(*outputMCTree, primary); 

        I3MCTree::iterator primary_in_input_tree  = GetMCTreeIterator(*inputMCTree,  primary);
        I3MCTree::iterator primary_in_output_tree = GetMCTreeIterator(*outputMCTree, primary);

        SliceMuonOrCopySubtree(*inputMCTree,
                               *MMCTrackList,
                               MMCTrackListIndex,
                               *outputMCTree,
                               primary_in_input_tree,
                               primary_in_output_tree);
    }
    
    
    
    // store the output I3MCTree
    if ((outputMCTreeName_=="") || (outputMCTreeName_==inputMCTreeName_)) {
        frame->Delete(inputMCTreeName_);
        frame->Put(inputMCTreeName_, outputMCTree);
    } else {
        frame->Put(outputMCTreeName_, outputMCTree);
    }
    
    // that's it!
    PushFrame(frame);
}
