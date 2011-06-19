/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3MuonSlicer.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 *
 *
 *  This file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *  
 */

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "clsim/I3MuonSlicer.h"

#include <boost/foreach.hpp>

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
    inline bool IsAnyOfType(const std::vector<I3Particle> &particles, I3Particle::ParticleType type)
    {
        BOOST_FOREACH(const I3Particle &particle, particles)
        {
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
    
    inline bool AreParticlesSortedInTime(const std::vector<I3Particle> &particles)
    {
        bool firstIt=true;
        double previousTime=NAN;
        
        BOOST_FOREACH(const I3Particle &particle, particles)
        {
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
    
    inline double GetTotalEnergyOfParticles(const std::vector<I3Particle> &particles,
                                            double tMin, double tMax)
    {
        double totalEnergy=0.;
        
        BOOST_FOREACH(const I3Particle &particle, particles)
        {
            if ((particle.GetTime() < tMin) || (particle.GetTime() > tMax))
                continue;

            totalEnergy+=particle.GetEnergy();
        }
        
        return totalEnergy;
    }
                                            
    
    void SliceMuonOrCopySubtree(const I3MCTree &inputTree,
                                const I3MMCTrackList &mmcTrackList,
                                const std::map<std::pair<uint64_t, int>, const I3MMCTrack *> &mmcTrackListIndex,
                                I3MCTree &outputTree,
                                const I3Particle &particle
                                )
    {
        const std::vector<I3Particle> daughters =
        I3MCTreeUtils::GetDaughters(inputTree, particle);

        // special treatment for muons with a length only
        if (((!isnan(particle.GetLength())) && (particle.GetLength() > 0.)) &&
            ((particle.GetType()==I3Particle::MuMinus) ||
            (particle.GetType()==I3Particle::MuPlus)))
        {
            // is any of the daughters a muon?
            if (IsAnyOfType(daughters, I3Particle::MuMinus) ||
                IsAnyOfType(daughters, I3Particle::MuPlus) ||
                IsAnyOfType(daughters, I3Particle::unknown))
            {
                log_fatal("It seems you either ran MMC with the \"-recc\" option or I3MuonSlicer has already been applied.");
            }

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
            if (!AreParticlesSortedInTime(daughters))
            {
                log_fatal("Muon daughters are not sorted in time (ascending).");
            }

            
            // correct values with information from I3Particle
            if ((Ei<=0.) || (isnan(Ei)))
            {
                Ei = particle.GetEnergy();
                ti = particle.GetTime();
            }
            if ((Ef<0.) || (isnan(Ef)))
            {
                Ef = 0.;
                tf = particle.GetTime() + particle.GetLength()/I3Constants::c;
            }

            if (isnan(ti)) log_fatal("t_initial is NaN");
            if (isnan(tf)) log_fatal("t_final is NaN");

            if (Ei<=0) 
            {
                if (daughters.size() > 0)
                    log_fatal("Muon with Energy==0 has children.");
            }
            else if (tf<ti)
            {
                log_warn("Muon stops before it starts.. (t_final==%fns < t_initial==%fns) (minorID=%i, majorID=%" PRIu64 ") ignoring.",
                          tf/I3Units::ns, ti/I3Units::ns, particle.GetMinorID(), particle.GetMajorID());
            }
            else if (tf==ti)
            {
                if (particle.GetLength()>0.)
                    log_warn("Particle has length but t_final==t_initial==%fns. (minorID=%i, majorID=%" PRIu64 ") ignoring.",
                              tf/I3Units::ns, particle.GetMinorID(), particle.GetMajorID());
            }
            else // (ti<tf) && (Ei>0.)
            {
                const double totalEnergyInCascades = 
                GetTotalEnergyOfParticles(daughters, ti, tf);
                const double dEdt_calc = ((Ef-Ei+totalEnergyInCascades)/(ti-tf));
                const double dEdt_max = (0.21+8.8e-3*log(Ei/I3Units::GeV)/log(10.))*(I3Units::GeV/I3Units::m)*I3Constants::c;  // for ice only at Ecut=500 MeV (stolen from PPC)
                const double dEdt = std::min(dEdt_calc,dEdt_max);
                
                // add all daughters to the muon (in the output tree)
                // while inserting the muon slices between them
                
                double currentEnergy = Ei;
                double currentTime = ti;
                
                BOOST_FOREACH(const I3Particle &daughter, daughters)
                {
                    if (isnan(daughter.GetTime())) continue;
                    if ((daughter.GetTime() < ti) || (daughter.GetTime() > tf))
                        continue;
                 
                    const double sliceDuration = daughter.GetTime()-currentTime;
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
                    muonSlice.SetType(particle.GetType());
                    muonSlice.SetShape(particle.GetShape());
                    muonSlice.SetFitStatus(I3Particle::NotSet);
                    muonSlice.SetLocationType(particle.GetLocationType());
                
                    I3MCTreeUtils::AppendChild(outputTree, particle, muonSlice);
                    I3MCTreeUtils::AppendChild(outputTree, particle, daughter);
                    
                    currentTime+=sliceDuration;
                    currentEnergy-=daughter.GetEnergy()+dEdt*sliceDuration;
                    
                    if (currentEnergy<0.)
                        log_fatal("Muon looses more energy than it has.");
                }

                // slice after last daughter
                const double sliceDuration = tf-currentTime;
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
                muonSlice.SetType(particle.GetType());
                muonSlice.SetShape(particle.GetShape());
                muonSlice.SetFitStatus(I3Particle::NotSet);
                muonSlice.SetLocationType(particle.GetLocationType());
                
                I3MCTreeUtils::AppendChild(outputTree, particle, muonSlice);
                
                return; // return here in order _not_ to do the default loop below
            }
            
        }

        // It's either something else or a muon without a length.
        // Add all the daughters and recurse.
        BOOST_FOREACH(const I3Particle &daughter, daughters)
        {
            I3MCTreeUtils::AppendChild(outputTree, particle, daughter);
            SliceMuonOrCopySubtree(inputTree,
                                   mmcTrackList,
                                   mmcTrackListIndex,
                                   outputTree,
                                   daughter);
        }
        
    }
    
    
}

#ifdef IS_Q_FRAME_ENABLED
void I3MuonSlicer::DAQ(I3FramePtr frame)
#else
void I3MuonSlicer::Physics(I3FramePtr frame)
#endif
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3MCTreeConstPtr inputMCTree = frame->Get<I3MCTreeConstPtr>(inputMCTreeName_);
    if (!inputMCTree) log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                                inputMCTreeName_.c_str());

    I3MMCTrackListConstPtr MMCTrackList = frame->Get<I3MMCTrackListConstPtr>(MMCTrackListName_);
    if (!MMCTrackList) log_fatal("Frame does not contain an I3MMCTrackList named \"%s\".",
                                 MMCTrackListName_.c_str());

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
        if ((primary.GetShape() != I3Particle::Primary) && (primary.GetShape() != I3Particle::Null))
            log_fatal("Input tree contains a particle with shape!=(Primary or Null) at its root. (shape=%s, type=%s)",
                      primary.GetShapeString().c_str(), primary.GetTypeString().c_str());
        
        // have to use I3TreeUtils here, I3MCTreeUtils::AddPrimary expects a non-const I3Particle..
        //I3MCTreeUtils::AddPrimary(outputMCTree, primary);
        I3TreeUtils::AddTopLevel<I3Particle>(*outputMCTree, primary); 

        SliceMuonOrCopySubtree(*inputMCTree,
                               *MMCTrackList,
                               MMCTrackListIndex,
                               *outputMCTree,
                               primary);
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
