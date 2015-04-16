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
 * @file I3PhotonToMCHitConverterForMultiPMT.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#define __STDC_FORMAT_MACROS 
#include <inttypes.h>

// Things from this module:
#include "clsim/dom/I3PhotonToMCHitConverterForMultiPMT.h"

#include "clsim/I3Photon.h"
#include "dataclasses/physics/I3MCHit.h"
#include "dataclasses/physics/I3MCTree.h"

// IceTray things:
#include "icetray/I3Logging.h"
#include "dataclasses/I3Units.h"
#include "dataclasses/I3Constants.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/geometry/I3Geometry.h"
#include "phys-services/I3RandomService.h"

#include <algorithm>

using namespace std;

I3_MODULE(I3PhotonToMCHitConverterForMultiPMT);

/************
 Constructor
 *************/
I3PhotonToMCHitConverterForMultiPMT::I3PhotonToMCHitConverterForMultiPMT(const I3Context& ctx) 
: I3ConditionalModule(ctx)
{
    AddOutBox("OutBox");
    
    log_debug("Enter I3PhotonToMCHitConverterForMultiPMT::I3PhotonToMCHitConverterForMultiPMT()");
    
    // Add all the parameters that can be defined in the control script:
    AddParameter("RandomService",
                 "A random number generating service (derived from I3RandomService).",
                 randomService_);
    
    outputMCHitSeriesMapName_="MCHitSeriesMultiOMMap";
    AddParameter("OutputMultiOMMCHitMapName",
                 "Name of the output I3MCHitSeriesMultiOMMap frame object",
                 outputMCHitSeriesMapName_);
    
    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object",
                 inputPhotonSeriesMapName_);
    
    MCTreeName_="I3MCTree";
    AddParameter("MCTreeName",
                 "Name of the I3MCTree frame object. All photon particle IDs are checked against this tree.",
                 MCTreeName_);
    
}

/**********
 Configure
 ***********/
void I3PhotonToMCHitConverterForMultiPMT::Configure() 
{
    log_trace("Entering Configure()");
    
    // Get the parameters from the control script:
    GetParameter("RandomService", randomService_);
    
    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputMultiOMMCHitMapName", outputMCHitSeriesMapName_);
    
    GetParameter("MCTreeName", MCTreeName_);
    
    if (!randomService_) {
        log_info("No random service provided as a parameter, trying to get one from the context..");
        randomService_ = context_.Get<I3RandomServicePtr>();
    }
    
    if (!randomService_)
        log_fatal("No random service provided using the \"RandomService\" parameter (and there is none installed on the context).");
    
}

namespace {
    inline int FindHitPMT(const I3Position &photonPos, 
                          const I3Position &omPos, 
                          const I3Direction &photonDir, 
                          const I3Orientation &omOrientation,
                          const I3OMTypeInfo &om_typeinfo,
                          double &pathLengthInOM,
                          I3Direction &rotatedPmtDir)
    {
        const double omRadius = om_typeinfo.GetSphereDiameter()*0.5;    // get the OM's outer diameter from the geometry
        const double omRadiusSquared = omRadius*omRadius;
        
        const double px=photonPos.GetX()-omPos.GetX();
        const double py=photonPos.GetY()-omPos.GetY();
        const double pz=photonPos.GetZ()-omPos.GetZ();
        const double pr2 = px*px + py*py + pz*pz;
        
        const double dx = photonDir.GetX();
        const double dy = photonDir.GetY();
        const double dz = photonDir.GetZ();
        
        // is photon entering?
        const double dot = px*dx + py*dy + pz*dz;
        if (dot > 0.) {
            log_debug("photon is leaving, dot=%f", dot);
            return -1;
        }
        
        // sanity check: are photons on the OM's surface?
        const double distFromDOMCenter = std::sqrt(pr2);
        if (std::abs(distFromDOMCenter - omRadius) > 3.*I3Units::cm) {
            log_warn("distance not %fmm.. it is %fmm (diff=%gmm)",
                     omRadius/I3Units::mm,
                     distFromDOMCenter/I3Units::mm,
                     (distFromDOMCenter-omRadius)/I3Units::mm);
        }
        
        pathLengthInOM = NAN;
        
        
        log_trace("OM orientation=(%f,%f,%f)", omOrientation.GetX(), omOrientation.GetY(), omOrientation.GetZ());
        
        int foundIntersection=-1;
        unsigned int numPMTs = om_typeinfo.GetNumPMTs();
        for (unsigned int pmtNum=0; pmtNum<numPMTs; ++pmtNum)
        {
            const I3PMTInfo &pmt_info = om_typeinfo.GetPMTInfo(pmtNum);
            const double pmtRadius = pmt_info.GetDiameter()/2.;
            const double pmtRadiusSquared = pmtRadius*pmtRadius;
            
            double nx = pmt_info.GetDirection().GetX();
            double ny = pmt_info.GetDirection().GetY();
            double nz = pmt_info.GetDirection().GetZ();
            
            log_trace(" PMT %u dir before = (%f,%f,%f)", pmtNum, nx, ny, nz);
            
            // rotate this vector into the OM's orientation
            omOrientation.RotVectorInPlace(nx,ny,nz);
            
            log_trace(" PMT %u dir after  = (%f,%f,%f)", pmtNum, nx, ny, nz);
            
            const double nl = nx*nx + ny*ny + nz*nz;
            if (fabs(nl - 1.) > 1e-6) log_fatal("INTERNAL ERROR: rotation does change vector length!");
            
            // find the intersection of the PMT's surface plane and the photon's path
            const double denom = dx*nx + dy*ny + dz*nz; // should be < 0., test that:
            
            if (denom>=1e-8) continue; // no intersection, photon is moving towards the PMT's back
            
            double ax = pmt_info.GetPosition().GetX();
            double ay = pmt_info.GetPosition().GetY();
            double az = pmt_info.GetPosition().GetZ();
            const double al_before = ax*ax + ay*ay + az*az;
            
            // rotate this vector into the OM's orientation
            omOrientation.RotVectorInPlace(ax,ay,az);
            
            // sanity check
            const double al_after = ax*ax + ay*ay + az*az;
            if (fabs(al_before-al_after) > 1e-6) log_fatal("INTERNAL ERROR: rotation does change vector length!");
            
            if (omRadiusSquared < al_after) log_fatal("OM sphere radius too small for this PMT! You will never get hits! Seems to be an error in your geometry definition!");
            
            const double mu = ((ax-px)*nx + (ay-py)*ny + (az-pz)*nz)/denom;
            
            if (mu < 0.) continue; // no intersection, photon is moving away from PMT
            
            // calculate the distance of the point of intersection
            // from the PMT position:
            const double distFromPMTCenterSquared = 
            (ax-px-mu*dx)*(ax-px-mu*dx) + 
            (ay-py-mu*dy)*(ay-py-mu*dy) + 
            (az-pz-mu*dz)*(az-pz-mu*dz);
            
            if (distFromPMTCenterSquared > pmtRadiusSquared) continue; // photon outside the PMT radius
            
            // there is an intersection with a pmt!
            if (foundIntersection >= 0) {
                log_warn("found another intersection! previousPMT=#%u, thisPMT=#%u", foundIntersection, pmtNum);
                if ((isnan(pathLengthInOM)) || (mu < pathLengthInOM))
                {
                    log_warn(" -> new intersection is closer than previous one. using it.");
                }
                else
                {
                    log_warn(" -> new intersection is further away than previous one. keeping old one.");
                    continue;
                }
                
            }
            
            foundIntersection = pmtNum;
            pathLengthInOM = mu;
            rotatedPmtDir.SetDir(nx, ny, nz);
        }
        
        return foundIntersection;
    }
    
    
}

namespace {
    // Return whether first element is greater than the second
    bool MCHitTimeLess(const I3MCHit &elem1, const I3MCHit &elem2)
    {
        return elem1.GetTime() < elem2.GetTime();
    }
    
}

/********
 Physics
 *********/
void I3PhotonToMCHitConverterForMultiPMT::DAQ(I3FramePtr frame)
{
    log_trace("Entering Physics()");
    
    // First we need to get our geometry
    const I3Geometry& geometry = frame->Get<I3Geometry>();
    
    // retrieve the MC track
    I3PhotonSeriesMapConstPtr input_hitmap = frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if (!input_hitmap) log_fatal("I3PhotonSeriesMap object with name \"%s\" not found in the frame!", inputPhotonSeriesMapName_.c_str());
    
    // create a new output hit map
    I3MCHitSeriesMultiOMMapPtr output_hitmap = I3MCHitSeriesMultiOMMapPtr(new I3MCHitSeriesMultiOMMap);
    
    // currently, the only reason we need the MCTree is that I3MCHit does
    // only allow setting the major/minor particle IDs using an existing
    // I3Particle instance with that ID combination.
    I3MCTreeConstPtr MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);
    if (!MCTree) log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                           MCTreeName_.c_str());
    
    // build an index into the I3MCTree
    std::map<std::pair<uint64_t, int>, const I3Particle *> mcTreeIndex;
    for (I3MCTree::iterator it = MCTree->begin();
         it != MCTree->end(); ++it)
    {
        const I3Particle &particle = *it;
        mcTreeIndex.insert(std::make_pair(std::make_pair(particle.GetMajorID(), particle.GetMinorID()), &particle));
    }
    
    
    // track the first and last non-noise hit times
    for (I3PhotonSeriesMap::const_iterator om_it = input_hitmap->begin(); om_it != input_hitmap->end(); ++om_it)
    {
        OMKey key = om_it->first;
        
        // Find the current OM in the geometry map
        I3OMGeoMap::const_iterator geo_it = geometry.omgeo.find(key);
        if (geo_it == geometry.omgeo.end())
            log_fatal("OM (%i/%u) not found in the current geometry map!", key.GetString(), key.GetOM());
        
        // check if any type information exists for this OM and retrieve it
        if (!geometry.ExistsOMTypeInfo(key))
            log_fatal("No type information found for OM (%i/%u)!", key.GetString(), key.GetOM());
        const I3OMTypeInfo &om_typeinfo = geometry.GetOMTypeInfo(key);
        
        // create an entry in the output hitmap
        I3MCHitSeriesMultiOMMap::iterator output_hitmap_it = (output_hitmap->insert( std::make_pair(key, I3MCHitSeriesMultiOM()) )).first;
        I3MCHitSeriesMultiOM &multiOMMap = output_hitmap_it->second;
        
        // loop over all hits on this OM and construct the I3MCHit series
        for (I3PhotonSeries::const_iterator hit_it = om_it->second.begin(); hit_it != om_it->second.end(); ++hit_it)
        {
            const I3Photon &photon = *hit_it;
            
            double pathLengthInsideOM=NAN;
            I3Direction rotatedPmtDir;
            int hitPmtNum = FindHitPMT(photon.GetPos(),
                                       geo_it->second.position,
                                       photon.GetDir(), 
                                       geo_it->second.orientation,
                                       om_typeinfo,
                                       pathLengthInsideOM,
                                       rotatedPmtDir);
            if (hitPmtNum < 0) continue; // no PMT hit
            
            // Calculate the time the photon travels inside the OM.
            // This assumes that the photon does only travel within
            // the OM's glass. The Gel/Air is not taken into account.
            // TODO: improve that, include gel/air.
            const double travelTimeInOM = pathLengthInsideOM / (I3Constants::c / om_typeinfo.GetGlassRefIndex());
            log_trace("travelTimeInOM=%fns, pathLengthInsideOM=%fmm", travelTimeInOM/I3Units::ns, pathLengthInsideOM/I3Units::mm);
            
            const I3PMTInfo &pmt_info = om_typeinfo[hitPmtNum];
            
            const double hit_angle = acos(- (rotatedPmtDir.GetX()*photon.GetDir().GetX() +
                                             rotatedPmtDir.GetY()*photon.GetDir().GetY() +
                                             rotatedPmtDir.GetZ()*photon.GetDir().GetZ()));
            
            if (hit_angle >= 90.*I3Units::deg) continue; // flat disc cannot be hit from behind
            
            // TODO: FIXME: pathLengthInsideOM should be only the length within the glass,
            // not the full length. The rest should be a length in the gel, which is assumed
            // to be fixed right now..
            const double glassGelSurvival_fac = om_typeinfo.GetGlassGelSurvivalProbability(photon.GetWavelength(), pathLengthInsideOM);
            //const double glassGelSurvival_fac = om_typeinfo.GetGlassGelSurvivalProbability(photon.GetWavelength());
            const double qe_fac = pmt_info.GetQuantumEfficiency(photon.GetWavelength()) * pmt_info.GetCollectionEfficiency();
            
            // the angular acceptance is defined in the ANTARES sense: multiply the acceptance with
            // the geometrical area of the PMT and get the "effective" area for the current angle.
            // We already do the tracking to each PMT, so a geometrical factor of cos(theta) is
            // already there from the probability to hit a PMT. (The projected area of a PMT
            // scales with cos(theta).) But this geometrical scaling is ALSO included in the 
            // angular acceptance factor from the geometry. So we have to get rid of that first.
            // This means that after the code knows that a PMT is hit, the geometrical acceptance
            // should be 1 if the acceptance factor from the geometry is cos(theta).
            const double ang_fac = pmt_info.GetAngularAcceptanceFactor(hit_angle)/std::fabs(std::cos(hit_angle));
            
            // calculate the measurement probability
            double measurement_prob = photon.GetWeight();
            measurement_prob *= glassGelSurvival_fac;
            measurement_prob *= qe_fac;
            measurement_prob *= ang_fac;
            
            if (measurement_prob > 1.)
                log_fatal("measurement_prob > 1 (it's %f): cannot continue. your hit weights are too high. (weight=%f, wlen=%fnm)",
                          measurement_prob, photon.GetWeight(),
                          photon.GetWavelength()/I3Units::nanometer);
            
            if (measurement_prob <= randomService_->Uniform()) continue;
            
            // find the particle
            std::map<std::pair<uint64_t, int>, const I3Particle *>::const_iterator it = 
            mcTreeIndex.find(std::make_pair(photon.GetParticleMajorID(), photon.GetParticleMinorID()));
            if (it==mcTreeIndex.end())
                log_fatal("Particle with id maj=%" PRIu64 ", min=%i does not exist in MC tree, but we have a photon that claims it was created by that particle..",
                          photon.GetParticleMajorID(), photon.GetParticleMinorID());
            const I3Particle &particle = *(it->second);
            
            
            // get the hit series into which we are going to insert the hit
            I3MCHitSeries &hitSeries = multiOMMap.insert(std::make_pair((unsigned int)hitPmtNum, I3MCHitSeries())).first->second;
            
            // add a new hit
            hitSeries.push_back(I3MCHit());
            I3MCHit &hit = hitSeries.back();
            
            // fill in all information
            hit.SetTime(photon.GetTime());
            hit.SetHitID(photon.GetID());
#ifdef I3MCHIT_WEIGHT_IS_DEPRECATED
            hit.SetNPE(1);
#else
            hit.SetWeight(1.0);
#endif
            hit.SetParticleID(particle);
            hit.SetCherenkovDistance(NAN);
            hit.SetHitSource(I3MCHit::SPE); // SPE for now, no afterpulses yet..
        }
        
        // remove a hit vector from the output map if it is empty
        if (multiOMMap.size()<=0) {output_hitmap->erase(output_hitmap_it); continue;}
        
        // sort the hit vectors for each PMT on this OM
        for (I3MCHitSeriesMultiOM::iterator pmtIt = multiOMMap.begin();
             pmtIt != multiOMMap.end(); ++pmtIt)
        {
            I3MCHitSeries &hitSeries = pmtIt->second;
            
            // now sort by time, regardless of particle ID
            std::sort(hitSeries.begin(), hitSeries.end(), MCHitTimeLess);
        }
        
    }
    
    
    // put the outputhitmap into the frame
    frame->Put(outputMCHitSeriesMapName_, output_hitmap);
    
    // Push that frame:
    PushFrame(frame,"OutBox"); 
    
    log_trace("Leaving Physics()");
}



/**********
 Finishing 
 ***********/
void I3PhotonToMCHitConverterForMultiPMT::Finish() 
{
    log_trace("Finish()");
}



/***********
 Destructor  
 ************/
I3PhotonToMCHitConverterForMultiPMT::~I3PhotonToMCHitConverterForMultiPMT() 
{
    
}

