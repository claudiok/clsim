/**
 * Copyright (c) 2012
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
 * @file I3PhotonToMCHitConverterForMDOMs.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

// Things from this module:
#include "clsim/dom/I3PhotonToMCHitConverterForMDOMs.h"

#include "clsim/I3Photon.h"
#include "dataclasses/physics/I3MCHit.h"
#include "dataclasses/physics/I3MCTree.h"

// IceTray things:
#include "icetray/I3Logging.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Constants.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/physics/I3Particle.h"
#include "phys-services/I3RandomService.h"

#include "dataclasses/geometry/I3ModuleGeo.h"
#include "dataclasses/geometry/I3OMGeo.h"
#include "dataclasses/I3Map.h"

#include <set>
#include <algorithm>
#include <boost/foreach.hpp>

using namespace std;

I3_MODULE(I3PhotonToMCHitConverterForMDOMs);

/************
 Constructor
 *************/
I3PhotonToMCHitConverterForMDOMs::I3PhotonToMCHitConverterForMDOMs(const I3Context& ctx) 
: I3ConditionalModule(ctx)
{
    AddOutBox("OutBox");
    
    log_debug("Enter I3PhotonToMCHitConverterForMDOMs::I3PhotonToMCHitConverterForMDOMs()");
    
    // Add all the parameters that can be defined in the control script:
    AddParameter("RandomService",
                 "A random number generating service (derived from I3RandomService).",
                 randomService_);
    
    outputMCHitSeriesMapName_="MCHitSeriesMap";
    AddParameter("OutputMCHitSeriesMapName",
                 "Name of the output I3MCHitSeriesMap frame object",
                 outputMCHitSeriesMapName_);
    
    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object",
                 inputPhotonSeriesMapName_);
    
    MCTreeName_="I3MCTree";
    AddParameter("MCTreeName",
                 "Name of the I3MCTree frame object. All photon particle IDs are checked against this tree.",
                 MCTreeName_);

    AddParameter("IgnoreSubdetectors",
                 "A list of subdetectors to ignore",
                 ignoreSubdetectors_);

    AddParameter("PMTWavelengthAcceptance",
                 "Wavelength acceptance of a PMT as a I3CLSimFunction object.",
                 pmtWavelengthAcceptance_);
    
    AddParameter("PMTAngularAcceptance",
                 "Angular acceptance of a PMT as a I3CLSimFunction object.",
                 pmtAngularAcceptance_);

    AddParameter("GlassAbsorptionLength",
                 "The absorption length of the DOM pressure housing glass.",
                 glassAbsorptionLength_);
    
    glassThickness_=NAN;
    AddParameter("GlassThickness",
                 "The thickness of the DOM pressure housing glass (the module assumes that "
                 "the rest of the path from the DOM surface to the PMT is in optical gel.",
                 glassThickness_);

    AddParameter("GelAbsorptionLength",
                 "The absorption length of the optical gel between the DOM sphere and the PMT.",
                 gelAbsorptionLength_);

}

/**********
 Configure
 ***********/
void I3PhotonToMCHitConverterForMDOMs::Configure() 
{
    log_trace("Entering Configure()");
    
    // Get the parameters from the control script:
    GetParameter("RandomService", randomService_);
    
    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputMCHitSeriesMapName", outputMCHitSeriesMapName_);
    
    GetParameter("MCTreeName", MCTreeName_);
    GetParameter("IgnoreSubdetectors", ignoreSubdetectors_);

    GetParameter("PMTWavelengthAcceptance", pmtWavelengthAcceptance_);
    GetParameter("PMTAngularAcceptance", pmtAngularAcceptance_);

    GetParameter("GlassAbsorptionLength", glassAbsorptionLength_);
    GetParameter("GlassThickness", glassThickness_);
    GetParameter("GelAbsorptionLength", gelAbsorptionLength_);

    if (!pmtWavelengthAcceptance_)
        log_fatal("The \"PMTWavelengthAcceptance\" parameter must not be empty.");
    if (!pmtAngularAcceptance_)
        log_fatal("The \"PMTAngularAcceptance\" parameter must not be empty.");

    if (!glassAbsorptionLength_)
        log_fatal("The \"GlassAbsorptionLength\" parameter must not be empty.");
    if (isnan(glassThickness_))
        log_fatal("The \"GlassThickness\" parameter must not be empty.");

    if (!gelAbsorptionLength_)
        log_fatal("The \"GelAbsorptionLength\" parameter must not be empty.");

    if (!randomService_) {
        log_info("No random service provided as a parameter, trying to get one from the context..");
        randomService_ = context_.Get<I3RandomServicePtr>();
    }
    
    if (!randomService_)
        log_fatal("No random service provided using the \"RandomService\" parameter (and there is none installed on the context).");
    
}

namespace {
    inline int FindHitPMT(const I3Position &photonPos, 
                          const I3Direction &photonDir, 
                          const ModuleKey &key,
                          const I3ModuleGeo &moduleGeo,
                          const I3OMGeoMap &omGeoMap,
                          const std::vector<unsigned char> &pmtNumbersToCheck,
                          double glassThickness,
                          double &pathLengthInOM,
                          double &pathLengthInGlass)
    {
        const double omRadius = moduleGeo.GetRadius();    // get the OM's outer diameter from the geometry
        const double omRadiusSquared = omRadius*omRadius;
        
        // photon position relative to DOM position
        double px=photonPos.GetX()-moduleGeo.GetPos().GetX();
        double py=photonPos.GetY()-moduleGeo.GetPos().GetY();
        double pz=photonPos.GetZ()-moduleGeo.GetPos().GetZ();
        double pr2 = px*px + py*py + pz*pz;
        
        const double dx = photonDir.GetX();
        const double dy = photonDir.GetY();
        const double dz = photonDir.GetZ();
        
        {
            // sanity check: are photons on the OM's surface?
            const double distFromDOMCenter = std::sqrt(pr2);
            if (std::abs(std::sqrt(pr2) - omRadius) > 3.*I3Units::cm) {
                log_warn("distance not %fmm.. it is %fmm (diff=%gmm). correcting to DOM radius.",
                         omRadius/I3Units::mm,
                         std::sqrt(pr2)/I3Units::mm,
                         (std::sqrt(pr2)-omRadius)/I3Units::mm);
            }
            
            // to make sure that photons start *exactly* on the outside of the sphere
            const double pr_scale = omRadius/distFromDOMCenter;

            px *= pr_scale;
            py *= pr_scale;
            pz *= pr_scale;
            pr2 = omRadiusSquared;
        }
        
        // is photon entering?
        const double dot = px*dx + py*dy + pz*dz;
        if (dot > 0.) {
            log_warn("photon is leaving, dot=%f", dot);
            return -1;
        }

        pathLengthInOM = NAN;
        pathLengthInGlass = NAN;
        double distFromPMTCenterSquaredFound = NAN;
        double pmtRadiusSquaredFound = NAN;
        
        log_trace("OM orientation=(%f,%f,%f)",
                  moduleGeo.GetOrientation().GetX(),
                  moduleGeo.GetOrientation().GetY(),
                  moduleGeo.GetOrientation().GetZ());
        
        int foundIntersection=-1;
        BOOST_FOREACH(unsigned char pmtNum, pmtNumbersToCheck)
        {
            const OMKey pmtKey(key.GetString(), key.GetOM(), pmtNum);
            I3OMGeoMap::const_iterator pmt_it = omGeoMap.find(pmtKey);
            if (pmt_it == omGeoMap.end())
                log_fatal("Internal error. OMKey(%i,%u,%u) should exist in I3OMGeoMap.",
                          pmtKey.GetString(), pmtKey.GetOM(),
                          static_cast<unsigned int>(pmtKey.GetPMT()));
            const I3OMGeo &pmtInfo = pmt_it->second;
            
            const double pmtArea = pmtInfo.area;
            if (isnan(pmtArea)) log_fatal("OMKey(%i,%u,%u) has NaN area!",
                                          pmtKey.GetString(), pmtKey.GetOM(),
                                          static_cast<unsigned int>(pmtKey.GetPMT()));
            
            // assume a flat, disc-shaped PMT window
            const double pmtRadiusSquared = pmtArea/M_PI;
            log_trace("pmtRadius=%fmm, pmtArea=%fmm^2", std::sqrt(pmtRadiusSquared)/I3Units::mm, pmtArea/I3Units::mm2);
            
            // this is already in the final coordinate frame
            double nx = pmtInfo.orientation.GetX();
            double ny = pmtInfo.orientation.GetY();
            double nz = pmtInfo.orientation.GetZ();
            log_trace(" PMT %u dir  = (%f,%f,%f)", pmtNum, nx, ny, nz);
            
            const double nl = nx*nx + ny*ny + nz*nz;
            if (fabs(nl - 1.) > 1e-6) log_fatal("INTERNAL ERROR: rotation does change vector length!");
            
            // find the intersection of the PMT's surface plane and the photon's path
            const double denom = dx*nx + dy*ny + dz*nz; // should be < 0., test that:
            
            if (denom>=1e-8) continue; // no intersection, photon is moving towards the PMT's back
            
            const double ax = pmtInfo.position.GetX() - moduleGeo.GetPos().GetX();
            const double ay = pmtInfo.position.GetY() - moduleGeo.GetPos().GetY();
            const double az = pmtInfo.position.GetZ() - moduleGeo.GetPos().GetZ();

            //double ar2 = ax*ax + ay*ay + az*az;
            //log_trace("n*a/|a|=%f", (nx*ax + ny*ay + nz*az)/std::sqrt(ar2));

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
            
            foundIntersection = static_cast<int>(pmtNum);
            pathLengthInOM = mu;
            distFromPMTCenterSquaredFound = distFromPMTCenterSquared;
            pmtRadiusSquaredFound = pmtRadiusSquared;
        }

        if (foundIntersection >= 0)
        {
            // calculate the path length inside the glass
            
            const double R = omRadius-glassThickness;
            const double denom = dot*dot - omRadiusSquared + R*R;
            if (denom < 0.) 
                log_fatal("Path never enters the DOM. This could never have hit a PMT. Yet, it hit PMT %i, denom=%f, dot=%f, pathLengthInOM=%fmm",
                          foundIntersection, denom, dot,
                          pathLengthInOM/I3Units::mm);
            
            const double denom_s = std::sqrt(denom);
            
            const double l1 = -dot+denom_s;
            const double l2 = -dot-denom_s;
            const double l = std::min(l1,l2);
            if (l < 0.) log_fatal("Photon leaves DOM.");
            
            const double p2x = px + dx*l;
            const double p2y = py + dy*l;
            const double p2z = pz + dz*l;
            
            if (std::abs(p2x*p2x + p2y*p2y + p2z*p2z - R*R) > 1e-5)
                log_fatal("Internal error.");
            
            pathLengthInGlass = l;
            
            if (pathLengthInGlass - glassThickness < -0.1*I3Units::mm) {
                const double hx = px+pathLengthInOM*dx;
                const double hy = py+pathLengthInOM*dy;
                const double hz = pz+pathLengthInOM*dz;
                const double hr = std::sqrt(hx*hx + hy*hy + hz*hz);
                log_error("hr=%fmm, Rinner=%fmm, Router=%fmm, distFromPMTCenter=%fmm, pmtRadius=%fmm",
                         hr/I3Units::mm,
                         (omRadius-glassThickness)/I3Units::mm,
                         omRadius/I3Units::mm,
                         std::sqrt(distFromPMTCenterSquaredFound)/I3Units::mm,
                         std::sqrt(pmtRadiusSquaredFound)/I3Units::mm);

                log_fatal("Internal error: pathLengthInGlass=%fmm < glassThickness=%fmm, pathLengthInOM=%fmm",
                          pathLengthInGlass/I3Units::mm,
                          glassThickness/I3Units::mm,
                          pathLengthInOM/I3Units::mm);
            } else if (pathLengthInGlass < glassThickness) {
                pathLengthInGlass = glassThickness;
            }
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
void I3PhotonToMCHitConverterForMDOMs::DAQ(I3FramePtr frame)
{
    log_trace("Entering Physics()");
    
    // First we need to get our geometry
    I3ModuleGeoMapConstPtr moduleGeoMap = frame->Get<I3ModuleGeoMapConstPtr>("I3ModuleGeoMap");
    I3OMGeoMapConstPtr omGeoMap = frame->Get<I3OMGeoMapConstPtr>("I3OMGeoMap");
    
    if (!moduleGeoMap) log_fatal("Frame does not have I3ModuleGeoMap");
    if (!omGeoMap) log_fatal("Frame does not have I3OMGeoMap");

    std::set<ModuleKey> ignoreModules;
    if (!ignoreSubdetectors_.empty()) {
        I3MapModuleKeyStringConstPtr subdetectors = frame->Get<I3MapModuleKeyStringConstPtr>("Subdetectors");
        if (!subdetectors) log_fatal("You chose to ignore certain subdetectors, but your geometry is missing a \"Subdetectors\" object!");
        
        for (I3MapModuleKeyString::const_iterator it=subdetectors->begin();
             it != subdetectors->end(); ++it)
        {
            bool shouldBeIgnored=false;
            BOOST_FOREACH(const std::string &sd, ignoreSubdetectors_)
            {
                if (it->second == sd) {
                    shouldBeIgnored=true;
                    break;
                }
            }
            
            if (shouldBeIgnored) {
                ignoreModules.insert(it->first);
            }
        }
    }
    
    // build a map of ModuleKey to vector<PMTNumber> for lookup
    std::map<ModuleKey, std::vector<unsigned char> > PMTsInModule;
    
    for (I3OMGeoMap::const_iterator pmt_it = omGeoMap->begin(); pmt_it != omGeoMap->end(); ++pmt_it)
    {
        const OMKey &omKey = pmt_it->first;
        const ModuleKey moduleKey(omKey.GetString(), omKey.GetOM());
        
        if (ignoreModules.count(moduleKey) > 0) continue; // ignore masked DOMs
        
        if (moduleGeoMap->find(moduleKey) == moduleGeoMap->end())
            log_fatal("OMKey(%i,%u,%u) does not have a corresponding ModuleKey(%i,%u)",
                      omKey.GetString(), omKey.GetOM(), static_cast<unsigned int>(omKey.GetPMT()),
                      moduleKey.GetString(), moduleKey.GetOM());
        
        // insert new vector or retrieve existing vector
        std::vector<unsigned char> &pmtsInThisModule =
        PMTsInModule.insert(std::make_pair(moduleKey, std::vector<unsigned char>())).first->second;
        
        pmtsInThisModule.push_back(omKey.GetPMT());
    }

    // are there any ModuleKeys without PMTs?
    for (I3ModuleGeoMap::const_iterator module_it = moduleGeoMap->begin(); module_it != moduleGeoMap->end(); ++module_it)
    {
        const ModuleKey &key = module_it->first;
        if (ignoreModules.count(key) > 0) continue; // ignore masked DOMs

        std::map<ModuleKey, std::vector<unsigned char> >::iterator it = PMTsInModule.find(key);
        if (it == PMTsInModule.end()) {
            log_warn("Your module ModuleKey(%i,%u) does not have any PMTs!",
                     key.GetString(), key.GetOM());
        }

        // insert an empty vector
        PMTsInModule.insert(std::make_pair(key, std::vector<unsigned char>()));
    }
    
    
    // retrieve the MC track
    I3PhotonSeriesMapConstPtr input_hitmap = frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if (!input_hitmap) log_fatal("I3PhotonSeriesMap object with name \"%s\" not found in the frame!", inputPhotonSeriesMapName_.c_str());
    
    // create a new output hit map
    I3MCHitSeriesMapPtr output_hitmap = I3MCHitSeriesMapPtr(new I3MCHitSeriesMap);
    
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
    
    
    for (I3PhotonSeriesMap::const_iterator om_it = input_hitmap->begin(); om_it != input_hitmap->end(); ++om_it)
    {
        const ModuleKey &key = om_it->first;
        if (ignoreModules.count(key) > 0) continue; // ignore masked DOMs

        // Find the current OM in the geometry map
        I3ModuleGeoMap::const_iterator geo_it = moduleGeoMap->find(key);
        if (geo_it == moduleGeoMap->end())
            log_fatal("Module (%i/%u) not found in the current geometry map!", key.GetString(), key.GetOM());
        const I3ModuleGeo &moduleGeo = geo_it->second;
        
        std::map<ModuleKey, std::vector<unsigned char> >::const_iterator pmt_index_it = PMTsInModule.find(key);
        if (pmt_index_it == PMTsInModule.end())
            log_fatal("Internal error Module (%i/%u) not found in the current geometry map!", key.GetString(), key.GetOM());
        const std::vector<unsigned char> &checkPMTNumbers = pmt_index_it->second;
        
        // loop over all hits on this OM and construct the I3MCHit series
        for (I3PhotonSeries::const_iterator hit_it = om_it->second.begin(); hit_it != om_it->second.end(); ++hit_it)
        {
            const I3Photon &photon = *hit_it;
            
            double pathLengthInsideOM=NAN;
            double pathLengthInsideGlass=NAN;
            I3Direction rotatedPmtDir;
            int hitPmtNum = FindHitPMT(photon.GetPos(),
                                       photon.GetDir(), 
                                       key,
                                       moduleGeo,
                                       *omGeoMap,
                                       checkPMTNumbers,
                                       glassThickness_,
                                       pathLengthInsideOM,
                                       pathLengthInsideGlass);
            if (hitPmtNum < 0) continue; // no PMT hit
            
            const OMKey pmtKey(key.GetString(), key.GetOM(), static_cast<unsigned char>(hitPmtNum));
            
            I3OMGeoMap::const_iterator pmt_it = omGeoMap->find(pmtKey);
            if (pmt_it == omGeoMap->end()) log_fatal("Internal error. OMKey(%i,%u,%u) does not exist.",
                                                     pmtKey.GetString(), pmtKey.GetOM(),
                                                     static_cast<unsigned int>(pmtKey.GetPMT()));
            const I3OMGeo &pmtGeo = pmt_it->second;
            const I3Direction pmtDir = pmtGeo.GetDirection();

            const double hit_cosangle = - (pmtDir.GetX()*photon.GetDir().GetX() +
                                           pmtDir.GetY()*photon.GetDir().GetY() +
                                           pmtDir.GetZ()*photon.GetDir().GetZ());
            if (hit_cosangle <= 0.) continue; // the flat disc window cannot be hit from behind

            const double wlen = photon.GetWavelength();

            //log_warn("corrected glass thickness: %fmm",
            //         pathLengthInsideGlass/I3Units::mm);

            if (pathLengthInsideOM < pathLengthInsideGlass) {
                log_fatal("A PMT window is inside the DOM glass sphere. pathLengthInsideOM=%fmm, pathLengthInsideGlass=%fmm, OMKey(%i,%u,%u)",
                          pathLengthInsideOM/I3Units::mm,
                          pathLengthInsideGlass/I3Units::mm,
                          pmtKey.GetString(), pmtKey.GetOM(),
                          static_cast<unsigned int>(pmtKey.GetPMT()));
            }
            
            const double gelThickness = pathLengthInsideOM-glassThickness_;
            log_trace("gelThickness=%fmm", gelThickness/I3Units::mm);
            
            const double glassGelSurvival_fac =
                std::exp(-glassThickness_/glassAbsorptionLength_->GetValue(wlen)
                         -gelThickness/gelAbsorptionLength_->GetValue(wlen));
            
            const double qe_fac = pmtWavelengthAcceptance_->GetValue(wlen);
            
            // the angular acceptance is defined in the usual way: multiply the acceptance with
            // the geometrical area of the PMT and get the "effective" area for the current angle.
            // We already do the tracking to each PMT, so a geometrical factor of cos(theta) is
            // already there from the probability to hit a PMT. (The projected area of a PMT
            // scales with cos(theta).) But this geometrical scaling is ALSO included in the 
            // angular acceptance factor from the geometry. So we have to get rid of that first.
            // This means that after the code knows that a PMT is hit, the geometrical acceptance
            // should be 1 if the acceptance factor from the geometry is cos(theta).
            const double ang_fac = pmtAngularAcceptance_->GetValue(hit_cosangle)/std::fabs(hit_cosangle);
            
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
            I3MCHitSeries &hitSeries = output_hitmap->insert(std::make_pair(pmtKey, I3MCHitSeries())).first->second;
            
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
            hit.SetHitSource(I3MCHit::SPE); // SPE for now, no afterpulses at this point
        }
        
    }

    // sort the hit vectors for each PMT
    for (I3MCHitSeriesMap::iterator pmtIt = output_hitmap->begin();
         pmtIt != output_hitmap->end(); ++pmtIt)
    {
        I3MCHitSeries &hitSeries = pmtIt->second;
        
        // now sort by time, regardless of particle ID
        std::sort(hitSeries.begin(), hitSeries.end(), MCHitTimeLess);
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
void I3PhotonToMCHitConverterForMDOMs::Finish() 
{
    log_trace("Finish()");
}



/***********
 Destructor  
 ************/
I3PhotonToMCHitConverterForMDOMs::~I3PhotonToMCHitConverterForMDOMs() 
{
    
}

