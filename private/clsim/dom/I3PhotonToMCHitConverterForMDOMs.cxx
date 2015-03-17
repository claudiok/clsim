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
#include "simclasses/I3MCPE.h"
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
    
    outputMCHitSeriesMapName_="MCPEHitSeriesMap";
    AddParameter("OutputMCPESeriesMapName",
                 "Name of the output I3MCPESeriesMap frame object",
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

    glassThickness_=1.4*I3Units::cm;
    AddParameter("GlassThickness",
                 "The thickness of the DOM pressure housing glass (the module assumes that "
                 "the rest of the path from the DOM surface to the PMT is in optical gel.",
                 glassThickness_);

    DOMOversizeFactor_=1.;
    AddParameter("DOMOversizeFactor",
                 "Specifiy the \"oversize factor\" (i.e. DOM radius scaling factor) you used during the CLSim run.\n"
                 "The photon arrival times will be corrected. In practice this means your large spherical DOMs will\n"
                 "become ellipsoids.",
                 DOMOversizeFactor_);

    DOMPancakeFactor_=1.;
    AddParameter("DOMPancakeFactor",
                 "Specifiy the \"pancake factor\" of a DOM. This is the factor a DOM has been *shrunk* again\n"
                 "(in the direction of the photon) after oversizing. You should set this to whatever\n"
                 "value you used during running I3CLSimModule. And most of the time this is the same as the\n"
                 "oversize factor.",
                 DOMPancakeFactor_);

#if 0
    AddParameter("GlassAbsorptionLength",
                 "The absorption length of the DOM pressure housing glass.",
                 glassAbsorptionLength_);
    


    AddParameter("GelAbsorptionLength",
                 "The absorption length of the optical gel between the DOM sphere and the PMT.",
                 gelAbsorptionLength_);
#endif
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
    GetParameter("OutputMCPESeriesMapName", outputMCHitSeriesMapName_);
    
    GetParameter("MCTreeName", MCTreeName_);
    GetParameter("IgnoreSubdetectors", ignoreSubdetectors_);

    GetParameter("PMTWavelengthAcceptance", pmtWavelengthAcceptance_);
    GetParameter("PMTAngularAcceptance", pmtAngularAcceptance_);

    GetParameter("GlassThickness", glassThickness_);

    GetParameter("DOMOversizeFactor", DOMOversizeFactor_);
    GetParameter("DOMPancakeFactor", DOMPancakeFactor_);
    
#if 0
    GetParameter("GlassAbsorptionLength", glassAbsorptionLength_);
    GetParameter("GelAbsorptionLength", gelAbsorptionLength_);
#endif

    if (!pmtWavelengthAcceptance_)
        log_fatal("The \"PMTWavelengthAcceptance\" parameter must not be empty.");
    if (!pmtAngularAcceptance_)
        log_fatal("The \"PMTAngularAcceptance\" parameter must not be empty.");

// NB: assume that pmtWavelengthAcceptance_ is a total efficiency calibration
#if 0
    if (!glassAbsorptionLength_)
        log_fatal("The \"GlassAbsorptionLength\" parameter must not be empty.");
    if (isnan(glassThickness_))
        log_fatal("The \"GlassThickness\" parameter must not be empty.");

    if (!gelAbsorptionLength_)
        log_fatal("The \"GelAbsorptionLength\" parameter must not be empty.");
#endif
    if (!randomService_) {
        log_info("No random service provided as a parameter, trying to get one from the context..");
        randomService_ = context_.Get<I3RandomServicePtr>();
    }
    
    if (!randomService_)
        log_fatal("No random service provided using the \"RandomService\" parameter (and there is none installed on the context).");
    
}

namespace {
    inline int FindHitPMT(const I3Photon &photon, 
                          const ModuleKey &key,
                          const I3ModuleGeo &moduleGeo,
                          const I3OMGeoMap &omGeoMap,
                          const std::vector<unsigned char> &pmtNumbersToCheck,
                          double glassThickness,
                          double DOMOversizeFactor,
                          double DOMPancakeFactor,
                          double &hitTime,
                          double &pathLengthInOM,
                          double &pathLengthInGlass)
    {
        const double omRadius = moduleGeo.GetRadius();    // get the OM's outer diameter from the geometry
        const double omRadiusSquared = omRadius*omRadius;
        
        // photon position relative to DOM position
        I3Position p(photon.GetPos() - moduleGeo.GetPos());
        const I3Direction &d = photon.GetDir();
        
        if (DOMOversizeFactor == 1.) {
            // sanity check: are photons on the OM's surface?
            const double distFromDOMCenter = p.Magnitude();
            if (std::abs(p.Magnitude() - omRadius) > 3.*I3Units::cm) {
                log_warn("distance not %fmm.. it is %fmm (diff=%gmm). correcting to DOM radius.",
                         omRadius/I3Units::mm,
                         p.Magnitude()/I3Units::mm,
                         (p.Magnitude()-omRadius)/I3Units::mm);
            }
            
            // to make sure that photons start *exactly* on the outside of the sphere
            p *= omRadius/distFromDOMCenter;
            
        } else {
            // The photon was recorded on the surface of an oversized, oblate
            // spheroid. Pull it to the surface of the physical OM sphere.

            // The orientation of the photon w.r.t. the OM. "Direction" points
            // along the photon direction, and "Right" away from the DOM center.
            // "Up" is a counterclockwise rotation around the OM.
            assert(d*p <= 0);
            const I3Orientation axis(d, I3Direction(-d.Cross(p)));
            const I3Position op(p);

            // Pull perpendicularly towards the module's center
            double odir = p*axis.GetDir();
            p -= (1.-1./DOMOversizeFactor)*(p*axis.GetRight())
                *axis.GetRight();
            // ensure that the photon is now in the cross-sectional area of the
            // physical OM
            assert(p*axis.GetRight() >= 0);
            assert(p*axis.GetRight() <= omRadius); 
            double ndir = p*axis.GetDir();
            assert(abs(odir-ndir) < I3Units::mm);
            
            // Pull parallelly towards the module's center
            double longStep = (1.-DOMPancakeFactor/DOMOversizeFactor)*(p*axis.GetDir());
            const I3Position opp(p);
            p -= longStep*axis.GetDir();
            // ensure that the photon is now on the surface of the physical OM
            assert(std::abs(p.Magnitude() - omRadius) < I3Units::mm);
            
            // Correct timing
            hitTime += longStep/photon.GetGroupVelocity();

        }
        
        // is photon entering?
        const double dot = p*d;
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
        int intersections = 0;
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
            I3Direction n(pmtInfo.orientation.GetDir());
            log_trace(" PMT %u dir  = (%f,%f,%f)", pmtNum, n.GetX(), n.GetY(), n.GetZ());
            
            // find the intersection of the PMT's surface plane and the photon's path
            const double denom = d*n; // should be < 0., test that:
            if (denom>=1e-8) continue; // no intersection, photon is moving towards the PMT's back
            
            const I3Position a(pmtInfo.position - moduleGeo.GetPos());
            assert(a.Magnitude() > 4*I3Units::cm);

            //double ar2 = ax*ax + ay*ay + az*az;
            //log_trace("n*a/|a|=%f", (nx*ax + ny*ay + nz*az)/std::sqrt(ar2));

            const double mu = (a-p)*n/denom;
            assert(std::isfinite(mu));
            
            if (mu < 0.) continue; // no intersection, photon is moving away from PMT
            
            // calculate the distance of the point of intersection
            // from the PMT position:
            const double distFromPMTCenterSquared = (a-p-mu*d).Mag2(); 
            assert(std::isfinite(distFromPMTCenterSquared));
            
            if (distFromPMTCenterSquared > pmtRadiusSquared) continue; // photon outside the PMT radius

            log_debug_stream("hit PMT#" << int(pmtNum) << " after " << (mu/I3Units::cm) << " cm");

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
            intersections++;
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
            
            const I3Position p2 = p + l*d;
            
            // Ensure that p2 is on the inner surface of the housing
            assert(std::abs(p2.Magnitude() - R) < I3Units::mm);
            
            pathLengthInGlass = l;
            
            if (pathLengthInGlass - glassThickness < -0.1*I3Units::mm) {
                const I3Position h = p + pathLengthInOM*d;
                log_error("hr=%fmm, Rinner=%fmm, Router=%fmm, distFromPMTCenter=%fmm, pmtRadius=%fmm",
                         h.Magnitude()/I3Units::mm,
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
            
            // Ensure that the PMT window is actually inside the housing
            assert(pathLengthInOM >= pathLengthInGlass);
        }
        
        // PMT windows should not shadow each other.
        assert(intersections <= 1);
        
        return foundIntersection;
    }
    
    
}

namespace {
    // Return whether first element is greater than the second
    bool MCHitTimeLess(const I3MCPE &elem1, const I3MCPE &elem2)
    {
        return elem1.time < elem2.time;
    }
    
    I3ParticleID GetParticleID(const I3Photon &p)
    {
        I3ParticleID id;
        
        id.majorID = p.GetParticleMajorID();
        id.minorID = p.GetParticleMinorID();
        
        return id;
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
    I3MCPESeriesMapPtr output_hitmap = I3MCPESeriesMapPtr(new I3MCPESeriesMap);
    
    // currently, the only reason we need the MCTree is that I3MCHit does
    // only allow setting the major/minor particle IDs using an existing
    // I3Particle instance with that ID combination.
    I3MCTreeConstPtr MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);
    if (!MCTree) log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                           MCTreeName_.c_str());
    
    for (I3PhotonSeriesMap::const_iterator om_it = input_hitmap->begin(); om_it != input_hitmap->end(); ++om_it)
    {
        const ModuleKey &key = om_it->first;
        if (ignoreModules.count(key) > 0) continue; // ignore masked DOMs

        // Find the current OM in the geometry map
        I3ModuleGeoMap::const_iterator geo_it = moduleGeoMap->find(key);
        if (geo_it == moduleGeoMap->end())
            log_fatal("Module (%i/%u) not found in the current geometry map!", key.GetString(), key.GetOM());
        const I3ModuleGeo &moduleGeo = geo_it->second;
        if (moduleGeo.GetModuleType() != I3ModuleGeo::mDOM)
            continue;
        assert(key.GetString() > 86);
        
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
            double hitTime=photon.GetTime();
            I3Direction rotatedPmtDir;
            int hitPmtNum = FindHitPMT(photon, 
                                       key,
                                       moduleGeo,
                                       *omGeoMap,
                                       checkPMTNumbers,
                                       glassThickness_,
                                       DOMOversizeFactor_,
                                       DOMPancakeFactor_,
                                       hitTime,
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

            const double hit_cosangle = -pmtGeo.GetDirection()*photon.GetDir();
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
            
#if 0
            const double glassGelSurvival_fac =
                std::exp(-glassThickness_/glassAbsorptionLength_->GetValue(wlen)
                         -gelThickness/gelAbsorptionLength_->GetValue(wlen));
#endif
            
            const double qe_fac = pmtWavelengthAcceptance_->GetValue(wlen);
            
            // the angular acceptance is defined in the usual way: multiply the acceptance with
            // the geometrical area of the PMT and get the "effective" area for the current angle.
            // We already do the tracking to each PMT, so a geometrical factor of cos(theta) is
            // already there from the probability to hit a PMT. (The projected area of a PMT
            // scales with cos(theta).) But this geometrical scaling is ALSO included in the 
            // angular acceptance factor from the geometry. So we have to get rid of that first.
            // This means that after the code knows that a PMT is hit, the geometrical acceptance
            // should be 1 if the acceptance factor from the geometry is cos(theta).
            // const double ang_fac = pmtAngularAcceptance_->GetValue(hit_cosangle)/std::fabs(hit_cosangle);
            
            // calculate the measurement probability
            double measurement_prob = photon.GetWeight();
            // measurement_prob *= glassGelSurvival_fac;
            measurement_prob *= qe_fac;
            // measurement_prob *= ang_fac;
            
            if (measurement_prob > 1.)
                log_fatal("measurement_prob > 1 (it's %f): cannot continue. your hit weights are too high. (weight=%f, wlen=%fnm)",
                          measurement_prob, photon.GetWeight(),
                          photon.GetWavelength()/I3Units::nanometer);
            
            if (measurement_prob <= randomService_->Uniform()) continue;
            
            // find the particle
            I3MCTree::const_iterator it = MCTree->find(GetParticleID(photon));
	    if (it==MCTree->end())
                log_fatal("Particle with id maj=%" PRIu64 ", min=%i does not exist in MC tree, but we have a photon that claims it was created by that particle..",
                          photon.GetParticleMajorID(), photon.GetParticleMinorID());
            const I3Particle &particle = *(it);
            
            
            // get the hit series into which we are going to insert the hit
            I3MCPESeries &hitSeries = output_hitmap->insert(std::make_pair(pmtKey, I3MCPESeries())).first->second;
            
            // add a new hit
            hitSeries.push_back(I3MCPE(particle));
            I3MCPE &hit = hitSeries.back();
            
            // fill in all information
            hit.time=hitTime;
            hit.npe=1;
        }
        
    }

    // sort the hit vectors for each PMT
    for (I3MCPESeriesMap::iterator pmtIt = output_hitmap->begin();
         pmtIt != output_hitmap->end(); ++pmtIt)
    {
        I3MCPESeries &hitSeries = pmtIt->second;
        
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

