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
 * $Id: I3PhotonToMCHitConverterForWOMs.cxx 149902 2016-09-09 03:15:41Z cweaver $
 *
 * @file I3PhotonToMCHitConverterForWOMs.cxx
 * @version $Revision: 149902 $
 * @date $Date: 2016-09-09 05:15:41 +0200 (Fri, 09 Sep 2016) $
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <algorithm>
#include <cmath>

//#include "clsim/dom/I3PhotonToMCHitConverterForWOMs.h"
#include "clsim/dom/I3PhotonToMCHitConverterForWOMs.h"

#include <boost/foreach.hpp>

#include "simclasses/I3Photon.h"

#include "phys-services/I3SummaryService.h"

#include "dataclasses/geometry/I3OMGeo.h"
#include "dataclasses/geometry/I3ModuleGeo.h"

#include "simclasses/I3MCPE.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/physics/I3ParticleID.h"

#include "dataclasses/I3Constants.h"

// The module
I3_MODULE(I3PhotonToMCHitConverterForWOMs);

I3PhotonToMCHitConverterForWOMs::I3PhotonToMCHitConverterForWOMs(const I3Context& context) 
: I3ConditionalModule(context)
{
    AddParameter("RandomService",
                 "A random number generating service (derived from I3RandomService).",
                 randomService_);

    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object.",
                 inputPhotonSeriesMapName_);

    outputMCHitSeriesMapName_="MCHitSeriesMap";
    AddParameter("OutputMCHitSeriesMapName",
                 "Name of the output I3MCPESeries frame object.",
                 outputMCHitSeriesMapName_);

    MCTreeName_="I3MCTree";
    AddParameter("MCTreeName",
                 "Name of the I3MCTree frame object. All photon particle IDs are checked against this tree.",
                 MCTreeName_);

    AddParameter("HeightAcceptance",
                 "Height acceptance of the WOM as I3tobedetermined object.",
                 HeightAcceptance_);

    AddParameter("PMTHeightAcceptance",
                 "PMT Height acceptance of the WOM as I3tobedetermined object.",
                 PMTHeightAcceptance_);

    AddParameter("AngularAcceptance",
                 "Angular acceptance of the WOM as a I3WlenDependedValue object.\n" 
                 "This includes Propagation from ice to glass to air to PMMA with Fresnel",
                 AngularAcceptance_);

    DOMOversizeFactor_=1.;
    AddParameter("DOMOversizeFactor",
                 "Specifiy the \"oversize factor\" (i.e. DOM radius scaling factor) you used during the CLSim run.\n"
                 "The photon arrival times will be corrected. In practice this means your large spherical DOMs will\n"
                 "become ellipsoids.",
                 DOMOversizeFactor_);

    DOMRadius_=0.114*I3Units::m; // 11.4 centimeter diameter
    AddParameter("DOMRadius",
                 "Specifiy the DOM radius. Do not include oversize factors here.",
                 DOMRadius_);

    OMHeight_=0.9*I3Units::m; // 90 centimeter WOM height
    AddParameter("OMHeight",
                 "Specifiy the OM hight in form of a cylindrical extension. Do not include oversize factors here.",
                 OMHeight_);

    AddParameter("GlassAbsorptionLength",
                 "Specifiy the absorptionlength of the glass used for the pressure vessel as a I3WlenDependedValue object..",
                 GlassAbsorptionLength_);

    GlassThickness_=0.01*I3Units::m; // I guessed about 1 centimeter
    AddParameter("GlassThickness",
                 "Specifiy the thickness of the glass used for the pressure vessel. A refractive index of 1.5 is assumed for refractive calculations.",
                 GlassThickness_);

    AddParameter("WLSPropagationEfficiency",
                 "Combined efficiency for photon capture of the WLS tube, propagation and the PMT for the WOM as a I3WlenDependedValue object.",
                 WLSPropagationEfficiency_);

    onlyWarnAboutInvalidPhotonPositions_=false;
    AddParameter("OnlyWarnAboutInvalidPhotonPositions",
                 "Make photon position/radius check a warning only (instead of a fatal condition)",
                 onlyWarnAboutInvalidPhotonPositions_);

    // add an outbox
    AddOutBox("OutBox");
    
    numGeneratedHits_=0;

}

I3PhotonToMCHitConverterForWOMs::~I3PhotonToMCHitConverterForWOMs()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3PhotonToMCHitConverterForWOMs::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("RandomService", randomService_);

    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputMCHitSeriesMapName", outputMCHitSeriesMapName_);

    GetParameter("MCTreeName", MCTreeName_);

    GetParameter("HeightAcceptance", HeightAcceptance_);
    GetParameter("PMTHeightAcceptance", PMTHeightAcceptance_);
    GetParameter("AngularAcceptance", AngularAcceptance_);

    GetParameter("DOMOversizeFactor", DOMOversizeFactor_);
    GetParameter("DOMRadius", DOMRadius_);
    GetParameter("OMHeight", OMHeight_);

    GetParameter("GlassAbsorptionLength", GlassAbsorptionLength_);
    GetParameter("GlassThickness", GlassThickness_);
    GetParameter("WLSPropagationEfficiency", WLSPropagationEfficiency_);

    GetParameter("OnlyWarnAboutInvalidPhotonPositions", onlyWarnAboutInvalidPhotonPositions_);


    //parameter check

    if (!HeightAcceptance_) 
        log_fatal("The \"HeightAcceptance\" parameter must be given!");

    if (!PMTHeightAcceptance_) 
        log_fatal("The \"PMTHeightAcceptance\" parameter must be given!");
    
    if (!AngularAcceptance_)
        log_fatal("The \"AngularAcceptance\" parameter must not be empty.");

    if (!WLSPropagationEfficiency_)
        log_fatal("The \"WLSPropagationEfficiency_\" parameter must not be empty.");

    if (!GlassAbsorptionLength_)
        log_fatal("The \"GlassAbsorptionLength_\" parameter must not be empty.");
        
    if (!HeightAcceptance_->HasNativeImplementation())
        log_fatal("The height acceptance function must have a native (i.e. non-OpenCL) implementation!");
    if (!PMTHeightAcceptance_->HasNativeImplementation())
        log_fatal("The PMT height acceptance function must have a native (i.e. non-OpenCL) implementation!");
    if (!AngularAcceptance_->HasNativeImplementation())
        log_fatal("The angular acceptance function must have a native (i.e. non-OpenCL) implementation!");
    if (!WLSPropagationEfficiency_->HasNativeImplementation())
        log_fatal("The WLS propagation efficiency function must have a native (i.e. non-OpenCL) implementation!");
    if (!GlassAbsorptionLength_->HasNativeImplementation())
        log_fatal("The glass absorption length function must have a native (i.e. non-OpenCL) implementation!");
    
    if (!randomService_) {
        log_info("No random service provided as a parameter, trying to get one from the context..");
        randomService_ = context_.Get<I3RandomServicePtr>();
    }
    
    if (!randomService_)
        log_fatal("No random service provided using the \"RandomService\" parameter (and there is none installed on the context).");
    
}


namespace {
    // Return whether first element is greater than the second
    bool MCPETimeLess(const I3MCPE &elem1, const I3MCPE &elem2)
    {
        return elem1.time < elem2.time;
    }
}

namespace {
    // Return whether first element is greater than the second
    bool MCHitTimeLess(const I3MCPE &elem1, const I3MCPE &elem2)
    {
        return elem1.time < elem2.time;
    }
}



void I3PhotonToMCHitConverterForWOMs::DAQ(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    // First we need to get our geometry
    I3OMGeoMapConstPtr omgeo = frame->Get<I3OMGeoMapConstPtr>("I3OMGeoMap");
    I3ModuleGeoMapConstPtr modulegeo = frame->Get<I3ModuleGeoMapConstPtr>("I3ModuleGeoMap");
    
    if (!omgeo)
        log_fatal("Missing geometry information! (No \"I3OMGeoMap\")");
    if (!modulegeo)
        log_fatal("Missing geometry information! (No \"I3ModuleGeoMap\")");

    I3PhotonSeriesMapConstPtr inputPhotonSeriesMap = frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if (!inputPhotonSeriesMap) {
        log_debug("Frame does not contain an I3PhotonSeriesMap named \"%s\".",
                  inputPhotonSeriesMapName_.c_str());
        
        // do nothing if there is no input data
        PushFrame(frame);
        return;
    }

    //get this part from GCD file in the future
    double OM_Radius = DOMOversizeFactor_*DOMRadius_;
    double OM_Height = DOMOversizeFactor_*OMHeight_;
    double WLS_decaytime = 5*I3Units::ns;

    const double distanceAccuracy = 3.*I3Units::cm;
    // currently, the only reason we need the MCTree is that I3MCPE does
    // only allow setting the major/minor particle IDs using an existing
    // I3Particle instance with that ID combination.
    I3MCTreeConstPtr MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);


    // allocate the output hitSeriesMap
    I3MCPESeriesMapPtr outputMCHitSeriesMap(new I3MCPESeriesMap());
    
    BOOST_FOREACH(const I3PhotonSeriesMap::value_type &it, *inputPhotonSeriesMap)
    {
        const ModuleKey &module_key = it.first;
        // assume this is IceCube (i.e. one PMT with index 0 per DOM)
        const OMKey key(module_key.GetString(), module_key.GetOM(), 0);        
        const I3PhotonSeries &photons = it.second;
        
        // Find the current OM in the omgeo map
        I3OMGeoMap::const_iterator geo_it = omgeo->find(key);
        if (geo_it == omgeo->end())
            log_fatal("OM (%i/%u%u) not found in the current geometry map!",
                      key.GetString(), key.GetOM(), static_cast<unsigned int>(key.GetPMT()));
        const I3OMGeo &om = geo_it->second;
        
        // Find the current OM in the module map
        I3ModuleGeoMap::const_iterator module_geo_it = modulegeo->find(module_key);
        if (module_geo_it == modulegeo->end())
            log_fatal("ModuleKey (%i/%u) not found in the current geometry map!",
                      module_key.GetString(), module_key.GetOM());
        const I3ModuleGeo &module = module_geo_it->second;
        
        // a pointer to the output vector. The vector will be allocated 
        // by the map, this is merely a pointer to it in case we have multiple
        // hits per OM.
        I3MCPESeries *hits = NULL;


        BOOST_FOREACH(const I3Photon &photon, photons)
        {
            double hitProbability = photon.GetWeight(); // Do not use weights with WOMs yet
            if (hitProbability < 0.) log_fatal("Photon with negative weight found.");
            if (hitProbability == 0.) continue;

            // refractive indexes, might be changed later by function(wvl)
            const double refractiveIndexIce = 1.33; 
            const double refractiveIndexGlass = 1.5;
            const double refractiveIndexPMMA = 1.5;

            const double dx=photon.GetDir().GetX();
            const double dy=photon.GetDir().GetY();
            const double dz=photon.GetDir().GetZ();
            const double px=module.GetPos().GetX()-photon.GetPos().GetX();
            const double py=module.GetPos().GetY()-photon.GetPos().GetY();
            const double pz=module.GetPos().GetZ()-photon.GetPos().GetZ();
            const double len2center = std::sqrt(px*px+py*py);

            // sanity check (is photon on OM surface +- 3 cm)
            bool isOnSurface = false;
            double minimalDistance=2*distanceAccuracy;
            if(std::abs(pz) < OM_Height/2. + distanceAccuracy){
                const double cylidnerDistance = std::abs(len2center-OM_Radius);
                minimalDistance = minimalDistance < cylidnerDistance ? minimalDistance : cylidnerDistance;
            }
            if( std::abs(pz) > OM_Height/2. - distanceAccuracy){
                const double pz2 = std::abs(pz)-OM_Height/2.;
                const double sphereDistance = std::abs(std::sqrt(px*px+py*py+pz2*pz2)-OM_Radius);
                minimalDistance = minimalDistance < sphereDistance ? minimalDistance : sphereDistance;
            }
            if( minimalDistance < distanceAccuracy){ 
                isOnSurface = true;
            }

            if(!isOnSurface){
              if (onlyWarnAboutInvalidPhotonPositions_) {
                        log_warn("Photon is not within %f mm of the OMs surface. DOMOversizeFactor=%f, DOMRadius=%fmm and %fmm, OMHeight=%fmm,  (OMKey=(%i,%u) (photon @ pos=(%g,%g,%g)m) (DOM @ pos=(%g,%g,%g)m, distance is cylinder %fmm or sphere %fmm)",
                                 distanceAccuracy/I3Units::mm,
                                 DOMOversizeFactor_,
                                 DOMRadius_/I3Units::mm,
                                 module.GetRadius()/I3Units::mm,
                                 OMHeight_/I3Units::mm,
                                 key.GetString(), key.GetOM(),
                                 photon.GetPos().GetX()/I3Units::m,
                                 photon.GetPos().GetY()/I3Units::m,
                                 photon.GetPos().GetZ()/I3Units::m,
                                 om.position.GetX()/I3Units::m,
                                 om.position.GetY()/I3Units::m,
                                 om.position.GetZ()/I3Units::m,
                                 std::abs(len2center-OM_Radius)/I3Units::mm,
                                 std::abs(std::sqrt(px*px+py*py+(std::abs(pz)-OM_Height/2.)*(std::abs(pz)-OM_Height/2.))-OM_Radius)/I3Units::mm
                                 );
                    } else {
                        log_fatal("Photon is not within %f mm of the OMs surface. DOMOversizeFactor=%f, DOMRadius=%fmm and %fmm, OMHeight=%fmm,  (OMKey=(%i,%u) (photon @ pos=(%g,%g,%g)m) (DOM @ pos=(%g,%g,%g)m, distance is cylinder %fmm or sphere %fmm)",
                                 distanceAccuracy/I3Units::mm,
                                 DOMOversizeFactor_,
                                 DOMRadius_/I3Units::mm,
                                 module.GetRadius()/I3Units::mm,
                                 OMHeight_/I3Units::mm,
                                 key.GetString(), key.GetOM(),
                                 photon.GetPos().GetX()/I3Units::m,
                                 photon.GetPos().GetY()/I3Units::m,
                                 photon.GetPos().GetZ()/I3Units::m,
                                 om.position.GetX()/I3Units::m,
                                 om.position.GetY()/I3Units::m,
                                 om.position.GetZ()/I3Units::m,
                                 std::abs(len2center-OM_Radius)/I3Units::mm,
                                 std::abs(std::sqrt(px*px+py*py+(std::abs(pz)-OM_Height/2.)*(std::abs(pz)-OM_Height/2.))-OM_Radius)/I3Units::mm
                                 );
                    }
            }


            // sanity check for unscattered photons: is their direction ok
            // w.r.t. the vector from emission to detection?
            if (photon.GetNumScattered()==0){
                double ppx = photon.GetPos().GetX()-photon.GetStartPos().GetX();
                double ppy = photon.GetPos().GetY()-photon.GetStartPos().GetY();
                double ppz = photon.GetPos().GetZ()-photon.GetStartPos().GetZ();
                const double ppl = std::sqrt(ppx*ppx + ppy*ppy + ppz*ppz);
                ppx/=ppl; ppy/=ppl; ppz/=ppl;
                const double cosang = dx*ppx + dy*ppy + dz*ppz;
                
                if ((cosang < 0.9) && (ppl>1.*I3Units::m)) {
                    log_fatal("unscattered photon direction is inconsistent: cos(ang)==%f, d=(%f,%f,%f), pp=(%f,%f,%f) pp_l=%f",
                              cosang,
                              dx, dy, dz,
                              ppx, ppy, ppz,
                              ppl);
                }
            }

            // Cutting away the end caps for now
            if( std::abs(pz) > OM_Height/2. ) continue;

            //get incidence angle 
            const double cosIncidenceAngle = (px*dx + py*dy)/len2center;
            //photonCosAngle = std::max(-1., std::min(1., photonCosAngle));

            // check sanity of photon direction
            if (cosIncidenceAngle < 0){
                    if (onlyWarnAboutInvalidPhotonPositions_) {
                        log_warn("Photon is directed away from the OMs surface. DOMOversizeFactor=%f, DOMRadius=%fmm, OMHeight=%fmm,  (OMKey=(%i,%u) (photon @ pos=(%g,%g,%g)m) (DOM @ pos=(%g,%g,%g)m)",
                                 DOMOversizeFactor_,
                                 DOMRadius_/I3Units::mm,
                                 OMHeight_/I3Units::mm,
                                 key.GetString(), key.GetOM(),
                                 photon.GetPos().GetX()/I3Units::m,
                                 photon.GetPos().GetY()/I3Units::m,
                                 photon.GetPos().GetZ()/I3Units::m,
                                 om.position.GetX()/I3Units::m,
                                 om.position.GetY()/I3Units::m,
                                 om.position.GetZ()/I3Units::m
                                 );
                    } else {
                        log_fatal("Photon is directed away from the OMs surface. DOMOversizeFactor=%f, DOMRadius=%fmm, OMHeight=%fmm,  (OMKey=(%i,%u) (photon @ pos=(%g,%g,%g)m) (DOM @ pos=(%g,%g,%g)m)",
                                 DOMOversizeFactor_,
                                 DOMRadius_/I3Units::mm,
                                 OMHeight_/I3Units::mm,
                                 key.GetString(), key.GetOM(),
                                 photon.GetPos().GetX()/I3Units::m,
                                 photon.GetPos().GetY()/I3Units::m,
                                 photon.GetPos().GetZ()/I3Units::m,
                                 om.position.GetX()/I3Units::m,
                                 om.position.GetY()/I3Units::m,
                                 om.position.GetZ()/I3Units::m
                                 );
                    }
            }


            //Check angular acceptance by Fresnel
            hitProbability *= AngularAcceptance_->GetValue(cosIncidenceAngle);// includes ice to glass to air to pmma
            if(hitProbability<=0.) continue; // avoid photon angles that are totally reflected 

            //Propagation through pressure glass vessel, approximated by flat surface
            const double cosAngleInGlass = std::cos(std::asin(std::sin(std::acos(cosIncidenceAngle))*refractiveIndexIce/refractiveIndexGlass)); // Snell's law
            hitProbability *= std::exp(-GlassThickness_/(cosAngleInGlass*GlassAbsorptionLength_->GetValue(photon.GetWavelength()))); // do this better at some point to include curvature of the pressure vessel
            
            //Add acceptance for WLS capture, shifting, tube capture, propagation and PMT acceptance
            hitProbability *= WLSPropagationEfficiency_->GetValue(photon.GetWavelength());

            //Add acceptance for Z-position along the WOM
            hitProbability *= HeightAcceptance_->GetValue(pz+OMHeight_/2.);

            if (hitProbability > 1.) {
                log_warn("hitProbability==%f > 1: your hit weights are too high. (Photon weight=%f)", hitProbability, photon.GetWeight());
                log_fatal("cannot continue.");
            }
            

            // does it survive?
            if (hitProbability <= randomService_->Uniform()){
                continue;
            }

            // find the particle
            const I3Particle *particle = NULL;
            
            if ((photon.GetParticleMajorID() != 0) && (photon.GetParticleMinorID() != 0))
            {
                // index (0,0) is used for flasher photons, set no hit particle for those
                if (!MCTree)
                    log_fatal_stream("Did not find a valid I3MCTree with name '" << MCTreeName_ << "'");
                particle = I3MCTreeUtils::GetParticlePtr(MCTree, photon.GetParticleID());
            }


            //determine pmt by use of PMTHeightAcceptance_
            int hitPmtNum=0;
            if (PMTHeightAcceptance_->GetValue(pz+OMHeight_/2.) <= randomService_->Uniform()){
                hitPmtNum=1;
            }

            
            const OMKey pmtKey(key.GetString(), key.GetOM(), static_cast<unsigned char>(hitPmtNum));
            I3MCPESeries &hitSeries = outputMCHitSeriesMap->insert(std::make_pair(pmtKey, I3MCPESeries())).first->second;
            
            // add a new hit
            if(particle)
                hitSeries.push_back(I3MCPE(*particle));
            else
                hitSeries.push_back(I3MCPE());
            I3MCPE &hit = hitSeries.back();
            
            // fill in all information // WLS_decaytime + z
            double timeOffset;
            double randVar=randomService_->Uniform();
            while(randVar==0){
                randVar=randomService_->Uniform();               
            }
            if(hitPmtNum==1)
                timeOffset=(pz+OMHeight_/2.)*refractiveIndexPMMA/I3Constants::c;
            else
                timeOffset=std::abs(-pz+OMHeight_/2.)*refractiveIndexPMMA/I3Constants::c;

            timeOffset+=-WLS_decaytime*std::log(randVar);
            hit.time=photon.GetTime()+timeOffset;
            hit.npe=1;   
        }
    }
    
    // sort the hit vectors for each PMT
    for (I3MCPESeriesMap::iterator pmtIt = outputMCHitSeriesMap->begin();
         pmtIt != outputMCHitSeriesMap->end(); ++pmtIt){
        I3MCPESeries &hitSeries = pmtIt->second;
        
        // now sort by time, regardless of particle ID
        std::sort(hitSeries.begin(), hitSeries.end(), MCHitTimeLess);
    }



    // store the output I3MCPESeriesMap
    frame->Put(outputMCHitSeriesMapName_, outputMCHitSeriesMap);
    
    // that's it!
    PushFrame(frame);
}

void I3PhotonToMCHitConverterForWOMs::Finish()
{
    // add some summary information to a potential I3SummaryService
    I3SummaryServicePtr summary = context_.Get<I3SummaryServicePtr>();
    if (summary) {
        const std::string prefix = "I3PhotonToMCHitConverterForWOMs_" + GetName() + "_";
        
        (*summary)[prefix+"NumGeneratedHits"] = numGeneratedHits_;
    }
    
}
