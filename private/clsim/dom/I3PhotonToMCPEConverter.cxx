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
 * @file I3PhotonToMCPEConverter.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <algorithm>

#include "clsim/dom/I3PhotonToMCPEConverter.h"

#include <boost/foreach.hpp>

#include "clsim/I3Photon.h"

#include "phys-services/I3SummaryService.h"

#include "dataclasses/geometry/I3OMGeo.h"
#include "dataclasses/geometry/I3ModuleGeo.h"

#include "simclasses/I3MCPE.h"
#include "dataclasses/physics/I3MCTree.h"

#include "dataclasses/I3Constants.h"

// The module
I3_MODULE(I3PhotonToMCPEConverter);

I3PhotonToMCPEConverter::I3PhotonToMCPEConverter(const I3Context& context) 
: I3ConditionalModule(context)
{
    AddParameter("RandomService",
                 "A random number generating service (derived from I3RandomService).",
                 randomService_);

    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object.",
                 inputPhotonSeriesMapName_);

    outputMCPESeriesMapName_="MCPESeriesMap";
    AddParameter("OutputMCPESeriesMapName",
                 "Name of the output I3MCPESeries frame object.",
                 outputMCPESeriesMapName_);

    MCTreeName_="I3MCTree";
    AddParameter("MCTreeName",
                 "Name of the I3MCTree frame object. All photon particle IDs are checked against this tree.",
                 MCTreeName_);

    AddParameter("WavelengthAcceptance",
                 "Wavelength acceptance of the (D)OM as a I3WlenDependedValue object.",
                 wavelengthAcceptance_);

    AddParameter("AngularAcceptance",
                 "Angular acceptance of the (D)OM as a I3WlenDependedValue object.",
                 angularAcceptance_);

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

    DOMRadiusWithoutOversize_=0.16510*I3Units::m; // 13 inch diameter
    AddParameter("DOMRadiusWithoutOversize",
                 "Specifiy the DOM radius. Do not include oversize factors here.",
                 DOMRadiusWithoutOversize_);

    defaultRelativeDOMEfficiency_=1.;
    AddParameter("DefaultRelativeDOMEfficiency",
                 "Default relative efficiency. This value is used if no entry is available from I3Calibration.",
                 defaultRelativeDOMEfficiency_);

    replaceRelativeDOMEfficiencyWithDefault_=false;
    AddParameter("ReplaceRelativeDOMEfficiencyWithDefault",
                 "Always use the default relative efficiency, ignore other values from I3Calibration.",
                 replaceRelativeDOMEfficiencyWithDefault_);

    ignoreDOMsWithoutDetectorStatusEntry_=false;
    AddParameter("IgnoreDOMsWithoutDetectorStatusEntry",
                 "Do not generate hits for OMKeys not found in the I3DetectorStatus.I3DOMStatusMap",
                 ignoreDOMsWithoutDetectorStatusEntry_);

    onlyWarnAboutInvalidPhotonPositions_=false;
    AddParameter("OnlyWarnAboutInvalidPhotonPositions",
                 "Make photon position/radius check a warning only (instead of a fatal condition)",
                 onlyWarnAboutInvalidPhotonPositions_);

    // add an outbox
    AddOutBox("OutBox");
    
    numGeneratedHits_=0;

}

I3PhotonToMCPEConverter::~I3PhotonToMCPEConverter()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3PhotonToMCPEConverter::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("RandomService", randomService_);

    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputMCPESeriesMapName", outputMCPESeriesMapName_);

    GetParameter("MCTreeName", MCTreeName_);

    GetParameter("WavelengthAcceptance", wavelengthAcceptance_);
    GetParameter("AngularAcceptance", angularAcceptance_);

    GetParameter("DOMOversizeFactor", DOMOversizeFactor_);
    GetParameter("DOMPancakeFactor", DOMPancakeFactor_);
    GetParameter("DOMRadiusWithoutOversize", DOMRadiusWithoutOversize_);

    GetParameter("DefaultRelativeDOMEfficiency", defaultRelativeDOMEfficiency_);
    GetParameter("ReplaceRelativeDOMEfficiencyWithDefault", replaceRelativeDOMEfficiencyWithDefault_);
    GetParameter("IgnoreDOMsWithoutDetectorStatusEntry", ignoreDOMsWithoutDetectorStatusEntry_);

    GetParameter("OnlyWarnAboutInvalidPhotonPositions", onlyWarnAboutInvalidPhotonPositions_);

    if (DOMOversizeFactor_ != DOMPancakeFactor_)
        log_warn("You chose \"DOMOversizeFactor\" and \"DOMPancakeFactor\" to be different. Be sure you know whot you are doing! You probably don't want this.");
    
    if (replaceRelativeDOMEfficiencyWithDefault_)
    {
        if (std::isnan(defaultRelativeDOMEfficiency_))
            log_fatal("You need to set \"DefaultRelativeDOMEfficiency\" to a value other than NaN if you enabled \"ReplaceRelativeDOMEfficiencyWithDefault\"");
    }

    if (defaultRelativeDOMEfficiency_<0.) 
        log_fatal("The \"DefaultRelativeDOMEfficiency\" parameter must not be < 0!");
    
    if (!wavelengthAcceptance_)
        log_fatal("The \"WavelengthAcceptance\" parameter must not be empty.");
    if (!angularAcceptance_)
        log_fatal("The \"AngularAcceptance\" parameter must not be empty.");
    
    if (!wavelengthAcceptance_->HasNativeImplementation())
        log_fatal("The wavelength acceptance function must have a native (i.e. non-OpenCL) implementation!");
    if (!angularAcceptance_->HasNativeImplementation())
        log_fatal("The angular acceptance function must have a native (i.e. non-OpenCL) implementation!");
    
    
    if (!randomService_) {
        log_info("No random service provided as a parameter, trying to get one from the context..");
        randomService_ = context_.Get<I3RandomServicePtr>();
    }
    
    if (!randomService_)
        log_fatal("No random service provided using the \"RandomService\" parameter (and there is none installed on the context).");
    
}

void I3PhotonToMCPEConverter::DetectorStatus(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
        
    I3DetectorStatusConstPtr detectorStatus = frame->Get<I3DetectorStatusConstPtr>("I3DetectorStatus");
    if (!detectorStatus)
        log_fatal("detector status frame does not have an I3DetectorStatus entry");
    
    // store it for later
    status_ = detectorStatus;

    PushFrame(frame);
}

void I3PhotonToMCPEConverter::Calibration(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3CalibrationConstPtr calibration = frame->Get<I3CalibrationConstPtr>("I3Calibration");
    if (!calibration)
        log_fatal("calibration frame does not have an I3Calibration entry");
    
    // store it for later
    calibration_ = calibration;
    
    PushFrame(frame);
}

namespace {
    // Return whether first element is greater than the second
    bool MCPETimeLess(const I3MCPE &elem1, const I3MCPE &elem2)
    {
        return elem1.time < elem2.time;
    }
}

#ifdef IS_Q_FRAME_ENABLED
void I3PhotonToMCPEConverter::DAQ(I3FramePtr frame)
#else
void I3PhotonToMCPEConverter::Physics(I3FramePtr frame)
#endif
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    // First we need to get our geometry
    I3OMGeoMapConstPtr omgeo = frame->Get<I3OMGeoMapConstPtr>("I3OMGeoMap");
    I3ModuleGeoMapConstPtr modulegeo = frame->Get<I3ModuleGeoMapConstPtr>("I3ModuleGeoMap");
    
    if (!omgeo)
        log_fatal("Missing geometry information! (No \"I3OMGeoMap\")");
    if (!modulegeo)
        log_fatal("Missing geometry information! (No \"I3ModuleGeoMap\")");
    
    if (!replaceRelativeDOMEfficiencyWithDefault_) {
        // no need to check for exitsing calibration frames if the efficiency
        // will be replaced with a default value anyway
        if (!calibration_)
            log_fatal("no Calibration frame yet, but received a Physics frame.");
    }

    if ((!status_) && (ignoreDOMsWithoutDetectorStatusEntry_))
        log_fatal("no DetectorStatus frame yet, but received a Physics frame.");
    
    I3PhotonSeriesMapConstPtr inputPhotonSeriesMap = frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if (!inputPhotonSeriesMap) {
        log_debug("Frame does not contain an I3PhotonSeriesMap named \"%s\".",
                  inputPhotonSeriesMapName_.c_str());
        
        // do nothing if there is no input data
        PushFrame(frame);
        return;
    }
    
    // currently, the only reason we need the MCTree is that I3MCPE does
    // only allow setting the major/minor particle IDs using an existing
    // I3Particle instance with that ID combination.
    I3MCTreeConstPtr MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);

    // allocate the output hitSeriesMap
    I3MCPESeriesMapPtr outputMCPESeriesMap(new I3MCPESeriesMap());
    
    std::map<std::pair<uint64_t, int>, const I3Particle *> mcTreeIndex;
    if (MCTree) {
        // build an index into the I3MCTree
        for (I3MCTree::iterator it = MCTree->begin();
             it != MCTree->end(); ++it)
        {
            const I3Particle &particle = *it;
            mcTreeIndex.insert(std::make_pair(std::make_pair(particle.GetMajorID(), particle.GetMinorID()), &particle));
        }
    }    
    
#if !defined(HAS_MULTIPMT_SUPPORT)
    // DOM is looking downwards
    const double DOMDir_x = 0.;
    const double DOMDir_y = 0.;
    const double DOMDir_z = -1.;
#endif
    
    BOOST_FOREACH(const I3PhotonSeriesMap::value_type &it, *inputPhotonSeriesMap)
    {
        const ModuleKey &module_key = it.first;
        // assume this is IceCube (i.e. one PMT with index 0 per DOM)
        const OMKey key(module_key.GetString(), module_key.GetOM(), 0);        
        const I3PhotonSeries &photons = it.second;

        if (ignoreDOMsWithoutDetectorStatusEntry_) {
            std::map<OMKey, I3DOMStatus>::const_iterator om_stat = status_->domStatus.find(key);
            if (om_stat==status_->domStatus.end()) continue; // ignore it
            if (om_stat->second.pmtHV==0.) continue; // ignore pmtHV==0
        }
        
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

        // this module assumes that all DOMs are IceCube-style with a single PMT per DOM
        if ((std::abs(om.position.GetX() - module.GetPos().GetX()) > .01*I3Units::mm) ||
            (std::abs(om.position.GetY() - module.GetPos().GetY()) > .01*I3Units::mm) ||
            (std::abs(om.position.GetZ() - module.GetPos().GetZ()) > .01*I3Units::mm))
            log_fatal("Module(%i/%u) has a PMT that is not in the center of the DOM!",
                      module_key.GetString(), module_key.GetOM());
        
        const I3Direction pmtDir = om.GetDirection();
        const I3Direction domDir = module.GetDir();
        
        const double DOMDir_x = pmtDir.GetX();
        const double DOMDir_y = pmtDir.GetY();
        const double DOMDir_z = pmtDir.GetZ();
        
        if ((std::abs(DOMDir_x - domDir.GetX()) > 1e-5) ||
            (std::abs(DOMDir_y - domDir.GetY()) > 1e-5) ||
            (std::abs(DOMDir_z - domDir.GetZ()) > 1e-5))
            log_fatal("PMT and DOM directions are not aligned!");
        
        // Find the current OM in the calibration map
        
        // relative DOM efficiency from calibration
        double efficiency_from_calibration=NAN;

        if (replaceRelativeDOMEfficiencyWithDefault_)
        {
            efficiency_from_calibration=defaultRelativeDOMEfficiency_;
        }
        else 
        {
            if (!calibration_) {
                if (isnan(defaultRelativeDOMEfficiency_)) {
                    log_fatal("There is no valid calibration! (Consider setting \"DefaultRelativeDOMEfficiency\" != NaN)");
                } else {
                    efficiency_from_calibration = defaultRelativeDOMEfficiency_;
                    log_debug("OM (%i/%u): efficiency_from_calibration=%g (default (I3Calibration not found))",
                              key.GetString(), key.GetOM(),
                              efficiency_from_calibration);
                }
            } else {
                std::map<OMKey, I3DOMCalibration>::const_iterator cal_it = calibration_->domCal.find(key);
                if (cal_it == calibration_->domCal.end()) {
                    if (isnan(defaultRelativeDOMEfficiency_)) {
                        log_fatal("OM (%i/%u) not found in the current calibration map! (Consider setting \"DefaultRelativeDOMEfficiency\" != NaN)", key.GetString(), key.GetOM());
                    } else {
                        efficiency_from_calibration = defaultRelativeDOMEfficiency_;
                        log_debug("OM (%i/%u): efficiency_from_calibration=%g (default (no calib))",
                                  key.GetString(), key.GetOM(),
                                  efficiency_from_calibration);
                    }
                } else {
                    const I3DOMCalibration &domCalibration = cal_it->second;
                    efficiency_from_calibration=domCalibration.GetRelativeDomEff();
                    
                    if (isnan(efficiency_from_calibration)) {
                        if (isnan(defaultRelativeDOMEfficiency_)) {
                            log_fatal("OM (%i/%u) found in the current calibration map, but it is NaN! (Consider setting \"DefaultRelativeDOMEfficiency\" != NaN)", key.GetString(), key.GetOM());
                        } else {                
                            efficiency_from_calibration = defaultRelativeDOMEfficiency_;
                            log_debug("OM (%i/%u): efficiency_from_calibration=%g (default (was: NaN))",
                                      key.GetString(), key.GetOM(),
                                      efficiency_from_calibration);
                        }
                    } else {
                        log_debug("OM (%i/%u): efficiency_from_calibration=%g",
                                  key.GetString(), key.GetOM(),
                                  efficiency_from_calibration);
                    }
                }
            }
        }
        
        // a pointer to the output vector. The vector will be allocated 
        // by the map, this is merely a pointer to it in case we have multiple
        // hits per OM.
        I3MCPESeries *hits = NULL;

        BOOST_FOREACH(const I3Photon &photon, photons)
        {
            double hitProbability = photon.GetWeight();
            if (hitProbability < 0.) log_fatal("Photon with negative weight found.");
            if (hitProbability == 0.) continue;

            const double dx=photon.GetDir().GetX();
            const double dy=photon.GetDir().GetY();
            const double dz=photon.GetDir().GetZ();
            const double px=om.position.GetX()-photon.GetPos().GetX();
            const double py=om.position.GetY()-photon.GetPos().GetY();
            const double pz=om.position.GetZ()-photon.GetPos().GetZ();
            const double pr2 = px*px + py*py + pz*pz;
            
            double photonCosAngle = -(dx * DOMDir_x +
                                      dy * DOMDir_y +
                                      dz * DOMDir_z);
            photonCosAngle = std::max(-1., std::min(1., photonCosAngle));
            
            const double distFromDOMCenter = std::sqrt(pr2);

            // do this only if DOMs are spherical
            if (DOMPancakeFactor_ == 1.)
            {
                // sanity check: are photons on the OM's surface?
                if (std::abs(distFromDOMCenter - DOMOversizeFactor_*DOMRadiusWithoutOversize_) > 3.*I3Units::cm) {
                    if (onlyWarnAboutInvalidPhotonPositions_) {
                        log_warn("distance not %f*%f=%fmm.. it is %fmm (diff=%gmm) (OMKey=(%i,%u) (photon @ pos=(%g,%g,%g)m) (DOM @ pos=(%g,%g,%g)m)",
                                 DOMOversizeFactor_,
                                 DOMRadiusWithoutOversize_/I3Units::mm,
                                 DOMOversizeFactor_*DOMRadiusWithoutOversize_/I3Units::mm,
                                 distFromDOMCenter/I3Units::mm,
                                 (distFromDOMCenter-DOMOversizeFactor_*DOMRadiusWithoutOversize_)/I3Units::mm,
                                 key.GetString(), key.GetOM(),
                                 photon.GetPos().GetX()/I3Units::m,
                                 photon.GetPos().GetY()/I3Units::m,
                                 photon.GetPos().GetZ()/I3Units::m,
                                 om.position.GetX()/I3Units::m,
                                 om.position.GetY()/I3Units::m,
                                 om.position.GetZ()/I3Units::m
                                 );
                    } else {
                        log_fatal("distance not %f*%f=%fmm.. it is %fmm (diff=%gmm) (OMKey=(%i,%u) (photon @ pos=(%g,%g,%g)m) (DOM @ pos=(%g,%g,%g)m)",
                                  DOMOversizeFactor_,
                                  DOMRadiusWithoutOversize_/I3Units::mm,
                                  DOMOversizeFactor_*DOMRadiusWithoutOversize_/I3Units::mm,
                                  distFromDOMCenter/I3Units::mm,
                                  (distFromDOMCenter-DOMOversizeFactor_*DOMRadiusWithoutOversize_)/I3Units::mm,
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
            }
            
            // sanity check for unscattered photons: is their direction ok
            // w.r.t. the vector from emission to detection?
            if (photon.GetNumScattered()==0)
            {
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
                              ppl
                              );
                }
            }
            
#ifndef NDEBUG
            const double photonAngle = std::acos(photonCosAngle);
            log_trace("Photon (lambda=%fnm, angle=%fdeg, dist=%fm) has weight %g",
                     photon.GetWavelength()/I3Units::nanometer,
                     photonAngle/I3Units::deg,
                     distFromDOMCenter/I3Units::m,
                     hitProbability);
#endif
            
            hitProbability *= wavelengthAcceptance_->GetValue(photon.GetWavelength());
            log_trace("After wlen acceptance: prob=%g (wlen acceptance is %f)",
                     hitProbability, wavelengthAcceptance_->GetValue(photon.GetWavelength()));

            hitProbability *= angularAcceptance_->GetValue(photonCosAngle);
            log_trace("After wlen&angular acceptance: prob=%g (angular acceptance is %f)",
                      hitProbability, angularAcceptance_->GetValue(photonCosAngle));

            hitProbability *= efficiency_from_calibration;
            log_trace("After efficiency from calibration: prob=%g (efficiency_from_calibration=%f)",
                      hitProbability, efficiency_from_calibration);

            if (hitProbability > 1.) {
                log_warn("hitProbability==%f > 1: your hit weights are too high. (hitProbability-1=%f)", hitProbability, hitProbability-1.);

                double hitProbability = photon.GetWeight();

                const double photonAngle = std::acos(photonCosAngle);
                log_warn("Photon (lambda=%fnm, angle=%fdeg, dist=%fm) has weight %g, 1/weight %g",
                         photon.GetWavelength()/I3Units::nanometer,
                         photonAngle/I3Units::deg,
                         distFromDOMCenter/I3Units::m,
                         hitProbability,
                         1./hitProbability);

                hitProbability *= wavelengthAcceptance_->GetValue(photon.GetWavelength());
                log_warn("After wlen acceptance: prob=%g (wlen acceptance is %f)",
                         hitProbability, wavelengthAcceptance_->GetValue(photon.GetWavelength()));

                hitProbability *= angularAcceptance_->GetValue(photonCosAngle);
                log_warn("After wlen&angular acceptance: prob=%g (angular acceptance is %f)",
                          hitProbability, angularAcceptance_->GetValue(photonCosAngle));

                hitProbability *= efficiency_from_calibration;
                log_warn("After efficiency from calibration: prob=%g (efficiency_from_calibration=%f)",
                          hitProbability, efficiency_from_calibration);
                
                log_fatal("cannot continue.");
            }
            
            // does it survive?
            if (hitProbability <= randomService_->Uniform()) continue;

            // find the particle
            const I3Particle *particle = NULL;
            
            if ((photon.GetParticleMajorID() != 0) || (photon.GetParticleMinorID() != 0))
            {
                // index (0,0) is used for flasher photons, set no hit particle for those
                std::map<std::pair<uint64_t, int>, const I3Particle *>::const_iterator it = 
                mcTreeIndex.find(std::make_pair(photon.GetParticleMajorID(), photon.GetParticleMinorID()));
                if (it==mcTreeIndex.end())
                    log_fatal("Particle with id maj=%" PRIu64 ", min=%i does not exist in MC tree, but we have a photon that claims it was created by that particle..",
                              photon.GetParticleMajorID(), photon.GetParticleMinorID());
                particle = it->second;
            }
        
            // allocate the output vector if not already done
            if (!hits) hits = &(outputMCPESeriesMap->insert(std::make_pair(key, I3MCPESeries())).first->second);

            // correct timing for oversized DOMs
            double correctedTime = photon.GetTime();
            {
                const double dot = px*dx + py*dy + pz*dz;
                const double bringForward = dot*(1.-DOMPancakeFactor_/DOMOversizeFactor_);
                correctedTime += bringForward/photon.GetGroupVelocity();
            }
            
            // add a new hit
			if(particle)
				hits->push_back(I3MCPE(*particle));
			else
				hits->push_back(I3MCPE());
            I3MCPE &hit = hits->back();
            
            // fill in all information
            hit.time=correctedTime;
            hit.npe=1;
        }
        
        if (hits) {
            // sort the photons in each hit series by time
            std::sort(hits->begin(), hits->end(), MCPETimeLess);
            
            // keep track of the number of hits generated
            numGeneratedHits_ += static_cast<uint64_t>(hits->size());
        }
    }
    
    // store the output I3MCPESeriesMap
    frame->Put(outputMCPESeriesMapName_, outputMCPESeriesMap);
    
    // that's it!
    PushFrame(frame);
}

void I3PhotonToMCPEConverter::Finish()
{
    // add some summary information to a potential I3SummaryService
    I3SummaryServicePtr summary = context_.Get<I3SummaryServicePtr>();
    if (summary) {
        const std::string prefix = "I3PhotonToMCPEConverter_" + GetName() + "_";
        
        (*summary)[prefix+"NumGeneratedHits"] = numGeneratedHits_;
    }
    
}
