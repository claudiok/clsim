/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3CLSimModule.cxx
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

#include "clsim/I3CLSimModule.h"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/variant/get.hpp>

#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

// The module
I3_MODULE(I3CLSimModule);

I3CLSimModule::I3CLSimModule(const I3Context& context) 
: I3ConditionalModule(context),
geometryIsConfigured_(false)
{
    // define parameters
    AddParameter("RandomService",
                 "A random number generating service (derived from I3RandomService).",
                 randomService_);

    AddParameter("MediumProperties",
                 "An instance of I3CLSimMediumProperties describing the ice/water properties.",
                 mediumProperties_);

    maxNumParallelEvents_=1000;
    AddParameter("MaxNumParallelEvents",
                 "Maximum number of events that will be processed by the GPU in parallel.",
                 MCTreeName_);

    MCTreeName_="I3MCTree";
    AddParameter("MCTreeName",
                 "Name of the I3MCTree frame object. All particles except neutrinos will be read from this tree.",
                 MCTreeName_);

    photonSeriesMapName_="PropagatedPhotons";
    AddParameter("PhotonSeriesMapName",
                 "Name of the I3CLSimPhotonSeriesMap frame object that will be written to the frame.",
                 photonSeriesMapName_);

    ignoreMuons_=false;
    AddParameter("IgnoreMuons",
                 "If set to True, muons will not be propagated.",
                 ignoreMuons_);


    // add an outbox
    AddOutBox("OutBox");

}

I3CLSimModule::~I3CLSimModule()
{
    log_trace("%s", __PRETTY_FUNCTION__);

}

namespace {
    I3CLSimStepToPhotonConverterOpenCLPtr initializeOpenCL(I3RandomServicePtr rng,
                                                     I3CLSimSimpleGeometryFromI3GeometryPtr geometry,
                                                     I3CLSimMediumPropertiesPtr medium)
    {
        I3CLSimStepToPhotonConverterOpenCLPtr conv(new I3CLSimStepToPhotonConverterOpenCL(rng, true));

        log_info("available OpenCL devices:");
        shared_ptr<const std::vector<std::pair<std::string, std::string> > >
        deviceList = conv->GetDeviceList();
        if (!deviceList) log_fatal("Internal error. GetDeviceList() returned NULL.");
        if (deviceList->size() <= 0) log_fatal("No OpenCL devices available!");
        
        typedef std::pair<std::string, std::string> stringPair_t;
        
        BOOST_FOREACH(const stringPair_t &device, *deviceList)
        {
            log_info("platform: \"%s\", device: \"%s\"",
                     device.first.c_str(), device.second.c_str());
        
        }

        // do a semi-smart device selection
        stringPair_t deviceToUse;
        
        // look for a "GeForce" device first
        std::vector<std::pair<std::string, std::string> > geForceDevices;
        BOOST_FOREACH(const stringPair_t &device, *deviceList)
        {
            std::string deviceNameLowercase = boost::to_lower_copy(device.second);
            
            if ( deviceNameLowercase.find("geforce") != deviceNameLowercase.npos )
                geForceDevices.push_back(device);
        }

        if (geForceDevices.size() > 0)
        {
            if (geForceDevices.size()>1)
            {
                deviceToUse=geForceDevices[0];
                log_info("You seem to have more than one GeForce GPU. Using the first one. (\"%s\")", 
                         deviceToUse.second.c_str());
            } else {
                deviceToUse=geForceDevices[0];
                log_info("You seem to have a GeForce GPU. (\"%s\")",
                         deviceToUse.second.c_str());
            }
        }
        else
        {
            log_info("No GeForce device found. Just using the first available one.");
            deviceToUse=(*deviceList)[0];
        }
        
        log_info(" -> using OpenCL device \"%s\" on platform \"%s\".",
                 deviceToUse.second.c_str(), deviceToUse.first.c_str());

        conv->SetDeviceName(deviceToUse.first, deviceToUse.second);
        
        conv->SetMediumProperties(medium);
        conv->SetGeometry(geometry);
        
        conv->Compile();
        //log_trace("%s", conv.GetFullSource().c_str());
        
        const std::size_t maxWorkgroupSize = conv->GetMaxWorkgroupSize();
        conv->SetWorkgroupSize(maxWorkgroupSize);

        const std::size_t workgroupSize = conv->GetWorkgroupSize();
        
        // use approximately 512000 work items, convert to a multiple of the workgroup size
        const std::size_t maxNumWorkitems = (512000/workgroupSize)*workgroupSize;
        conv->SetMaxNumWorkitems(maxNumWorkitems);

        log_info("maximum workgroup size is %zu", maxWorkgroupSize);
        log_info("configured workgroup size is %zu", workgroupSize);
        log_info("maximum number of work items is %zu", maxNumWorkitems);

        conv->Initialize();
        
        return conv;
    }

    I3CLSimParticleToStepConverterGeant4Ptr initializeGeant4(I3RandomServicePtr rng,
                                                             I3CLSimMediumPropertiesPtr medium,
                                                             I3CLSimStepToPhotonConverterOpenCLPtr openCLconverter,
                                                             bool multiprocessor=false)
    {
        I3CLSimParticleToStepConverterGeant4Ptr conv
        (
         new I3CLSimParticleToStepConverterGeant4(rng->Integer(900000000)) // TODO: eventually Geant4 should use the IceTray rng!!
        );
        
        conv->SetMediumProperties(medium);
        conv->SetMaxBunchSize(openCLconverter->GetMaxNumWorkitems());
        conv->SetBunchSizeGranularity(openCLconverter->GetWorkgroupSize());
        
        //conv->SetElectronPositronMinEnergyForSecondary(1.*I3Units::GeV);
        //conv->SetElectronPositronMaxEnergyForSecondary(10.*I3Units::TeV);
        
        conv->Initialize();
        
        return conv;
    }

}

void I3CLSimModule::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("RandomService", randomService_);
    GetParameter("MediumProperties", mediumProperties_);
    GetParameter("MaxNumParallelEvents", maxNumParallelEvents_);
    GetParameter("MCTreeName", MCTreeName_);
    GetParameter("PhotonSeriesMapName", photonSeriesMapName_);
    GetParameter("IgnoreMuons", ignoreMuons_);
    
    
    if (!randomService_) log_fatal("You have to specify the \"RandomService\" parameter!");
    if (!mediumProperties_) log_fatal("You have to specify the \"MediumProperties\" parameter!");
    if (maxNumParallelEvents_ <= 0) log_fatal("Values <= 0 are invalid for the \"MaxNumParallelEvents\" parameter!");
    
    currentParticleCacheIndex_ = 0;
    geometryIsConfigured_ = false;
    totalSimulatedEnergyForFlush_ = 0.;
    totalNumParticlesForFlush_ = 0;
    
    photonSeriesMapConverter_ = I3CLSimPhotonSeriesToPhotonSeriesMapConverterPtr(new I3CLSimPhotonSeriesToPhotonSeriesMapConverter());
    
}

void I3CLSimModule::Geometry(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    if (geometryIsConfigured_)
        log_fatal("This module currently supports only a single geometry per input file.");
    
    log_info("Retrieving geometry..");
    I3GeometryConstPtr geometryObject = frame->Get<I3GeometryConstPtr>();
    if (!geometryObject) log_fatal("Geometry frame does not have an I3Geometry object!");
    
    log_info("Converting geometry..");
    geometry_ = I3CLSimSimpleGeometryFromI3GeometryPtr
    (
     new I3CLSimSimpleGeometryFromI3Geometry(43.18*I3Units::cm/2., geometryObject)
    );

    
    log_info("Initializing CLSim..");
    // initialize OpenCL
    openCLStepsToPhotonsConverter_ = initializeOpenCL(randomService_,
                                                      geometry_,
                                                      mediumProperties_);
    
    log_info("Initializing Geant4..");
    // initialize Geant4 (will set bunch sizes according to the OpenCL settings)
    geant4ParticleToStepsConverter_ = initializeGeant4(randomService_,
                                                       mediumProperties_,
                                                       openCLStepsToPhotonsConverter_,
                                                       false); // the multiprocessor version is not yet safe to use

    
    log_info("Initializing photon map converter..");
    photonSeriesMapConverter_->SetGeometry(geometry_);
    
    log_info("Initialization complete.");
    geometryIsConfigured_=true;

    PushFrame(frame);
}

namespace {
    static inline I3CLSimStepSeriesConstPtr ResultAsStepPtr(const I3CLSimParticleToStepConverter::ConversionResult_t &res)
    {
        try {
            return boost::get<I3CLSimStepSeriesConstPtr>(res);
        } catch(boost::bad_get &) {
            return I3CLSimStepSeriesConstPtr();            
        }
    }

    static inline bool ResultIsStepPtr(const I3CLSimParticleToStepConverter::ConversionResult_t &res)
    {
        try {
            I3CLSimStepSeriesConstPtr tmp = boost::get<I3CLSimStepSeriesConstPtr>(res);
            return true;
        } catch(boost::bad_get &) {
            return false;            
        }
    }

    // this is specific to the OpenCL simulator.
    // TODO: I don't like this, it's using a magic number..
    static inline OMKey OMKeyFromOpenCLSimDummy(int32_t stringAndDomID)
    {
        int32_t stringID;
        uint32_t domID;
        
        if (stringAndDomID >= 0)
        {
            stringID = stringAndDomID / 1000;
            domID = stringAndDomID % 1000;
        }
        else
        {
            // convention for negative string IDs
            stringID = -((-stringAndDomID) / 1000);
            domID = (-stringAndDomID) % 1000;
        }

        // now the OMKey can be generated
        return OMKey(stringID, domID);
    }
    
}

void I3CLSimModule::AddPhotonsToFrames(const I3CLSimPhotonSeries &photons)
{
    if (photonsForFrameList_.size() != frameList_.size())
        log_fatal("Internal error: cache sizes differ. (1)");
    if (photonsForFrameList_.size() != currentPhotonIdForFrame_.size())
        log_fatal("Internal error: cache sizes differ. (2)");
    
    BOOST_FOREACH(const I3CLSimPhoton &photon, photons)
    {
        // find identifier in particle cache
        std::map<uint32_t, particleCacheEntry>::iterator it = particleCache_.find(photon.identifier);
        if (it == particleCache_.end())
            log_fatal("Internal error: unknown particle id from OpenCL: %" PRIu32,
                      photon.identifier);
        const particleCacheEntry &cacheEntry = it->second;

        if (cacheEntry.frameListEntry >= photonsForFrameList_.size())
            log_fatal("Internal error: particle cache entry uses invalid frame cache position");
        
        //I3FramePtr &frame = frameList_[cacheEntry.frameListEntry];
        I3PhotonSeriesMap &outputPhotonMap = *(photonsForFrameList_[cacheEntry.frameListEntry]);

        // get the current photon id
        int32_t &currentPhotonId = currentPhotonIdForFrame_[cacheEntry.frameListEntry];
        
        // generate the OMKey
        const OMKey key = OMKeyFromOpenCLSimDummy(photon.dummy);
        
        // this either inserts a new vector or retrieves an existing one
        I3PhotonSeries &outputPhotonSeries = outputPhotonMap.insert(std::make_pair(key, I3PhotonSeries())).first->second;
        
        // append a new I3Photon to the list
        outputPhotonSeries.push_back(I3Photon());
        
        // get a reference to the new photon
        I3Photon &outputPhoton = outputPhotonSeries.back();

        // fill the photon data
        outputPhoton.SetTime(photon.GetTime());
        outputPhoton.SetID(currentPhotonId); // per-frame ID for every photon
        outputPhoton.SetWeight(photon.GetWeight());
        outputPhoton.SetParticleMinorID(cacheEntry.particleMinorID);
        outputPhoton.SetParticleMajorID(cacheEntry.particleMajorID);
        outputPhoton.SetCherenkovDist(photon.GetCherenkovDist());
        outputPhoton.SetCherenkovTime(NAN); // TODO: calculate cherenkov distance from wavelength
        outputPhoton.SetWavelength(photon.GetWavelength());
        outputPhoton.SetNumScattered(photon.GetNumScatters());

        outputPhoton.SetPos(I3Position(photon.GetPosX(), photon.GetPosY(), photon.GetPosZ()));
        {
            I3Direction outDir;
            outDir.SetThetaPhi(photon.GetDirTheta(), photon.GetDirPhi());
            outputPhoton.SetDir(outDir);
        }
        
        currentPhotonId++;
    }
}

void I3CLSimModule::FlushFrameCache()
{
    log_info("Flushing frame cache..");
    
    geant4ParticleToStepsConverter_->EnqueueBarrier();
    uint32_t numBunchesSentToOpenCL = 0;
    
    for (;;)
    {
        I3CLSimParticleToStepConverter::ConversionResult_t result;
        
        if (!geant4ParticleToStepsConverter_->MoreStepsAvailable())
        {
            log_debug("No steps are available right now");
            if (!geant4ParticleToStepsConverter_->BarrierActive())
            {
                log_debug("No barrier is active, we are done");
                break; // nothing to read and no barrier. We are done
            }
            
            // no steps available, but barrier is still enqueued. Geant4
            // must still be working. Try to get steps, but timeout after a while.
            log_debug("Waiting until steps become available or a timeout occurs");

            result = geant4ParticleToStepsConverter_->GetConversionResult(0.1*I3Units::s); // time in seconds
            
            if (!ResultAsStepPtr(result))
            {
                log_debug("Timeout while waiting for steps, trying again");
                continue;
            }
            
            log_debug("Steps retrieved");
        }
        else
        {
            // there are steps available, fetch them!
            log_debug("Steps available, retrieving them");
            result = geant4ParticleToStepsConverter_->GetConversionResult();
            log_debug("Steps retrieved");
        }

        if (!ResultIsStepPtr(result))
        {
            std::pair<uint32_t, I3ParticleConstPtr> particleTuple = boost::get<std::pair<uint32_t, I3ParticleConstPtr> >(result);

            log_info("Got a secondary with E = %fGeV (ignoring it: needs to be fixed!)",
                     particleTuple.second->GetEnergy()/I3Units::GeV);
            continue;
        }

        I3CLSimStepSeriesConstPtr steps = ResultAsStepPtr(result);
        if (!steps) log_fatal("Got NULL I3CLSimStepSeriesConstPtr from Geant4.");
        
        log_info("Got %zu steps from Geant4, sending them to OpenCL",
                 steps->size());
        
        openCLStepsToPhotonsConverter_->EnqueueSteps(steps, numBunchesSentToOpenCL);
        numBunchesSentToOpenCL++;
    }

    log_info("Geant4 finished, retrieving results from GPU..");

    std::size_t totalNumOutPhotons=0;
    
    for (uint32_t i=0;i<numBunchesSentToOpenCL;++i)
    {
        //typedef std::pair<uint32_t, I3CLSimPhotonSeriesPtr> ConversionResult_t;
        I3CLSimStepToPhotonConverter::ConversionResult_t res =
        openCLStepsToPhotonsConverter_->GetConversionResult();
        if (!res.second) log_fatal("Internal error: received NULL photon series from OpenCL.");

        // convert to I3Photons and add to their respective frames
        AddPhotonsToFrames(*(res.second));
        
        totalNumOutPhotons += res.second->size();
    }
    
    log_info("results fetched from OpenCL.");
    
    log_info("Got %zu photons in total during flush.", totalNumOutPhotons);
    
    log_info("finished.");
    
    for (std::size_t identifier=0;identifier<frameList_.size();++identifier)
    {
        log_info("putting photons into frame %zu...", identifier);
        frameList_[identifier]->Put(photonSeriesMapName_, photonsForFrameList_[identifier]);
        
        log_info("pushing frame number %zu...", identifier);
        PushFrame(frameList_[identifier]);
    }
    
    // empty the frame cache
    particleCache_.clear();
    frameList_.clear();
    photonsForFrameList_.clear();
    currentPhotonIdForFrame_.clear();
}

void I3CLSimModule::Physics(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    if (!geometryIsConfigured_)
        log_fatal("Received Physics frame before Geometry frame");
    
    I3MCTreeConstPtr MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);
    if (!MCTree) log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                           MCTreeName_.c_str());
    
    frameList_.push_back(frame);
    photonsForFrameList_.push_back(I3PhotonSeriesMapPtr(new I3PhotonSeriesMap()));
    currentPhotonIdForFrame_.push_back(0);
    std::size_t currentFrameListIndex = frameList_.size()-1;
    
    for (I3MCTree::iterator it = MCTree->begin();
         it != MCTree->end(); ++it)
    {
        const I3Particle &particle = *it;
        
        // In-ice particles only
        if (particle.GetLocationType() != I3Particle::InIce) continue;

        // check particle type
        const bool isMuon = (particle.GetType() == I3Particle::MuMinus) || (particle.GetType() == I3Particle::MuPlus);
        const bool isNeutrino = particle.IsNeutrino();

        // always ignore neutrinos
        if (isNeutrino) continue;

        // ignore muons if requested
        if ((ignoreMuons_) && (isMuon)) continue;
        
        totalSimulatedEnergyForFlush_ += particle.GetEnergy();
        totalNumParticlesForFlush_++;
        geant4ParticleToStepsConverter_->EnqueueParticle(particle, currentParticleCacheIndex_);

        if (particleCache_.find(currentParticleCacheIndex_) != particleCache_.end())
            log_fatal("Internal error. Particle cache index already used.");
        
        particleCacheEntry &cacheEntry = 
        particleCache_.insert(std::make_pair(currentParticleCacheIndex_, particleCacheEntry())).first->second;
        
        cacheEntry.frameListEntry = currentFrameListIndex;
        cacheEntry.particleMajorID = particle.GetMajorID();
        cacheEntry.particleMinorID = particle.GetMinorID();
        
        // make a new index. This will eventually overflow,
        // but at that time, index 0 should be unused again.
        ++currentParticleCacheIndex_;
    }
    
    if (frameList_.size() >= maxNumParallelEvents_)
    {
        log_info("Flushing results for a total energy of %fGeV for %" PRIu64 " particles",
                 totalSimulatedEnergyForFlush_/I3Units::GeV, totalNumParticlesForFlush_);
                 
        totalSimulatedEnergyForFlush_=0.;
        totalNumParticlesForFlush_=0;
        
        FlushFrameCache();
    }
    
    
}

void I3CLSimModule::Finish()
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    log_info("Flushing results for a total energy of %fGeV for %" PRIu64 " particles",
             totalSimulatedEnergyForFlush_/I3Units::GeV, totalNumParticlesForFlush_);

    totalSimulatedEnergyForFlush_=0.;
    totalNumParticlesForFlush_=0;
    
    FlushFrameCache();
    
    log_info("Flushing I3Tray..");
    Flush();
}


