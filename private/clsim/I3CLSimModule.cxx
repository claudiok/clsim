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

#include "clsim/I3CLSimWlenDependentValueConstant.h"

#include "clsim/I3CLSimParticleToStepConverterGeant4.h"

#include "clsim/I3CLSimModuleHelper.h"

#include <limits>

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

    generateCherenkovPhotonsWithoutDispersion_=true;
    AddParameter("GenerateCherenkovPhotonsWithoutDispersion",
                 "The wavelength of the generated Cherenkov photons will be generated\n"
                 "according to a spectrum without dispersion. This does not change\n"
                 "the total number of photons, only the distribution of wavelengths.",
                 generateCherenkovPhotonsWithoutDispersion_);

    AddParameter("WavelengthGenerationBias",
                 "An instance of I3CLSimWlenDependentValue describing the reciprocal weight a photon gets assigned as a function of its wavelength.\n"
                 "You can set this to the wavelength depended acceptance of your DOM to pre-scale the number of generated photons.",
                 wavelengthGenerationBias_);

    AddParameter("MediumProperties",
                 "An instance of I3CLSimMediumProperties describing the ice/water properties.",
                 mediumProperties_);

    AddParameter("ParameterizationList",
                 "An instance I3CLSimParticleParameterizationSeries specifying the fast simulation parameterizations to be used.",
                 parameterizationList_);

    maxNumParallelEvents_=1000;
    AddParameter("MaxNumParallelEvents",
                 "Maximum number of events that will be processed by the GPU in parallel.",
                 maxNumParallelEvents_);

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

    openCLPlatformName_="";
    AddParameter("OpenCLPlatformName",
                 "Name of the OpenCL platform. Leave empty for auto-selection.",
                 openCLPlatformName_);

    openCLDeviceName_="";
    AddParameter("OpenCLDeviceName",
                 "Name of the OpenCL device. Leave empty for auto-selection.",
                 openCLDeviceName_);

    openCLUseNativeMath_=false;
    AddParameter("OpenCLUseNativeMath",
                 "Use native math instructions in the OpenCL kernel. Has proven not to work\n"
                 "correctly with kernels running on Intel CPUs. Seems to work correctly on\n"
                 "Nvidia GPUs (may speed up things, but make sure it does not change your\n"
                 "results).",
                 openCLUseNativeMath_);

    openCLApproximateNumberOfWorkItems_=512000;
    AddParameter("OpenCLApproximateNumberOfWorkItems",
                 "The approximate number of work items per block. Larger numbers (e.g. 512000)\n"
                 "are ok for dedicated GPGPU cards, but try to keep the number lower if you also\n"
                 "use your GPU for display. Your display may freeze if you use the card interactively\n"
                 "and this number is too high.",
                 openCLApproximateNumberOfWorkItems_);

    DOMRadius_=0.16510*I3Units::m; // 13 inch diameter
    AddParameter("DOMRadius",
                 "The DOM radius used during photon tracking.",
                 DOMRadius_);

    ignoreNonIceCubeOMNumbers_=false;
    AddParameter("IgnoreNonIceCubeOMNumbers",
                 "Ignore string numbers < 1 and OM numbers > 60. (AMANDA and IceTop)",
                 ignoreNonIceCubeOMNumbers_);

    geant4PhysicsListName_="QGSP_BERT";
    AddParameter("Geant4PhysicsListName",
                 "Geant4 physics list name. Examples are \"QGSP_BERT_EMV\" and \"QGSP_BERT\"",
                 geant4PhysicsListName_);

    geant4MaxBetaChangePerStep_=I3CLSimParticleToStepConverterGeant4::default_maxBetaChangePerStep;
    AddParameter("Geant4MaxBetaChangePerStep",
                 "Maximum change of beta=v/c per Geant4 step.",
                 geant4MaxBetaChangePerStep_);

    geant4MaxNumPhotonsPerStep_=I3CLSimParticleToStepConverterGeant4::default_maxNumPhotonsPerStep;
    AddParameter("Geant4MaxNumPhotonsPerStep",
                 "Approximate maximum number of Cherenkov photons generated per step by Geant4.",
                 geant4MaxNumPhotonsPerStep_);

    statisticsName_="";
    AddParameter("StatisticsName",
                 "Collect statistics in this frame object (e.g. number of photons generated or reaching the DOMs)",
                 statisticsName_);


    // add an outbox
    AddOutBox("OutBox");

}

I3CLSimModule::~I3CLSimModule()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    StopThread();
    
}

void I3CLSimModule::StartThread()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    if (threadObj_) {
        log_debug("Thread is already running. Not starting a new one.");
        return;
    }
    
    log_trace("Thread not running. Starting a new one..");

    // clear statistics counters
    photonNumGeneratedPerParticle_.clear();
    photonWeightSumGeneratedPerParticle_.clear();
    
    // re-set flags
    threadStarted_=false;
    threadFinishedOK_=false;
    
    threadObj_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&I3CLSimModule::Thread_starter, this)));
    
    // wait for startup
    {
        boost::unique_lock<boost::mutex> guard(threadStarted_mutex_);
        for (;;)
        {
            if (threadStarted_) break;
            threadStarted_cond_.wait(guard);
        }
    }        
    
    log_trace("Thread started..");

}

void I3CLSimModule::StopThread()
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    if (threadObj_)
    {
        if (threadObj_->joinable())
        {
            log_trace("Stopping the thread..");
            threadObj_->interrupt();
            threadObj_->join(); // wait for it indefinitely
            log_trace("thread stopped.");
        } 
        else
        {
            log_trace("Thread did already finish, deleting reference.");
        }
        threadObj_.reset();
    }
}


void I3CLSimModule::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("RandomService", randomService_);

    GetParameter("GenerateCherenkovPhotonsWithoutDispersion", generateCherenkovPhotonsWithoutDispersion_);
    GetParameter("WavelengthGenerationBias", wavelengthGenerationBias_);

    GetParameter("MediumProperties", mediumProperties_);
    GetParameter("MaxNumParallelEvents", maxNumParallelEvents_);
    GetParameter("MCTreeName", MCTreeName_);
    GetParameter("PhotonSeriesMapName", photonSeriesMapName_);
    GetParameter("IgnoreMuons", ignoreMuons_);
    GetParameter("ParameterizationList", parameterizationList_);

    GetParameter("OpenCLPlatformName", openCLPlatformName_);
    GetParameter("OpenCLDeviceName", openCLDeviceName_);
    GetParameter("OpenCLUseNativeMath", openCLUseNativeMath_);
    GetParameter("OpenCLApproximateNumberOfWorkItems", openCLApproximateNumberOfWorkItems_);

    GetParameter("DOMRadius", DOMRadius_);
    GetParameter("IgnoreNonIceCubeOMNumbers", ignoreNonIceCubeOMNumbers_);

    GetParameter("Geant4PhysicsListName", geant4PhysicsListName_);
    GetParameter("Geant4MaxBetaChangePerStep", geant4MaxBetaChangePerStep_);
    GetParameter("Geant4MaxNumPhotonsPerStep", geant4MaxNumPhotonsPerStep_);

    GetParameter("StatisticsName", statisticsName_);
    collectStatistics_ = (statisticsName_!="");
    
    if (!wavelengthGenerationBias_) {
        wavelengthGenerationBias_ = I3CLSimWlenDependentValueConstantConstPtr(new I3CLSimWlenDependentValueConstant(1.));
    }

    if (!randomService_) log_fatal("You have to specify the \"RandomService\" parameter!");
    if (!mediumProperties_) log_fatal("You have to specify the \"MediumProperties\" parameter!");
    if (maxNumParallelEvents_ <= 0) log_fatal("Values <= 0 are invalid for the \"MaxNumParallelEvents\" parameter!");

    // fill wavelengthGenerator_
    wavelengthGenerator_ =
    I3CLSimModuleHelper::makeWavelengthGenerator
    (wavelengthGenerationBias_,
     generateCherenkovPhotonsWithoutDispersion_,
     mediumProperties_);
    
    currentParticleCacheIndex_ = 0;
    geometryIsConfigured_ = false;
    totalSimulatedEnergyForFlush_ = 0.;
    totalNumParticlesForFlush_ = 0;
    
    if (parameterizationList_.size() > 0) {
        log_info("Using the following parameterizations:");
        
        BOOST_FOREACH(const I3CLSimParticleParameterization &parameterization, parameterizationList_)
        {
            I3Particle tmpParticle;
            tmpParticle.SetType(parameterization.forParticleType);
            
            log_info("  * particle=%s, from energy=%fGeV, to energy=%fGeV",
                     tmpParticle.GetTypeString().c_str(),
                     parameterization.fromEnergy/I3Units::GeV,
                     parameterization.toEnergy/I3Units::GeV);
            
        }
        
    }
    
}


/**
 * This thread takes care of passing steps from Geant4 to OpenCL
 */

bool I3CLSimModule::Thread(boost::this_thread::disable_interruption &di)
{
    // do some setup while the main thread waits..
    numBunchesSentToOpenCL_=0;
    
    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(threadStarted_mutex_);
        threadStarted_=true;
    }
    threadStarted_cond_.notify_all();

    // the main thread is running again
    
    uint32_t counter=0;
    
    for (;;)
    {
        // retrieve steps from Geant4
        I3CLSimStepSeriesConstPtr steps;
        bool barrierWasJustReset=false;
        
        {
            boost::this_thread::restore_interruption ri(di);
            try {
                steps = geant4ParticleToStepsConverter_->GetConversionResultWithBarrierInfo(barrierWasJustReset);
            } catch(boost::thread_interrupted &i) {
                return false;
            }
        }
        
        if (!steps) 
        {
            log_debug("Got NULL I3CLSimStepSeriesConstPtr from Geant4.");
        }
        else if (steps->empty())
        {
            log_warn("Got 0 steps from Geant4, nothing to do for OpenCL.");
        }
        else
        {
            log_warn("Got %zu steps from Geant4, sending them to OpenCL",
                     steps->size());

            
            // collect statistics if requested
            if (collectStatistics_)
            {
                BOOST_FOREACH(const I3CLSimStep &step, *steps)
                {
                    const uint32_t particleID = step.identifier;
                    
                    (photonNumGeneratedPerParticle_.insert(std::make_pair(particleID, 0)).first->second)+=step.numPhotons;
                    (photonWeightSumGeneratedPerParticle_.insert(std::make_pair(particleID, 0.)).first->second)+=static_cast<double>(step.numPhotons)*step.weight;
                }
            }

            // send to OpenCL
            {
                boost::this_thread::restore_interruption ri(di);
                try {
                    openCLStepsToPhotonsConverter_->EnqueueSteps(steps, counter);
                } catch(boost::thread_interrupted &i) {
                    return false;
                }
            }
            
            ++numBunchesSentToOpenCL_;
            ++counter; // this may overflow, but it is not used for anything important/unique
        }
        
        if (barrierWasJustReset) {
            log_trace("Geant4 barrier has been reached. Exiting thread.");
            break;
        }
    }
    
    return true;
}

void I3CLSimModule::Thread_starter()
{
    // do not interrupt this thread by default
    boost::this_thread::disable_interruption di;
    
    try {
        threadFinishedOK_ = Thread(di);
        if (!threadFinishedOK_)
        {
            log_trace("thread exited due to interruption.");
        }
    } catch(...) { // any exceptions?
        log_warn("thread died unexpectedly..");
        
        throw;
    }
    
    log_debug("thread exited.");
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
    if (ignoreNonIceCubeOMNumbers_) 
    {    
        geometry_ = I3CLSimSimpleGeometryFromI3GeometryPtr
        (
         new I3CLSimSimpleGeometryFromI3Geometry(DOMRadius_, geometryObject,
                                                 1,                                     // ignoreStringIDsSmallerThan
                                                 std::numeric_limits<int32_t>::max(),   // ignoreStringIDsLargerThan
                                                 1,                                     // ignoreDomIDsSmallerThan
                                                 60)                                    // ignoreDomIDsLargerThan
        );
    }
    else
    {
        geometry_ = I3CLSimSimpleGeometryFromI3GeometryPtr
        (
         new I3CLSimSimpleGeometryFromI3Geometry(DOMRadius_, geometryObject,
                                                 std::numeric_limits<int32_t>::min(),   // ignoreStringIDsSmallerThan
                                                 std::numeric_limits<int32_t>::max(),   // ignoreStringIDsLargerThan
                                                 std::numeric_limits<uint32_t>::min(),  // ignoreDomIDsSmallerThan
                                                 std::numeric_limits<uint32_t>::max())  // ignoreDomIDsLargerThan
        );
    }
    
    log_info("Initializing CLSim..");
    // initialize OpenCL
    openCLStepsToPhotonsConverter_ =
    I3CLSimModuleHelper::initializeOpenCL(openCLPlatformName_,
                                          openCLDeviceName_,
                                          randomService_,
                                          geometry_,
                                          mediumProperties_,
                                          wavelengthGenerationBias_,
                                          wavelengthGenerator_,
                                          openCLApproximateNumberOfWorkItems_,
                                          openCLUseNativeMath_);
    
    log_info("Initializing Geant4..");
    // initialize Geant4 (will set bunch sizes according to the OpenCL settings)
    geant4ParticleToStepsConverter_ =
    I3CLSimModuleHelper::initializeGeant4(randomService_,
                                          mediumProperties_,
                                          wavelengthGenerationBias_,
                                          openCLStepsToPhotonsConverter_,
                                          parameterizationList_,
                                          geant4PhysicsListName_,
                                          geant4MaxBetaChangePerStep_,
                                          geant4MaxNumPhotonsPerStep_,
                                          false); // the multiprocessor version is not yet safe to use

    
    log_info("Initialization complete.");
    geometryIsConfigured_=true;

    PushFrame(frame);
}

namespace {
    static inline OMKey OMKeyFromOpenCLSimIDs(int16_t stringID, uint16_t domID)
    {
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
        const OMKey key = OMKeyFromOpenCLSimIDs(photon.stringID, photon.omID);
        
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
        outputPhoton.SetWavelength(photon.GetWavelength());
        outputPhoton.SetGroupVelocity(photon.GetGroupVelocity());
        outputPhoton.SetNumScattered(photon.GetNumScatters());

        outputPhoton.SetPos(I3Position(photon.GetPosX(), photon.GetPosY(), photon.GetPosZ()));
        {
            I3Direction outDir;
            outDir.SetThetaPhi(photon.GetDirTheta(), photon.GetDirPhi());
            outputPhoton.SetDir(outDir);
        }

        outputPhoton.SetStartTime(photon.GetStartTime());

        outputPhoton.SetStartPos(I3Position(photon.GetStartPosX(), photon.GetStartPosY(), photon.GetStartPosZ()));
        {
            I3Direction outStartDir;
            outStartDir.SetThetaPhi(photon.GetStartDirTheta(), photon.GetStartDirPhi());
            outputPhoton.SetStartDir(outStartDir);
        }

        if (collectStatistics_)
        {
            // collect statistics
            (photonNumAtOMPerParticle_.insert(std::make_pair(photon.identifier, 0)).first->second)++;
            (photonWeightSumAtOMPerParticle_.insert(std::make_pair(photon.identifier, 0.)).first->second)+=photon.GetWeight();
        }
        
        currentPhotonId++;
    }
    
}

void I3CLSimModule::FlushFrameCache()
{
    log_debug("Flushing frame cache..");
    
    // start the connector thread if necessary
    if (!threadObj_) {
        log_trace("No thread found running during FlushFrameCache(), starting one.");
        StartThread();
    }

    // tell the Geant4 converter to not accept any new data until it is finished 
    // with its current work.
    geant4ParticleToStepsConverter_->EnqueueBarrier();
    
    // At this point the thread should keep on passing steps from Geant4 to OpenCL.
    // As soon as the barrier for Geant4 is no longer active and all data has been
    // sent to OpenCL, the thread will finish. We wait for that and then receive
    // data from the infinite OpenCL queue.
    if (!threadObj_->joinable())
        log_fatal("Thread should be joinable at this point!");
        
    log_debug("Waiting for thread..");
    threadObj_->join(); // wait for it indefinitely
    StopThread(); // stop it completely
    if (!threadFinishedOK_) log_fatal("Thread was aborted or failed.");
    log_debug("thread finished.");

    
    log_info("Geant4 finished, retrieving results from GPU..");

    std::size_t totalNumOutPhotons=0;
    
    photonNumAtOMPerParticle_.clear();
    photonWeightSumAtOMPerParticle_.clear();
    
    for (uint64_t i=0;i<numBunchesSentToOpenCL_;++i)
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
    
    log_warn("Got %zu photons in total during flush.", totalNumOutPhotons);

    if (collectStatistics_)
    {
        std::vector<I3CLSimEventStatisticsPtr> eventStatisticsForFrame;
        for (std::size_t i=0;i<frameList_.size();++i) {
            eventStatisticsForFrame.push_back(I3CLSimEventStatisticsPtr(new I3CLSimEventStatistics()));
        }

        
        // generated photons (count)
        for(std::map<uint32_t, uint64_t>::const_iterator it=photonNumGeneratedPerParticle_.begin();
            it!=photonNumGeneratedPerParticle_.end();++it)
        {
            // find identifier in particle cache
            std::map<uint32_t, particleCacheEntry>::iterator it_cache = particleCache_.find(it->first);
            if (it_cache == particleCache_.end())
                log_fatal("Internal error: unknown particle id from Geant4: %" PRIu32,
                          it_cache->first);
            const particleCacheEntry &cacheEntry = it_cache->second;
            
            if (cacheEntry.frameListEntry >= eventStatisticsForFrame.size())
                log_fatal("Internal error: particle cache entry uses invalid frame cache position");
            
            eventStatisticsForFrame[cacheEntry.frameListEntry]->AddNumPhotonsGeneratedWithWeights(it->second, 0.,
                                                                                                  cacheEntry.particleMajorID,
                                                                                                  cacheEntry.particleMinorID);
        }

        // generated photons (weight sum)
        for(std::map<uint32_t, double>::const_iterator it=photonWeightSumGeneratedPerParticle_.begin();
            it!=photonWeightSumGeneratedPerParticle_.end();++it)
        {
            // find identifier in particle cache
            std::map<uint32_t, particleCacheEntry>::iterator it_cache = particleCache_.find(it->first);
            if (it_cache == particleCache_.end())
                log_fatal("Internal error: unknown particle id from Geant4: %" PRIu32,
                          it_cache->first);
            const particleCacheEntry &cacheEntry = it_cache->second;
            
            if (cacheEntry.frameListEntry >= eventStatisticsForFrame.size())
                log_fatal("Internal error: particle cache entry uses invalid frame cache position");

            eventStatisticsForFrame[cacheEntry.frameListEntry]->AddNumPhotonsGeneratedWithWeights(0, it->second,
                                                                                                  cacheEntry.particleMajorID,
                                                                                                  cacheEntry.particleMinorID);
        }

        
        // photons @ DOMs(count)
        for(std::map<uint32_t, uint64_t>::const_iterator it=photonNumAtOMPerParticle_.begin();
            it!=photonNumAtOMPerParticle_.end();++it)
        {
            // find identifier in particle cache
            std::map<uint32_t, particleCacheEntry>::iterator it_cache = particleCache_.find(it->first);
            if (it_cache == particleCache_.end())
                log_fatal("Internal error: unknown particle id from Geant4: %" PRIu32,
                          it_cache->first);
            const particleCacheEntry &cacheEntry = it_cache->second;
            
            if (cacheEntry.frameListEntry >= eventStatisticsForFrame.size())
                log_fatal("Internal error: particle cache entry uses invalid frame cache position");
            
            eventStatisticsForFrame[cacheEntry.frameListEntry]->AddNumPhotonsAtDOMsWithWeights(it->second, 0.,
                                                                                               cacheEntry.particleMajorID,
                                                                                               cacheEntry.particleMinorID);
        }
        
        // photons @ DOMs (weight sum)
        for(std::map<uint32_t, double>::const_iterator it=photonWeightSumAtOMPerParticle_.begin();
            it!=photonWeightSumAtOMPerParticle_.end();++it)
        {
            // find identifier in particle cache
            std::map<uint32_t, particleCacheEntry>::iterator it_cache = particleCache_.find(it->first);
            if (it_cache == particleCache_.end())
                log_fatal("Internal error: unknown particle id from Geant4: %" PRIu32,
                          it_cache->first);
            const particleCacheEntry &cacheEntry = it_cache->second;
            
            if (cacheEntry.frameListEntry >= eventStatisticsForFrame.size())
                log_fatal("Internal error: particle cache entry uses invalid frame cache position");
            
            eventStatisticsForFrame[cacheEntry.frameListEntry]->AddNumPhotonsAtDOMsWithWeights(0, it->second,
                                                                                               cacheEntry.particleMajorID,
                                                                                               cacheEntry.particleMinorID);
        }

        
        // store statistics to frame
        for (std::size_t i=0;i<frameList_.size();++i)
        {
            frameList_[i]->Put(statisticsName_, eventStatisticsForFrame[i]);
        }
    }    
    
    // not needed anymore
    photonWeightSumGeneratedPerParticle_.clear();
    photonNumGeneratedPerParticle_.clear();
    photonNumAtOMPerParticle_.clear();
    photonWeightSumAtOMPerParticle_.clear();

    
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

namespace {
    bool ParticleHasMuonDaughter(const I3Particle &particle, const I3MCTree &mcTree)
    {
        const std::vector<I3Particle> daughters =
        I3MCTreeUtils::GetDaughters(mcTree, particle);

        BOOST_FOREACH(const I3Particle &daughter, daughters)
        {
            if ((daughter.GetType()==I3Particle::MuMinus) ||
                (daughter.GetType()==I3Particle::MuPlus))
                return true;
        }
        return false;
    }
}

void I3CLSimModule::Physics(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    if (!geometryIsConfigured_)
        log_fatal("Received Physics frame before Geometry frame");

    // start the connector thread if necessary
    if (!threadObj_) {
        log_debug("No thread found running during Physics(), starting one.");
        StartThread();
    }
    
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

        // mmc-icetray currently stores continuous loss entries as "unknown"
        const bool isContinuousLoss = (particle.GetType() == -1111) || (particle.GetType() == I3Particle::unknown); // special magic number used by MMC
        
        // ignore continuous loss entries
        if (isContinuousLoss) {
            log_debug("ignored a continuous loss I3MCTree entry");
            continue;
        }
        
        // always ignore neutrinos
        if (isNeutrino) continue;

        // ignore muons if requested
        if ((ignoreMuons_) && (isMuon)) continue;
        
        // ignore muons with muons as child particles
        // -> those already ran through MMC(-recc) or
        // were sliced with I3MuonSlicer. Only add their
        // children.
        if (!ignoreMuons_) {
            if (ParticleHasMuonDaughter(particle, *MCTree))
                continue;
        }
        
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
        
        log_warn("============== CACHE FLUSHED ================");
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


