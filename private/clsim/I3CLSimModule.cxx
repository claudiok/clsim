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
 * @file I3CLSimModule.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "clsim/I3CLSimModule.h"

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/variant/get.hpp>
#include <boost/math/common_factor_rt.hpp>

#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "simclasses/I3CompressedPhoton.h"

#include "clsim/function/I3CLSimFunctionConstant.h"

#include "clsim/I3CLSimLightSource.h"
#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

#include "clsim/I3CLSimModuleHelper.h"

#include <limits>
#include <set>
#include <deque>
#include <cmath>


namespace {
    class ScopedGILRelease
    {
    public:
        inline ScopedGILRelease()
        {
            m_thread_state = PyEval_SaveThread();
        }
        
        inline ~ScopedGILRelease()
        {
            PyEval_RestoreThread(m_thread_state);
            m_thread_state = NULL;
        }
        
    private:
        PyThreadState *m_thread_state;
    };    
}

// The module
I3_MODULE(I3CLSimModule<I3PhotonSeriesMap>);
I3_MODULE(I3CLSimModule<I3CompressedPhotonSeriesMap>);

template <typename OutputMapType>
I3CLSimModule<OutputMapType>::I3CLSimModule(const I3Context& context) 
: I3ConditionalModule(context),
geometryIsConfigured_(false)
{
    // define parameters
    workOnTheseStops_.clear();
    workOnTheseStops_.push_back(I3Frame::DAQ);
    AddParameter("WorkOnTheseStops",
                 "Work on MCTrees found in the stream types (\"stops\") specified in this list",
                 workOnTheseStops_);

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
                 "An instance of I3CLSimFunction describing the reciprocal weight a photon gets assigned as a function of its wavelength.\n"
                 "You can set this to the wavelength depended acceptance of your DOM to pre-scale the number of generated photons.",
                 wavelengthGenerationBias_);

    AddParameter("MediumProperties",
                 "An instance of I3CLSimMediumProperties describing the ice/water properties.",
                 mediumProperties_);

    AddParameter("SpectrumTable",
                 "All spectra that could be requested by an I3CLSimStep.\n"
                 "If set to NULL/None, only spectrum #0 (Cherenkov photons) will be available.",
                 spectrumTable_);

    AddParameter("ParameterizationList",
                 "An instance I3CLSimLightSourceParameterizationSeries specifying the fast simulation parameterizations to be used.",
                 parameterizationList_);

    maxNumParallelEvents_=1000;
    AddParameter("MaxNumParallelEvents",
                 "Maximum number of events that will be processed by the GPU in parallel.",
                 maxNumParallelEvents_);
    
    totalEnergyToProcess_=0.;
    AddParameter("TotalEnergyToProcess",
                 "Maximum energy that will be processed by the GPU in parallel.",
                 totalEnergyToProcess_);
    
    MCTreeName_="I3MCTree";
    AddParameter("MCTreeName",
                 "Name of the I3MCTree frame object. All particles except neutrinos will be read from this tree.",
                 MCTreeName_);

    flasherPulseSeriesName_="";
    AddParameter("FlasherPulseSeriesName",
                 "Name of the I3CLSimFlasherPulseSeries frame object. Flasher pulses will be read from this object.\n"
                 "Set this to the empty string to disable flashers.",
                 flasherPulseSeriesName_);

    photonSeriesMapName_="PropagatedPhotons";
    AddParameter("PhotonSeriesMapName",
                 "Name of the I3CLSimPhotonSeriesMap frame object that will be written to the frame.",
                 photonSeriesMapName_);

    omKeyMaskName_="";
    AddParameter("OMKeyMaskName",
                 "Name of a I3VectorOMKey or I3VectorModuleKey with masked DOMs. DOMs in this list will not record I3Photons.",
                 omKeyMaskName_);

    ignoreMuons_=false;
    AddParameter("IgnoreMuons",
                 "If set to True, muons will not be propagated.",
                 ignoreMuons_);

    AddParameter("OpenCLDeviceList",
                 "A vector of I3CLSimOpenCLDevice objects, describing the devices to be used for simulation.",
                 openCLDeviceList_);

    DOMRadius_=0.16510*I3Units::m; // 13 inch diameter
    AddParameter("DOMRadius",
                 "The DOM radius used during photon tracking.",
                 DOMRadius_);

    DOMOversizeFactor_=1.; // no oversizing
    AddParameter("DOMOversizeFactor",
                 "Specifiy the \"oversize factor\" (i.e. DOM radius scaling factor).",
                 DOMOversizeFactor_);

    pancakeFactor_=1.; // no oversizing
    AddParameter("DOMPancakeFactor",
                 "Sets the \"pancake\" factor for DOMs. For standard\n"
                 "oversized-DOM simulations, this should be the\n"
                 "radius oversizing factor. This will flatten the\n"
                 "DOM in the direction parallel to the photon.\n"
                 "The DOM will have a pancake-like shape, elongated\n"
                 "in the directions perpendicular to the photon direction.\n"
                 "\n"
                 "The DOM radius (supplied by the geometry) must also include\n"
                 "the oversizing factor. Most of the time you will want to set this\n"
                 "to be the same as the \"DOMOversizeFactor\".",
                 pancakeFactor_);

    ignoreNonIceCubeOMNumbers_=false;
    AddParameter("IgnoreNonIceCubeOMNumbers",
                 "Ignore string numbers < 1 and OM numbers > 60. (AMANDA and IceTop)",
                 ignoreNonIceCubeOMNumbers_);

    geant4PhysicsListName_="QGSP_BERT";
    AddParameter("Geant4PhysicsListName",
                 "Geant4 physics list name. Examples are \"QGSP_BERT_EMV\" and \"QGSP_BERT\"",
                 geant4PhysicsListName_);

    geant4MaxBetaChangePerStep_=I3CLSimLightSourceToStepConverterGeant4::default_maxBetaChangePerStep;
    AddParameter("Geant4MaxBetaChangePerStep",
                 "Maximum change of beta=v/c per Geant4 step.",
                 geant4MaxBetaChangePerStep_);

    geant4MaxNumPhotonsPerStep_=I3CLSimLightSourceToStepConverterGeant4::default_maxNumPhotonsPerStep;
    AddParameter("Geant4MaxNumPhotonsPerStep",
                 "Approximate maximum number of Cherenkov photons generated per step by Geant4.",
                 geant4MaxNumPhotonsPerStep_);

    statisticsName_="";
    AddParameter("StatisticsName",
                 "Collect statistics in this frame object (e.g. number of photons generated or reaching the DOMs)",
                 statisticsName_);

    AddParameter("IgnoreStrings",
                 "Ignore all OMKeys with these string IDs",
                 ignoreStrings_);

    AddParameter("IgnoreDomIDs",
                 "Ignore all OMKeys with these DOM IDs",
                 ignoreDomIDs_);

    AddParameter("IgnoreSubdetectors",
                 "Ignore all OMKeys with these subdetector names",
                 ignoreSubdetectors_);

    splitGeometryIntoPartsAcordingToPosition_=false;
    AddParameter("SplitGeometryIntoPartsAcordingToPosition",
                 "If you have a geometry with multiple OMs per floor (e.g. Antares or KM3NeT-tower-like),\n"
                 "it will internally get split into subdetectors to save memory on the GPU. This is 100% transparent.\n"
                 "By default the split is according to the \"floor index\" of an OM on a floor. If you enable this\n"
                 "option, the split will also be done according to the x-y projected positions of the OMs per string.\n"
                 "This may be necessary for \"tower\" geometries.",
                 splitGeometryIntoPartsAcordingToPosition_);

    useHardcodedDeepCoreSubdetector_=false;
    AddParameter("UseHardcodedDeepCoreSubdetector",
                 "Split off DeepCore as its own two subdetectors (upper and lower part).\n"
                 "This may save constant memory on your GPU.\n"
                 "Assumes that strings [79..86] are DeepCore strings with an upper part at z>-30m\n"
                 "and a lower part at z<-30m.",
                 useHardcodedDeepCoreSubdetector_);

    
    enableDoubleBuffering_=false;
    AddParameter("EnableDoubleBuffering",
                 "Disables or enables double-buffered GPU usage. Double buffering will use\n"
                 "two command queues and two sets of input and output buffers in order to transfer\n"
                 "data to the GPU while a kernel is executing on the other buffer.\n"
                 "This has been observed to yield empty results results on older drivers for the nVidia\n"
                 "architecture, so it is disabled by default.\n"
                 "\n"
                 "Before enabling this for a certain driver/hardware combination, make sure that both correct results\n"
                 "are returned. Most of the time the second buffer results are always empty, so this error should be\n"
                 "easy to observe.",
                 enableDoubleBuffering_);

    doublePrecision_=false;
    AddParameter("DoublePrecision",
                 "Enables double-precision support in the kernel. This slows down calculations and\n"
                 "requires more memory. The performance hit is minimal on CPUs but up to an order\n"
                 "of magnitude on GPUs.",
                 doublePrecision_);

    stopDetectedPhotons_=true;
    AddParameter("StopDetectedPhotons",
                 "Configures behaviour for photons that hit a DOM. If this is true (the default)\n"
                 "photons will be stopped once they hit a DOM. If this is false, they continue to\n"
                 "propagate.",
                 stopDetectedPhotons_);

    saveAllPhotons_=false;
    AddParameter("SaveAllPhotons",
                 "Saves all photons, even if they don't hit a DOM. Cannot be used with \"StopDetectedPhotons\".",
                 saveAllPhotons_);

    saveAllPhotonsPrescale_=0.01;
    AddParameter("SaveAllPhotonsPrescale",
                 "Sets the prescale factor of photons being generated in \"saveAllPhotons\" mode.\n"
                 "Only this fraction of photons is actually generated.",
                 saveAllPhotonsPrescale_);

    fixedNumberOfAbsorptionLengths_=NAN;
    AddParameter("FixedNumberOfAbsorptionLengths",
                 "Sets the number of absorption lengths each photon should be propagated. If set to NaN (the default),\n"
                 "the number is sampled from an exponential distribution. (This is what you want for \"normal\" propagation.)\n"
                 "Use this override for table-making.",
                 fixedNumberOfAbsorptionLengths_);

    photonHistoryEntries_=0;
    AddParameter("PhotonHistoryEntries",
                 "Sets the number of scattering step positions that are saved for a photon hitting\n"
                 "a DOM. The last N photons are saved if there are more scattering points than available entries.",
                 photonHistoryEntries_);

    limitWorkgroupSize_=0;
    AddParameter("LimitWorkgroupSize",
                 "Limits the maximum OpenCL workgroup size (the number of bunches to be processed in parallel).\n"
                 "If set to zero (the default) the largest possible workgroup size will be chosen.",
                 limitWorkgroupSize_);

    closestDOMDistanceCutoff_=300.*I3Units::m;
    AddParameter("ClosestDOMDistanceCutoff",
                 "Do not even start light from sources that do not have any DOMs closer to\n"
                 "to them than this distance.",
                 closestDOMDistanceCutoff_);

    // add an outbox
    AddOutBox("OutBox");

    frameListPhysicsFrameCounter_=0;
}

template <typename OutputMapType>
I3CLSimModule<OutputMapType>::~I3CLSimModule()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    StopThread();
    
}

template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::StartThread()
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
    
    threadObj_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&I3CLSimModule<OutputMapType>::Thread_starter, this)));
    
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

template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::StopThread()
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


template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("WorkOnTheseStops", workOnTheseStops_);
    workOnTheseStops_set_ = std::set<I3Frame::Stream>(workOnTheseStops_.begin(), workOnTheseStops_.end());
    
    GetParameter("RandomService", randomService_);
    if (!randomService_) {
        log_info("Getting the default random service from the context..");
        randomService_ = context_.Get<I3RandomServicePtr>();
        if (!randomService_) 
            log_fatal("You have to specify the \"RandomService\" parameter or add a I3RandomServiceFactor using tray.AddService()!");
    }
    
    GetParameter("GenerateCherenkovPhotonsWithoutDispersion", generateCherenkovPhotonsWithoutDispersion_);
    GetParameter("WavelengthGenerationBias", wavelengthGenerationBias_);

    GetParameter("MediumProperties", mediumProperties_);
    GetParameter("SpectrumTable", spectrumTable_);

    GetParameter("MaxNumParallelEvents", maxNumParallelEvents_);
    GetParameter("TotalEnergyToProcess", totalEnergyToProcess_);
    GetParameter("MCTreeName", MCTreeName_);
    GetParameter("FlasherPulseSeriesName", flasherPulseSeriesName_);
    GetParameter("PhotonSeriesMapName", photonSeriesMapName_);
    GetParameter("OMKeyMaskName", omKeyMaskName_);
    GetParameter("IgnoreMuons", ignoreMuons_);
    GetParameter("ParameterizationList", parameterizationList_);

    GetParameter("OpenCLDeviceList", openCLDeviceList_);

    GetParameter("DOMRadius", DOMRadius_);
    GetParameter("DOMOversizeFactor", DOMOversizeFactor_);
    GetParameter("DOMPancakeFactor", pancakeFactor_);

    GetParameter("IgnoreNonIceCubeOMNumbers", ignoreNonIceCubeOMNumbers_);

    GetParameter("Geant4PhysicsListName", geant4PhysicsListName_);
    GetParameter("Geant4MaxBetaChangePerStep", geant4MaxBetaChangePerStep_);
    GetParameter("Geant4MaxNumPhotonsPerStep", geant4MaxNumPhotonsPerStep_);

    GetParameter("StatisticsName", statisticsName_);
    collectStatistics_ = (statisticsName_!="");
    
    GetParameter("IgnoreStrings", ignoreStrings_);
    GetParameter("IgnoreDomIDs", ignoreDomIDs_);
    GetParameter("IgnoreSubdetectors", ignoreSubdetectors_);

    GetParameter("SplitGeometryIntoPartsAcordingToPosition", splitGeometryIntoPartsAcordingToPosition_);

    GetParameter("UseHardcodedDeepCoreSubdetector", useHardcodedDeepCoreSubdetector_);

    GetParameter("EnableDoubleBuffering", enableDoubleBuffering_);
    GetParameter("DoublePrecision", doublePrecision_);
    GetParameter("StopDetectedPhotons", stopDetectedPhotons_);
    GetParameter("SaveAllPhotons", saveAllPhotons_);
    GetParameter("SaveAllPhotonsPrescale", saveAllPhotonsPrescale_);

    GetParameter("FixedNumberOfAbsorptionLengths", fixedNumberOfAbsorptionLengths_);

    GetParameter("PhotonHistoryEntries", photonHistoryEntries_);

    GetParameter("LimitWorkgroupSize", limitWorkgroupSize_);

    GetParameter("ClosestDOMDistanceCutoff", closestDOMDistanceCutoff_);

    if (pancakeFactor_ != DOMOversizeFactor_) {
        log_warn("***** You set the \"DOMOversizeFactor\" to a different value than the \"DOMPancakeFactor\". Be sure you know what you are doing!");
    }
    
    if ((saveAllPhotons_) && (stopDetectedPhotons_)) {
        log_fatal("The \"SaveAllPhotons\" option cannot be used when \"StopDetectedPhotons\" is active.");
    }
    
    if ((flasherPulseSeriesName_=="") && (MCTreeName_==""))
        log_fatal("You need to set at least one of the \"MCTreeName\" and \"FlasherPulseSeriesName\" parameters.");
    
    if (!wavelengthGenerationBias_) {
        wavelengthGenerationBias_ = I3CLSimFunctionConstantConstPtr(new I3CLSimFunctionConstant(1.));
    }

    if (!mediumProperties_) log_fatal("You have to specify the \"MediumProperties\" parameter!");

    if ((totalEnergyToProcess_ > 0) && (!std::isnan(totalEnergyToProcess_)))    
    {
        log_warn("Total Energy to Process mode! MaxNumParallelEvents is set to 1! "
                 "CLSim is going to figure out the number of frames to process "
                 "from the light deposited in each frame!");
        maxNumParallelEvents_=1;
        maxNumParallelEventsSecondFlush_=1;
    }
    if ((maxNumParallelEvents_ <= 0) && (totalEnergyToProcess_ <= 0)) 
        log_fatal("Values <= 0 are invalid for both the \"MaxNumParallelEvents\" and \"TotalEnergyToProcess\" parameter!");
    // maxNumParallelEvents_ is the number of frames buffered by this module.
    // Since we use double-buffering, divide the number by 2.
    maxNumParallelEvents_ /= 2;
    if (maxNumParallelEvents_==0) maxNumParallelEvents_=1;
    maxNumParallelEventsSecondFlush_ = maxNumParallelEvents_;
    

    if (openCLDeviceList_.empty()) 
        log_fatal("You have to provide at least one OpenCL device using the \"OpenCLDeviceList\" parameter.");
    
    // fill wavelengthGenerators_[0] (index 0 is the Cherenkov generator)
    wavelengthGenerators_.clear();
    wavelengthGenerators_.push_back(I3CLSimModuleHelper::makeCherenkovWavelengthGenerator
                                    (wavelengthGenerationBias_,
                                     generateCherenkovPhotonsWithoutDispersion_,
                                     mediumProperties_
                                    )
                                   );
    
    if ((spectrumTable_) && (spectrumTable_->size() > 1)) {
        // a spectrum table has been configured and it contains more than the
        // default Cherenkov spectrum at index #0.
        
        for (std::size_t i=1;i<spectrumTable_->size();++i)
        {
            wavelengthGenerators_.push_back(I3CLSimModuleHelper::makeWavelengthGenerator
                                            ((*spectrumTable_)[i],
                                             wavelengthGenerationBias_,
                                             mediumProperties_
                                             )
                                            );
        }
        
        log_info("%zu additional (non-Cherenkov) wavelength generators (spectra) have been configured.",
                 spectrumTable_->size()-1);
    }

    currentParticleCacheIndex_ = 1;
    geometryIsConfigured_ = false;
    totalSimulatedEnergyForFlush_ = 0.;
    totalSimulatedEnergy_ = 0;
    totalSimulatedEnergySecondFlush_ = 0;
    totalNumParticlesForFlush_ = 0;
    
    if (parameterizationList_.size() > 0) {
        log_info("Using the following parameterizations:");
        
        BOOST_FOREACH(const I3CLSimLightSourceParameterization &parameterization, parameterizationList_)
        {
            I3Particle tmpParticle;
#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
            tmpParticle.SetType(parameterization.forParticleType);
#else
            tmpParticle.SetPdgEncoding(parameterization.forPdgEncoding);
#endif

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

template <typename OutputMapType>
bool I3CLSimModule<OutputMapType>::Thread(boost::this_thread::disable_interruption &di)
{
    // do some setup while the main thread waits..
    numBunchesSentToOpenCL_.assign(openCLStepsToPhotonsConverters_.size(), 0);
    
    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(threadStarted_mutex_);
        threadStarted_=true;
    }
    threadStarted_cond_.notify_all();

    // the main thread is running again
    
    uint32_t counter=0;
    std::size_t lastDeviceIndexToUse=0;
    
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
            log_debug("Got 0 steps from Geant4, nothing to do for OpenCL.");
        }
        else
        {
            log_debug("Got %zu steps from Geant4, sending them to OpenCL",
                     steps->size());

            
            // collect statistics if requested
            if (collectStatistics_)
            {
                BOOST_FOREACH(const I3CLSimStep &step, *steps)
                {
                    const uint32_t particleID = step.identifier;

                    // skip dummy steps
                    if ((step.weight<=0.) || (step.numPhotons<=0)) continue;
                    
                    // sanity check
                    if (particleID==0) log_fatal("particleID==0, this should not happen (this index is never used)");
                    
                    (photonNumGeneratedPerParticle_.insert(std::make_pair(particleID, 0)).first->second)+=step.numPhotons;
                    (photonWeightSumGeneratedPerParticle_.insert(std::make_pair(particleID, 0.)).first->second)+=static_cast<double>(step.numPhotons)*step.weight;
                }
            }

            // determine which OpenCL device to use
            std::vector<std::size_t> fillLevels(openCLStepsToPhotonsConverters_.size());
            for (std::size_t i=0;i<openCLStepsToPhotonsConverters_.size();++i)
            {
                fillLevels[i]=openCLStepsToPhotonsConverters_[i]->QueueSize();
            }
            
            std::size_t minimumFillLevel = fillLevels[0];
            for (std::size_t i=1;i<openCLStepsToPhotonsConverters_.size();++i)
            {
                if (fillLevels[i] < minimumFillLevel)
                {
                    minimumFillLevel=fillLevels[i];
                }
            }
            
            std::size_t deviceIndexToUse=lastDeviceIndexToUse;
            do {
                ++deviceIndexToUse;
                if (deviceIndexToUse>=fillLevels.size()) deviceIndexToUse=0;
            } while (fillLevels[deviceIndexToUse] != minimumFillLevel);
            lastDeviceIndexToUse=deviceIndexToUse;
            
            // send to OpenCL
            {
                boost::this_thread::restore_interruption ri(di);
                try {
                    openCLStepsToPhotonsConverters_[deviceIndexToUse]->EnqueueSteps(steps, counter);
                } catch(boost::thread_interrupted &i) {
                    return false;
                }
            }
            
            ++numBunchesSentToOpenCL_[deviceIndexToUse];
            ++counter; // this may overflow, but it is not used for anything important/unique
        }
        
        if (barrierWasJustReset) {
            log_trace("Geant4 barrier has been reached. Exiting thread.");
            break;
        }
    }
    
    return true;
}

template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::Thread_starter()
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


template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::DigestGeometry(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    if (geometryIsConfigured_)
        log_fatal("This module currently supports only a single geometry per input file.");
    
    //log_debug("Retrieving geometry..");
    //I3GeometryConstPtr geometryObject = frame->Get<I3GeometryConstPtr>();
    //if (!geometryObject) log_fatal("Geometry frame does not have an I3Geometry object!");
    
    log_debug("Converting geometry..");
    
    std::set<int> ignoreStringsSet(ignoreStrings_.begin(), ignoreStrings_.end());
    std::set<unsigned int> ignoreDomIDsSet(ignoreDomIDs_.begin(), ignoreDomIDs_.end());
    std::set<std::string> ignoreSubdetectorsSet(ignoreSubdetectors_.begin(), ignoreSubdetectors_.end());
    
    if (ignoreNonIceCubeOMNumbers_) 
    {    
        geometry_ = I3CLSimSimpleGeometryFromI3GeometryPtr
        (
         new I3CLSimSimpleGeometryFromI3Geometry(DOMRadius_,
                                                 DOMOversizeFactor_,
                                                 frame,
                                                 ignoreStringsSet,
                                                 ignoreDomIDsSet,
                                                 ignoreSubdetectorsSet,
                                                 1,                                     // ignoreStringIDsSmallerThan
                                                 std::numeric_limits<int32_t>::max(),   // ignoreStringIDsLargerThan
                                                 1,                                     // ignoreDomIDsSmallerThan
                                                 60,                                    // ignoreDomIDsLargerThan
                                                 splitGeometryIntoPartsAcordingToPosition_,
                                                 useHardcodedDeepCoreSubdetector_)
        );
    }
    else
    {
        geometry_ = I3CLSimSimpleGeometryFromI3GeometryPtr
        (
         new I3CLSimSimpleGeometryFromI3Geometry(DOMRadius_,
                                                 DOMOversizeFactor_,
                                                 frame,
                                                 ignoreStringsSet,
                                                 ignoreDomIDsSet,
                                                 ignoreSubdetectorsSet,
                                                 std::numeric_limits<int32_t>::min(),   // ignoreStringIDsSmallerThan
                                                 std::numeric_limits<int32_t>::max(),   // ignoreStringIDsLargerThan
                                                 std::numeric_limits<uint32_t>::min(),  // ignoreDomIDsSmallerThan
                                                 std::numeric_limits<uint32_t>::max(),  // ignoreDomIDsLargerThan
                                                 splitGeometryIntoPartsAcordingToPosition_,
                                                 useHardcodedDeepCoreSubdetector_)
        );
    }
    
    if (closestDOMDistanceCutoff_ >= 0 and std::isfinite(closestDOMDistanceCutoff_)) {
        std::vector<I3Position> doms;
        const std::vector<double>& x = geometry_->GetPosXVector();
        const std::vector<double>& y = geometry_->GetPosYVector();
        const std::vector<double>& z = geometry_->GetPosZVector();
        
        for (std::size_t i=0;i<geometry_->size();++i)
            doms.emplace_back(I3Position(x[i],y[i],z[i]));
        detectorHull_.reset(new I3Surfaces::ExtrudedPolygon(doms, closestDOMDistanceCutoff_));
    }
    
    log_info("Initializing CLSim..");
    // initialize OpenCL converters
    openCLStepsToPhotonsConverters_.clear();
    
    uint64_t granularity=0;
    uint64_t maxBunchSize=0;
    
    BOOST_FOREACH(const I3CLSimOpenCLDevice &openCLdevice, openCLDeviceList_)
    {
#ifdef I3_LOG4CPLUS_LOGGING
        LOG_IMPL(INFO, " -> platform: %s device: %s",
                 openCLdevice.GetPlatformName().c_str(), openCLdevice.GetDeviceName().c_str());
#else
        log_info(" -> platform: %s device: %s",
                 openCLdevice.GetPlatformName().c_str(), openCLdevice.GetDeviceName().c_str());
#endif
        
        I3CLSimStepToPhotonConverterOpenCLPtr openCLStepsToPhotonsConverter =
        I3CLSimModuleHelper::initializeOpenCL(openCLdevice,
                                              randomService_,
                                              geometry_,
                                              mediumProperties_,
                                              wavelengthGenerationBias_,
                                              wavelengthGenerators_,
                                              enableDoubleBuffering_,
                                              doublePrecision_,
                                              stopDetectedPhotons_,
                                              saveAllPhotons_,
                                              saveAllPhotonsPrescale_,
                                              fixedNumberOfAbsorptionLengths_,
                                              pancakeFactor_,
                                              photonHistoryEntries_,
                                              limitWorkgroupSize_);
        if (!openCLStepsToPhotonsConverter)
            log_fatal("Could not initialize OpenCL!");
        
        if (openCLStepsToPhotonsConverter->GetWorkgroupSize()==0)
            log_fatal("Internal error: converter.GetWorkgroupSize()==0.");
        if (openCLStepsToPhotonsConverter->GetMaxNumWorkitems()==0)
            log_fatal("Internal error: converter.GetMaxNumWorkitems()==0.");
        
        openCLStepsToPhotonsConverters_.push_back(openCLStepsToPhotonsConverter);
        
        if (granularity==0) {
            granularity = openCLStepsToPhotonsConverter->GetWorkgroupSize();
        } else {
            // least common multiple
            const uint64_t currentGranularity = openCLStepsToPhotonsConverter->GetWorkgroupSize();
            const uint64_t newGranularity = boost::math::lcm(currentGranularity, granularity);
            
            if (newGranularity != granularity) {
#ifdef I3_LOG4CPLUS_LOGGING
                LOG_IMPL(INFO, "new OpenCL device work group size is not compatible (%" PRIu64 "), changing granularity from %" PRIu64 " to %" PRIu64,
                         currentGranularity, granularity, newGranularity);
#else
                log_info("new OpenCL device work group size is not compatible (%" PRIu64 "), changing granularity from %" PRIu64 " to %" PRIu64,
                         currentGranularity, granularity, newGranularity);
#endif
            }
            
            granularity=newGranularity;
        }

        if (maxBunchSize==0) {
            maxBunchSize = openCLStepsToPhotonsConverter->GetMaxNumWorkitems();
        } else {
            const uint64_t currentMaxBunchSize = openCLStepsToPhotonsConverter->GetMaxNumWorkitems();
            const uint64_t newMaxBunchSize = std::min(maxBunchSize, currentMaxBunchSize);
            const uint64_t newMaxBunchSizeWithGranularity = newMaxBunchSize - newMaxBunchSize%granularity;

            if (newMaxBunchSizeWithGranularity != maxBunchSize)
            {
#ifdef I3_LOG4CPLUS_LOGGING
                LOG_IMPL(INFO, "maximum bunch size decreased from %" PRIu64 " to %" PRIu64 " because of new devices maximum request of %" PRIu64 " and a granularity of %" PRIu64,
                         maxBunchSize, newMaxBunchSizeWithGranularity, currentMaxBunchSize, granularity);
#else
                log_info("maximum bunch size decreased from %" PRIu64 " to %" PRIu64 " because of new devices maximum request of %" PRIu64 " and a granularity of %" PRIu64,
                         maxBunchSize, newMaxBunchSizeWithGranularity, currentMaxBunchSize, granularity);
#endif
            }

            if (newMaxBunchSizeWithGranularity==0)
                log_fatal("maximum bunch sizes are incompatible with kernel work group sizes.");
            
            maxBunchSize = newMaxBunchSizeWithGranularity;
        }
        
    }
    
    
    log_info("Initializing Geant4..");
    // initialize Geant4 (will set bunch sizes according to the OpenCL settings)
    geant4ParticleToStepsConverter_ =
    I3CLSimModuleHelper::initializeGeant4(randomService_,
                                          mediumProperties_,
                                          wavelengthGenerationBias_,
                                          granularity,
                                          maxBunchSize,
                                          parameterizationList_,
                                          geant4PhysicsListName_,
                                          geant4MaxBetaChangePerStep_,
                                          geant4MaxNumPhotonsPerStep_,
                                          false); // the multiprocessor version is not yet safe to use

    
    log_info("Initialization complete.");
    geometryIsConfigured_=true;
}

namespace {
    static inline ModuleKey ModuleKeyFromOpenCLSimIDs(int16_t stringID, uint16_t domID)
    {
        return ModuleKey(stringID, domID);
    }
}

namespace {

typedef typename I3CLSimModule<I3PhotonSeriesMap>::particleCacheEntry particleCacheEntry;

void EmitPhoton(const I3CLSimPhoton &photon, uint32_t currentPhotonId,
    double timeShift, uint32_t particleMinorID, uint64_t particleMajorID,
    I3PhotonSeries &outputPhotonSeries)
{
    // append a new I3Photon to the list
    outputPhotonSeries.push_back(I3Photon());
    
    // get a reference to the new photon
    I3Photon &outputPhoton = outputPhotonSeries.back();

    // fill the photon data
    outputPhoton.SetTime(photon.GetTime() + timeShift);
    outputPhoton.SetID(currentPhotonId); // per-frame ID for every photon
    outputPhoton.SetWeight(photon.GetWeight());
    outputPhoton.SetParticleMinorID(particleMinorID);
    outputPhoton.SetParticleMajorID(particleMajorID);
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

    outputPhoton.SetStartTime(photon.GetStartTime() + timeShift);

    outputPhoton.SetStartPos(I3Position(photon.GetStartPosX(), photon.GetStartPosY(), photon.GetStartPosZ()));
    {
        I3Direction outStartDir;
        outStartDir.SetThetaPhi(photon.GetStartDirTheta(), photon.GetStartDirPhi());
        outputPhoton.SetStartDir(outStartDir);
    }

    outputPhoton.SetDistanceInAbsorptionLengths(photon.GetDistInAbsLens());
}

void EmitPhoton(const I3CLSimPhoton &photon, uint32_t currentPhotonId,
    double timeShift, uint32_t particleMinorID, uint64_t particleMajorID,
    I3CompressedPhotonSeries &outputPhotonSeries)
{
    // append a new I3Photon to the list
    outputPhotonSeries.push_back(I3Photon());
    
    // get a reference to the new photon
    I3CompressedPhoton &outputPhoton = outputPhotonSeries.back();

    // fill the photon data
    outputPhoton.SetTime(photon.GetTime() + timeShift);
    outputPhoton.SetWeight(photon.GetWeight());
    outputPhoton.SetParticleMinorID(particleMinorID);
    outputPhoton.SetParticleMajorID(particleMajorID);
    outputPhoton.SetWavelength(photon.GetWavelength());
    outputPhoton.SetGroupVelocity(photon.GetGroupVelocity());

    outputPhoton.SetPos(I3Position(photon.GetPosX(), photon.GetPosY(), photon.GetPosZ()));
    {
        I3Direction outDir;
        outDir.SetThetaPhi(photon.GetDirTheta(), photon.GetDirPhi());
        outputPhoton.SetDir(outDir);
    }
}


void AddHistoryEntries(const I3CLSimPhotonHistory &photonHistory, I3PhotonSeries &outputPhotonSeries)
{
    // get a reference to the new photon
    I3Photon &outputPhoton = outputPhotonSeries.back();
    
    if (photonHistory.size() > outputPhoton.GetNumScattered())
        log_fatal("Logic error: photonHistory.size() [==%zu] > photon.GetNumScatters() [==%zu]",
                  photonHistory.size(), static_cast<std::size_t>(outputPhoton.GetNumScattered()));
    
    for (std::size_t j=0;j<photonHistory.size();++j)
    {
        outputPhoton.AppendToIntermediatePositionList(I3Position( photonHistory.GetX(j), photonHistory.GetY(j), photonHistory.GetZ(j) ),
                                                      photonHistory.GetDistanceInAbsorptionLengths(j)
                                                     );
    }
}

void AddHistoryEntries(const I3CLSimPhotonHistory &photonHistory, I3CompressedPhotonSeries &outputPhotonSeries)
{
    
}

template <typename OutputMapType>
void AddPhotonsToFrames(const I3CLSimPhotonSeries &photons,
                        I3CLSimPhotonHistorySeriesConstPtr photonHistories,
                        const std::vector<boost::shared_ptr<OutputMapType> > &photonsForFrameList_,
                        std::vector<int32_t> &currentPhotonIdForFrame_,
                        const std::vector<I3FramePtr> &frameList_,
                        const std::map<uint32_t, typename I3CLSimModule<OutputMapType>::particleCacheEntry> &particleCache_,
                        const std::vector<std::set<ModuleKey> > &maskedOMKeys_,
                        bool collectStatistics_,
                        std::map<uint32_t, uint64_t> &photonNumAtOMPerParticle,
                        std::map<uint32_t, double> &photonWeightSumAtOMPerParticle
                        )
{
    if (photonsForFrameList_.size() != frameList_.size())
        log_fatal("Internal error: cache sizes differ. (1)");
    if (photonsForFrameList_.size() != currentPhotonIdForFrame_.size())
        log_fatal("Internal error: cache sizes differ. (2)");
    
    if (photonHistories) {
        if (photonHistories->size() != photons.size())
        {
            log_fatal("Error: photon history vector size (%zu) != photon vector size (%zu)",
                      photonHistories->size(), photons.size());
        }
    }
    
    typedef OutputMapType PhotonSeriesMap;
    typedef typename PhotonSeriesMap::mapped_type PhotonSeries;
    typedef typename PhotonSeries::value_type Photon;
    typedef typename I3CLSimModule<OutputMapType>::particleCacheEntry particleCacheEntry;
    for (std::size_t i=0;i<photons.size();++i)
    {
        const I3CLSimPhoton &photon = photons[i];
        
        // find identifier in particle cache
        auto it = particleCache_.find(photon.identifier);
        if (it == particleCache_.end())
            log_fatal("Internal error: unknown particle id from OpenCL: %" PRIu32,
                      photon.identifier);
        const particleCacheEntry &cacheEntry = it->second;

        if (cacheEntry.frameListEntry >= photonsForFrameList_.size())
            log_fatal("Internal error: particle cache entry uses invalid frame cache position");
        
        //I3FramePtr &frame = frameList_[cacheEntry.frameListEntry];
        PhotonSeriesMap &outputPhotonMap = *(photonsForFrameList_[cacheEntry.frameListEntry]);

        // get the current photon id
        int32_t &currentPhotonId = currentPhotonIdForFrame_[cacheEntry.frameListEntry];
        
        // generate the OMKey
        const ModuleKey key = ModuleKeyFromOpenCLSimIDs(photon.stringID, photon.omID);
        
        // get the OMKey mask
        const std::set<ModuleKey> &keyMask = maskedOMKeys_[cacheEntry.frameListEntry];

        if (keyMask.count(key) > 0) continue; // ignore masked DOMs
        
        // this either inserts a new vector or retrieves an existing one
        PhotonSeries &outputPhotonSeries = outputPhotonMap.insert(std::make_pair(key, PhotonSeries())).first->second;
        
        EmitPhoton(photon, currentPhotonId, cacheEntry.timeShift,
            cacheEntry.particleMinorID, cacheEntry.particleMajorID,
            outputPhotonSeries);
        
        if (photonHistories) {
            const I3CLSimPhotonHistory &photonHistory = (*photonHistories)[i];
            AddHistoryEntries(photonHistory, outputPhotonSeries);
        }
        
        if (collectStatistics_)
        {
            // collect statistics
            (photonNumAtOMPerParticle.insert(std::make_pair(photon.identifier, 0)).first->second)++;
            (photonWeightSumAtOMPerParticle.insert(std::make_pair(photon.identifier, 0.)).first->second)+=photon.GetWeight();
        }
        
        currentPhotonId++;
    }
    
}

} // namespace

template <typename OutputMapType>
std::size_t I3CLSimModule<OutputMapType>::FlushFrameCache()
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
        
    {
        // allow other threads to access python
        ScopedGILRelease scopedGIL;

        log_debug("Waiting for thread..");
        threadObj_->join(); // wait for it indefinitely
        StopThread(); // stop it completely
        if (!threadFinishedOK_) log_fatal("Thread was aborted or failed.");
        log_debug("thread finished.");
    }

    // swap all frame cache objects with local versions

    std::map<uint32_t, uint64_t> photonNumGeneratedPerParticle_old;
    std::map<uint32_t, double> photonWeightSumGeneratedPerParticle_old;
    photonNumGeneratedPerParticle_old.swap(photonNumGeneratedPerParticle_);
    photonWeightSumGeneratedPerParticle_old.swap(photonWeightSumGeneratedPerParticle_);

    std::vector<boost::shared_ptr<OutputMapType> > photonsForFrameList_old;
    std::vector<int32_t> currentPhotonIdForFrame_old;
    std::vector<I3FramePtr> frameList_old;
    std::map<uint32_t, particleCacheEntry> particleCache_old;
    std::vector<std::set<ModuleKey> > maskedOMKeys_old;
    std::vector<bool> frameIsBeingWorkedOn_old;

    photonsForFrameList_old.swap(photonsForFrameList_);
    currentPhotonIdForFrame_old.swap(currentPhotonIdForFrame_);
    frameList_old.swap(frameList_);
    particleCache_old.swap(particleCache_);
    maskedOMKeys_old.swap(maskedOMKeys_);
    frameIsBeingWorkedOn_old.swap(frameIsBeingWorkedOn_);

    bool startThreadLater = false;

    // at this point, if we have frames in the secondary cache, 
    // immediately push them to Geant4 so it can start working on them
    for (;;)
    {
        if (frameList2_.empty()) break;
        if (frameList_.size() >= maxNumParallelEvents_) break;

        DigestOtherFrame(frameList2_.front(), false); // do not start the connector thread right now
        startThreadLater=true;
        frameList2_.pop_front();
    }

    // now wait for OpenCL to finish; retrieve results
    std::map<uint32_t, uint64_t> photonNumAtOMPerParticle;
    std::map<uint32_t, double> photonWeightSumAtOMPerParticle;

    std::deque<I3CLSimStepToPhotonConverter::ConversionResult_t> res_list;
    for (std::size_t deviceIndex=0;deviceIndex<numBunchesSentToOpenCL_.size();++deviceIndex)
    {
        log_debug("Geant4 finished, retrieving results from GPU %zu..", deviceIndex);

        for (uint64_t i=0;i<numBunchesSentToOpenCL_[deviceIndex];++i)
        {
            I3CLSimStepToPhotonConverter::ConversionResult_t res =
            openCLStepsToPhotonsConverters_[deviceIndex]->GetConversionResult();
            if (!res.photons) log_fatal("Internal error: received NULL photon series from OpenCL.");

            res_list.push_back(res);
        }
    }
    
    log_debug("results fetched from OpenCL.");

    // new frames were already sent to Geant4, we can re-start the thread right now
    // since we are done with fetching results from OpenCL    
    if (startThreadLater) {
        // start the connector thread if necessary
        if (!threadObj_) {
            log_debug("Delayed start of Thread().");
            StartThread();
        }
    }

    log_debug("Adding photons to frame.");
    std::size_t totalNumOutPhotons=0;

    while (!res_list.empty()) 
    {
        const I3CLSimStepToPhotonConverter::ConversionResult_t &res =
            res_list.front();

        // convert to I3Photons and add to their respective frames
        AddPhotonsToFrames(*(res.photons), res.photonHistories,
                           photonsForFrameList_old,
                           currentPhotonIdForFrame_old,
                           frameList_old,
                           particleCache_old,
                           maskedOMKeys_old,
                           collectStatistics_,
                           photonNumAtOMPerParticle,
                           photonWeightSumAtOMPerParticle
                           );
        
        totalNumOutPhotons += res.photons->size();

        res_list.pop_front();
    }


    log_debug("Got %zu photons in total during flush.", totalNumOutPhotons);

    if (collectStatistics_)
    {
        std::vector<I3CLSimEventStatisticsPtr> eventStatisticsForFrame;
        for (std::size_t i=0;i<frameList_old.size();++i) {
            if (frameIsBeingWorkedOn_old[i]) {
                eventStatisticsForFrame.push_back(I3CLSimEventStatisticsPtr(new I3CLSimEventStatistics()));
            } else {
                eventStatisticsForFrame.push_back(I3CLSimEventStatisticsPtr()); // NULL pointer for non-physics(/DAQ)-frames
            }
        }

        
        // generated photons (count)
        for(std::map<uint32_t, uint64_t>::const_iterator it=photonNumGeneratedPerParticle_old.begin();
            it!=photonNumGeneratedPerParticle_old.end();++it)
        {
            // find identifier in particle cache
            typename std::map<uint32_t, particleCacheEntry>::iterator it_cache = particleCache_old.find(it->first);
            if (it_cache == particleCache_old.end())
                log_error("Internal error: unknown particle id from Geant4: %" PRIu32,
                          it->first);
            const particleCacheEntry &cacheEntry = it_cache->second;
            
            if (cacheEntry.frameListEntry >= eventStatisticsForFrame.size())
                log_fatal("Internal error: particle cache entry uses invalid frame cache position");
            
            eventStatisticsForFrame[cacheEntry.frameListEntry]->AddNumPhotonsGeneratedWithWeights(it->second, 0.,
                                                                                                  cacheEntry.particleMajorID,
                                                                                                  cacheEntry.particleMinorID);
        }

        // generated photons (weight sum)
        for(std::map<uint32_t, double>::const_iterator it=photonWeightSumGeneratedPerParticle_old.begin();
            it!=photonWeightSumGeneratedPerParticle_old.end();++it)
        {
            // find identifier in particle cache
            typename std::map<uint32_t, particleCacheEntry>::iterator it_cache = particleCache_old.find(it->first);
            if (it_cache == particleCache_old.end())
                log_fatal("Internal error: unknown particle id from Geant4: %" PRIu32,
                          it->first);
            const particleCacheEntry &cacheEntry = it_cache->second;
            
            if (cacheEntry.frameListEntry >= eventStatisticsForFrame.size())
                log_fatal("Internal error: particle cache entry uses invalid frame cache position");

            eventStatisticsForFrame[cacheEntry.frameListEntry]->AddNumPhotonsGeneratedWithWeights(0, it->second,
                                                                                                  cacheEntry.particleMajorID,
                                                                                                  cacheEntry.particleMinorID);
        }

        
        // photons @ DOMs(count)
        for(std::map<uint32_t, uint64_t>::const_iterator it=photonNumAtOMPerParticle.begin();
            it!=photonNumAtOMPerParticle.end();++it)
        {
            // find identifier in particle cache
            typename std::map<uint32_t, particleCacheEntry>::iterator it_cache = particleCache_old.find(it->first);
            if (it_cache == particleCache_old.end())
                log_fatal("Internal error: unknown particle id from Geant4: %" PRIu32,
                          it->first);
            const particleCacheEntry &cacheEntry = it_cache->second;
            
            if (cacheEntry.frameListEntry >= eventStatisticsForFrame.size())
                log_fatal("Internal error: particle cache entry uses invalid frame cache position");
            
            eventStatisticsForFrame[cacheEntry.frameListEntry]->AddNumPhotonsAtDOMsWithWeights(it->second, 0.,
                                                                                               cacheEntry.particleMajorID,
                                                                                               cacheEntry.particleMinorID);
        }
        
        // photons @ DOMs (weight sum)
        for(std::map<uint32_t, double>::const_iterator it=photonWeightSumAtOMPerParticle.begin();
            it!=photonWeightSumAtOMPerParticle.end();++it)
        {
            // find identifier in particle cache
            typename std::map<uint32_t, particleCacheEntry>::iterator it_cache = particleCache_old.find(it->first);
            if (it_cache == particleCache_old.end())
                log_fatal("Internal error: unknown particle id from Geant4: %" PRIu32,
                          it->first);
            const particleCacheEntry &cacheEntry = it_cache->second;
            
            if (cacheEntry.frameListEntry >= eventStatisticsForFrame.size())
                log_fatal("Internal error: particle cache entry uses invalid frame cache position");
            
            eventStatisticsForFrame[cacheEntry.frameListEntry]->AddNumPhotonsAtDOMsWithWeights(0, it->second,
                                                                                               cacheEntry.particleMajorID,
                                                                                               cacheEntry.particleMinorID);
        }

        
        // store statistics to frame
        for (std::size_t i=0;i<frameList_old.size();++i)
        {
            if (frameIsBeingWorkedOn_old[i]) {
                frameList_old[i]->Put(statisticsName_, eventStatisticsForFrame[i]);
            }
        }
    }    
    
    
    log_debug("finished.");
    
    std::size_t framesPushed=0;
    for (std::size_t identifier=0;identifier<frameList_old.size();++identifier)
    {
        if (frameIsBeingWorkedOn_old[identifier]) {
            log_debug("putting photons into frame %zu...", identifier);
            frameList_old[identifier]->Put(photonSeriesMapName_, photonsForFrameList_old[identifier]);
        }
        
        log_debug("pushing frame number %zu...", identifier);
        PushFrame(frameList_old[identifier]);
        ++framesPushed;
    }
    
    return framesPushed;
}

namespace {
    bool ParticleHasMuonDaughter(const I3MCTree::const_iterator &particle_it, const I3MCTree &mcTree)
    {
        BOOST_FOREACH( const I3Particle & daughter, mcTree.children(*particle_it))
        {
            if ((daughter.GetType()==I3Particle::MuMinus) ||
                (daughter.GetType()==I3Particle::MuPlus))
                return true;
        }
        return false;
    }

}

template <typename OutputMapType>
bool I3CLSimModule<OutputMapType>::ShouldDoProcess(I3FramePtr frame)
{
    return true;
}

template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::Process()
{
    I3FramePtr frame = PopFrame();
    if (!frame) return;
    
    if (frame->GetStop() == I3Frame::Geometry)
    {
        // special handling for Geometry frames
        // these will trigger a full re-initialization of OpenCL
        // (currently a second Geometry frame triggers a fatal error
        // in DigestGeometry()..)

        DigestGeometry(frame);
        PushFrame(frame);
        return;
    }
    
    // if the cache is empty and the frame stop is not Physics/DAQ, we can immediately push it
    // (and not add it to the cache)
    if ((frameList_.empty()) && (workOnTheseStops_set_.count(frame->GetStop()) == 0) )
    {
        PushFrame(frame);
        return;
    }
    
    if ((totalEnergyToProcess_ > 0) && (!std::isnan(totalEnergyToProcess_)))
    {
        double totalLightEnergyInFrame = GetLightSourceEnergy(frame);
        if (totalSimulatedEnergy_ + totalSimulatedEnergySecondFlush_ + totalLightEnergyInFrame < totalEnergyToProcess_ / 2.)
        {
            maxNumParallelEvents_++;
            totalSimulatedEnergy_ += totalLightEnergyInFrame;
        }
        else if (totalSimulatedEnergy_ + totalSimulatedEnergySecondFlush_ + totalLightEnergyInFrame < totalEnergyToProcess_)
        {
            maxNumParallelEventsSecondFlush_++;
            totalSimulatedEnergySecondFlush_ += totalLightEnergyInFrame;
        }
        log_debug("Energy in Frame = %f GeV", totalLightEnergyInFrame);
    }
    
    // it's either Physics or something else..
    if (frameListPhysicsFrameCounter_ < maxNumParallelEvents_)
    {
        // we currently treat physics and other frames/empty Physics
        // frames the same
        //const bool isPhysicsFrame =
        DigestOtherFrame(frame);
        frameListPhysicsFrameCounter_++;
    }
    else if (frameListPhysicsFrameCounter_ < maxNumParallelEvents_ + maxNumParallelEventsSecondFlush_) // maxNumParallelEvents_*2) 
    {
        // keep a second buffer so we have it available once 
        // the first buffer has finished processing
        frameList2_.push_back(frame);
        frameListPhysicsFrameCounter_++;
    }


    if (frameListPhysicsFrameCounter_ >= maxNumParallelEvents_ + maxNumParallelEventsSecondFlush_) //maxNumParallelEvents_*2) 
    {
        log_debug("Flushing results for a total energy of %f GeV for %" PRIu64 " particles",
                 totalSimulatedEnergyForFlush_/I3Units::GeV, totalNumParticlesForFlush_);
             
        totalSimulatedEnergyForFlush_= 0.;
        totalNumParticlesForFlush_=0;
    
        // this will finish processing the first buffer and
        // push all its frames
        const std::size_t framesPushed =
            FlushFrameCache();
        frameListPhysicsFrameCounter_ -= framesPushed;

        // now immediately start processing frames from the second buffer
        for (std::size_t i=0;i<frameList2_.size();++i) {
            DigestOtherFrame(frameList2_[i]);
        }
        frameList2_.clear();
        if ((totalEnergyToProcess_ > 0) && (!std::isnan(totalEnergyToProcess_)))
        {
            totalSimulatedEnergy_ = 0.;
            maxNumParallelEvents_ = 1;
            std::swap(totalSimulatedEnergy_, totalSimulatedEnergySecondFlush_);
            std::swap(maxNumParallelEvents_, maxNumParallelEventsSecondFlush_);
        }

            
    
        log_debug("============== CACHE FLUSHED ================");
    }
}

template <typename OutputMapType>
double I3CLSimModule<OutputMapType>::GetLightSourceEnergy(I3FramePtr frame)
{
    I3MCTreeConstPtr MCTree;
    I3CLSimFlasherPulseSeriesConstPtr flasherPulses;
    double totalLightSourceEnergy = 0;
    
    if (MCTreeName_ != "")
        MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);
    if (flasherPulseSeriesName_ != "")
        flasherPulses = frame->Get<I3CLSimFlasherPulseSeriesConstPtr>(flasherPulseSeriesName_);
    if (!MCTree) 
    {   
        log_warn("Frame will not be processed cause MCTree not present.");
        return totalLightSourceEnergy;
    }
    if (flasherPulses) 
        log_fatal("Flashers! Cannot calculate how much energy is deposited in detector. "
                  "Set MaxNumParallelEvents > 0 and totalEnergyToProcess_ = 0");
    
    std::deque<I3CLSimLightSource> lightSources;
    std::deque<double> timeOffsets;
    if (MCTree) ConvertMCTreeToLightSources(*MCTree, lightSources, timeOffsets);
    
    for (std::size_t i=0;i<lightSources.size();++i)
    {
        const I3CLSimLightSource &lightSource = lightSources[i];
        if (lightSource.GetType() == I3CLSimLightSource::Particle)
        {
            
            const I3Particle &particle = lightSource.GetParticle();
            if (particle.GetType() == I3Particle::MuMinus || particle.GetType() == I3Particle::MuPlus)
            {
                // This is the upper estimate of the muon energy loss due to ionization
                // We use the upper estimate to be conservative about how much 
                // energy ends up in the detetor.
                // It is taken from I3MuonSlicer Line 306. It originally comes from PPC.
                totalLightSourceEnergy += (0.21+8.8e-3*log(particle.GetEnergy()/I3Units::GeV)/log(10.))*(I3Units::GeV/I3Units::m) * particle.GetLength();
            }
            else
            {
                totalLightSourceEnergy += particle.GetEnergy();
            }
            // log_info_stream(particle);
        }
        
    }
    
    return totalLightSourceEnergy;
}

template <typename OutputMapType>
bool I3CLSimModule<OutputMapType>::DigestOtherFrame(I3FramePtr frame, bool startThread)
{
    log_trace("%s", __PRETTY_FUNCTION__);
     
    frameList_.push_back(frame);
    photonsForFrameList_.push_back(boost::make_shared<OutputMapType>());
    currentPhotonIdForFrame_.push_back(0);
    std::size_t currentFrameListIndex = frameList_.size()-1;
    maskedOMKeys_.push_back(std::set<ModuleKey>()); // insert an empty ModuleKey mask
    
    // check if we got a geometry before starting to work
    if (!geometryIsConfigured_)
        log_fatal("Received Physics frame before Geometry frame");
    
    if (startThread)
    {
        // start the connector thread if necessary
        if (!threadObj_) {
            log_debug("No thread found running during Physics(), starting one.");
            StartThread();
        }
    }
    
    // a few cases where we don't work with the frame:
    
    //// not our designated Stop
    if (workOnTheseStops_set_.count(frame->GetStop()) == 0) {
        // nothing to do for this frame, it is chached, however
        frameIsBeingWorkedOn_.push_back(false); // do not touch this frame, just push it later on
        return false;
    }
    
    // should we process it? (conditional module)
    const bool shouldDoProcess_fromConditionalModule =
    I3ConditionalModule::ShouldDoProcess(frame);
    
    if (!shouldDoProcess_fromConditionalModule) {
        frameIsBeingWorkedOn_.push_back(false); // do not touch this frame, just push it later on
        return false;
    }
    
    // does it include some work?
    I3MCTreeConstPtr MCTree;
    I3CLSimFlasherPulseSeriesConstPtr flasherPulses;
    
    if (MCTreeName_ != "")
        MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);
    if (flasherPulseSeriesName_ != "")
        flasherPulses = frame->Get<I3CLSimFlasherPulseSeriesConstPtr>(flasherPulseSeriesName_);

    if ((!MCTree) && (!flasherPulses)) {
        // ignore frames without any MCTree and/or Flashers
        frameIsBeingWorkedOn_.push_back(false); // do not touch this frame, just push it later on
        return false;
    }

    I3VectorOMKeyConstPtr omKeyMask;
    I3VectorModuleKeyConstPtr moduleKeyMask;
    if (omKeyMaskName_ != "") {
        omKeyMask = frame->Get<I3VectorOMKeyConstPtr>(omKeyMaskName_);
        
        if (!omKeyMask) {
            moduleKeyMask = frame->Get<I3VectorModuleKeyConstPtr>(omKeyMaskName_);
        }
    }

    
    // work with this frame!
    frameIsBeingWorkedOn_.push_back(true); // this frame will receive results (->Put() will be called later)
    
    std::deque<I3CLSimLightSource> lightSources;
    std::deque<double> timeOffsets;
    if (MCTree) ConvertMCTreeToLightSources(*MCTree, lightSources, timeOffsets);
    if (flasherPulses) ConvertFlasherPulsesToLightSources(*flasherPulses, lightSources, timeOffsets);
    
    // support both vectors of OMKeys and vectors of ModuleKeys
    
    if (omKeyMask) {
        // assign the current OMKey mask if there is one
        BOOST_FOREACH(const OMKey &key, *omKeyMask) {
            maskedOMKeys_.back().insert(ModuleKey(key.GetString(), key.GetOM()));
        }
    }
    
    if (moduleKeyMask) {
        // assign the current ModuleKey mask if there is one
        BOOST_FOREACH(const ModuleKey &key, *moduleKeyMask) {
            maskedOMKeys_.back().insert(key);
        }
    }
   
    for (std::size_t i=0;i<lightSources.size();++i)
    {
        const I3CLSimLightSource &lightSource = lightSources[i];
        const double timeOffset = timeOffsets[i];

        if (lightSource.GetType() == I3CLSimLightSource::Particle)
        {
            const I3Particle &particle = lightSource.GetParticle();
            
            totalSimulatedEnergyForFlush_ += particle.GetEnergy();
            totalNumParticlesForFlush_++;
        }
        
        geant4ParticleToStepsConverter_->EnqueueLightSource(lightSource, currentParticleCacheIndex_);

        if (particleCache_.find(currentParticleCacheIndex_) != particleCache_.end())
            log_fatal("Internal error. Particle cache index already used.");
        
        particleCacheEntry &cacheEntry = 
        particleCache_.insert(std::make_pair(currentParticleCacheIndex_, particleCacheEntry())).first->second;
        
        cacheEntry.frameListEntry = currentFrameListIndex;
        cacheEntry.timeShift = timeOffset;
        if (lightSource.GetType() == I3CLSimLightSource::Particle) {
            cacheEntry.particleMajorID = lightSource.GetParticle().GetMajorID();
            cacheEntry.particleMinorID = lightSource.GetParticle().GetMinorID();
        } else {
            cacheEntry.particleMajorID = 0; // flashers, etc. do get ID 0,0
            cacheEntry.particleMinorID = 0;
        }
        
        // make a new index. This will eventually overflow,
        // but at that time, index 0 should be unused again.
        ++currentParticleCacheIndex_;
        if (currentParticleCacheIndex_==0) ++currentParticleCacheIndex_; // never use index==0
    }
    
    lightSources.clear();

    return true;
}

template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::Finish()
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    log_debug("Flushing results for a total energy of %fGeV for %" PRIu64 " particles",
             totalSimulatedEnergyForFlush_/I3Units::GeV, totalNumParticlesForFlush_);

    totalSimulatedEnergyForFlush_=0.;
    totalNumParticlesForFlush_=0;

    std::size_t framesPushed = 0;
    while (frameListPhysicsFrameCounter_ > 0) {
        framesPushed = FlushFrameCache();
        frameListPhysicsFrameCounter_ -= framesPushed;

        log_info("Flushing I3Tray..");
        Flush();

        // finish frames in 2nd buffer if there are any
        if (frameList2_.size()>0) {
            for (std::size_t i=0;i<frameList2_.size();++i) {
                DigestOtherFrame(frameList2_[i]);
            }
            frameList2_.clear();

            framesPushed = FlushFrameCache();
            frameListPhysicsFrameCounter_ -= framesPushed;

            log_info("Flushing I3Tray (again)..");
            Flush();
        }
    }

    log_info("I3CLSimModule is done.");

    // add some summary information to a potential I3SummaryService
    I3MapStringDoublePtr summary = context_.Get<I3MapStringDoublePtr>("I3SummaryService");
    if (summary) {
        const std::string prefix = "I3CLSimModule_" + GetName() + "_";
        
        for (std::size_t i=0; i<openCLStepsToPhotonsConverters_.size(); ++i)
        {
            const std::string postfix = (openCLStepsToPhotonsConverters_.size()==1)?"":"_"+boost::lexical_cast<std::string>(i);
            
            const double totalNumPhotonsGenerated = openCLStepsToPhotonsConverters_[i]->GetTotalNumPhotonsGenerated();
            const double totalDeviceTime = static_cast<double>(openCLStepsToPhotonsConverters_[i]->GetTotalDeviceTime())*I3Units::ns;
            const double totalHostTime = static_cast<double>(openCLStepsToPhotonsConverters_[i]->GetTotalHostTime())*I3Units::ns;
            
            (*summary)[prefix+"TotalDeviceTime"           +postfix] = totalDeviceTime;
            (*summary)[prefix+"TotalHostTime"             +postfix] = totalHostTime;
            (*summary)[prefix+"NumKernelCalls"            +postfix] = openCLStepsToPhotonsConverters_[i]->GetNumKernelCalls();
            (*summary)[prefix+"TotalNumPhotonsGenerated"  +postfix] = totalNumPhotonsGenerated;
            (*summary)[prefix+"TotalNumPhotonsAtDOMs"     +postfix] = openCLStepsToPhotonsConverters_[i]->GetTotalNumPhotonsAtDOMs();
            
            (*summary)[prefix+"AverageDeviceTimePerPhoton"+postfix] = totalDeviceTime/totalNumPhotonsGenerated;
            (*summary)[prefix+"AverageHostTimePerPhoton"  +postfix] = totalHostTime/totalNumPhotonsGenerated;
            (*summary)[prefix+"DeviceUtilization"         +postfix] = totalDeviceTime/totalHostTime;
        }
        
    }

}




//////////////

template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::ConvertMCTreeToLightSources(const I3MCTree &mcTree,
                                                std::deque<I3CLSimLightSource> &lightSources,
                                                std::deque<double> &timeOffsets)
{
    for (I3MCTree::const_iterator particle_it = mcTree.begin();
         particle_it != mcTree.end(); ++particle_it)
    {
        const I3Particle &particle_ref = *particle_it;

        // In-ice particles only
        if (particle_ref.GetLocationType() != I3Particle::InIce) continue;
        
        // ignore particles with shape "Dark"
        if (particle_ref.GetShape() == I3Particle::Dark) continue;

        // skip primaries that are clearly outside the ice 
        // (those are probably cosmic rays that get marked as "InIce"
        // by ucr-icetray)
        if (particle_ref.GetShape() == I3Particle::Primary) {
            if (particle_ref.GetZ() >= mediumProperties_->GetAirZCoord()) continue;
        }

        // check particle type
        const bool isMuon = (particle_ref.GetType() == I3Particle::MuMinus) || (particle_ref.GetType() == I3Particle::MuPlus);
        const bool isNeutrino = particle_ref.IsNeutrino();
        const bool isTrack = particle_ref.IsTrack();
        

        // mmc-icetray currently stores continuous loss entries as "unknown"
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        // The ContinuousEnergyLoss type is only defined in the I3Particle version
        // that also introduced "I3PARTICLE_SUPPORTS_PDG_ENCODINGS". In order
        // to make clsim work on previous versions, disable the check for ContinuousEnergyLoss.
        // (the way mmc-icetray is implemented, it would show up as "unknown" anyway.)

        const bool isContinuousLoss = (particle_ref.GetType() == I3Particle::unknown) ||
                                      (particle_ref.GetType() == I3Particle::ContinuousEnergyLoss);
#else
        const bool isContinuousLoss = (particle_ref.GetType() == I3Particle::unknown);
#endif
        
        // ignore continuous loss entries
        if (isContinuousLoss) {
            log_debug("ignored a continuous loss I3MCTree entry");
            continue;
        }
        
        // always ignore neutrinos
        if (isNeutrino) continue;
        
        // ignore muons if requested
        if ((ignoreMuons_) && (isMuon)) continue;
        
        if (!isTrack && detectorHull_)
        {
            auto intersection = detectorHull_->GetIntersection(particle_ref.GetPos(), particle_ref.GetDir());
            
            if (!(intersection.first <= 0) || !(intersection.second >= 0))
            {
                log_debug("Ignored a non-track that is >%fm away from the detector hull",
                          closestDOMDistanceCutoff_);
                continue;
            }
        }
        
        // make a copy of the particle, we may need to change its length
        I3Particle particle = particle_ref;
        
        if (isTrack && detectorHull_)
        {
            bool nostart = false;
            bool nostop = false;
            double particleLength = particle.GetLength();
            
            if (std::isnan(particleLength)) {
                // assume infinite track (starting at given position)
                nostop = true;
            } else if (particleLength < 0.) {
                log_warn("got track with negative length. assuming it starts at given position.");
                nostop = true;
            } else if (particleLength == 0.){
                // zero length: starting track
                nostop = true;
            }
            
            auto intersection = detectorHull_->GetIntersection(particle_ref.GetPos(), particle_ref.GetDir());
            // Skip the track if it does not intersect the hull at all,
            // stops before it reaches the hull, or starts after the hull
            if ((std::isnan(intersection.first) && std::isnan(intersection.second))
                || (!nostop && (intersection.first > particleLength))
                || (intersection.second < 0))
            {
                log_debug("Ignored a track that is always at least %fm away from the closest DOM.",
                          closestDOMDistanceCutoff_);
                continue;
            }
        }
        
        // ignore muons with muons as child particles
        // -> those already ran through MMC(-recc) or
        // were sliced with I3MuonSlicer. Only add their
        // children.
        if ((!ignoreMuons_) && (isMuon)) {
            if (ParticleHasMuonDaughter(particle_it, mcTree)) {
                log_warn("particle has muon as daughter(s) but is not \"Dark\". Strange. Ignoring.");
                continue;
            }
        }
        
        // simulate the particle around time 0, add the offset later
        const double particleTime = particle.GetTime();
        particle.SetTime(0.);
        
        lightSources.push_back(I3CLSimLightSource(particle));
        timeOffsets.push_back(particleTime);
    }
    
    
}


template <typename OutputMapType>
void I3CLSimModule<OutputMapType>::ConvertFlasherPulsesToLightSources(const I3CLSimFlasherPulseSeries &flasherPulses,
                                                       std::deque<I3CLSimLightSource> &lightSources,
                                                       std::deque<double> &timeOffsets)
{
    BOOST_FOREACH(const I3CLSimFlasherPulse &flasherPulse, flasherPulses)
    {
        lightSources.push_back(I3CLSimLightSource(flasherPulse));
        timeOffsets.push_back(0.);
    }
}
