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
 * @file I3CLSimModule.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMMODULE_H_INCLUDED
#define I3CLSIMMODULE_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "dataclasses/I3Vector.h"
#include "dataclasses/physics/I3MCTree.h"

#include "icetray/OMKey.h"

#include "dataclasses/ModuleKey.h"

#include "phys-services/I3RandomService.h"

#include "clsim/I3CLSimOpenCLDevice.h"

#include "clsim/random_value/I3CLSimRandomValue.h"
#include "clsim/function/I3CLSimFunction.h"

#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimSpectrumTable.h"

#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"
#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

#include "clsim/I3CLSimLightSourceParameterization.h"

#include "clsim/I3Photon.h"

#include "clsim/I3CLSimPhotonHistory.h"
#include "clsim/I3CLSimEventStatistics.h"

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include <vector>
#include <set>
#include <map>
#include <string>



/**
 * @brief
 */
class I3CLSimModule : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3CLSimModule(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    virtual ~I3CLSimModule();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
    /**
     * The module needs full control of all streams
     * -> implement Process()!
     */
    virtual void Process();

    /**
     * Warn the user if the module is aborted prematurely
     */
    virtual void Finish();
    

    /**
     * Hack to allow buffering. This ShouldDoProcess overrides
     * the I3ConditionalModule implementation and always returns
     * true. The original ShouldDoProcess() is used in our
     * Process() in order mark frames in case they should not be
     * touched. They will be put in the buffer in any case, however.
     *
     * This should retain frame ordering.
     */
    virtual bool ShouldDoProcess(I3FramePtr frame);
    
private:
    /**
     * The module needs to process Physics frames
     */
    bool DigestOtherFrame(I3FramePtr frame, bool startThread=true);
    
    /**
     * The module needs to process Geometry frames
     */
    void DigestGeometry(I3FramePtr frame);

    // parameters
    
    /// Parameter: work on MCTrees found in the stream types ("stops") specified in this list
    std::vector<I3Frame::Stream> workOnTheseStops_;
    std::set<I3Frame::Stream> workOnTheseStops_set_; // will be converted to a set internally
    
    /// Parameter: An instance I3CLSimLightSourceParameterizationSeries specifying the fast simulation parameterizations to be used.
    I3CLSimLightSourceParameterizationSeries parameterizationList_;
    
    /// Parameter: A random number generating service (derived from I3RandomService).
    I3RandomServicePtr randomService_;

    /// Parameter: The wavelength of the generated Cherenkov photons will be generated
    ///            according to a spectrum without dispersion. This does not change
    ///            the total number of photons, only the distribution of wavelengths.
    bool generateCherenkovPhotonsWithoutDispersion_;
    
    /// Parameter: An instance of I3CLSimFunction describing the reciprocal weight a photon gets assigned as a function of its wavelength.
    ///            You can set this to the wavelength depended acceptance of your DOM to pre-scale the number of generated photons.
    I3CLSimFunctionConstPtr wavelengthGenerationBias_;

    /// Parameter: An instance of I3CLSimMediumProperties describing the ice/water properties.
    I3CLSimMediumPropertiesConstPtr mediumProperties_;

    /// Parameter: All spectra that could be requested by an I3CLSimStep.
    /// If set to NULL/None, only spectrum #0 (Cherenkov photons) will be available.
    I3CLSimSpectrumTableConstPtr spectrumTable_;

    /// Parameter: Maximum number of events that will be processed by the GPU in parallel.
    unsigned int maxNumParallelEvents_;

    /// Parameter: A vector of I3CLSimOpenCLDevice objects, describing the devices to be used for simulation.
    I3CLSimOpenCLDeviceSeries openCLDeviceList_;

    /// Parameter: The DOM radius used during photon tracking.
    double DOMRadius_;

    /// Parameter: Specifiy the "oversize factor" (i.e. DOM radius scaling factor).
    double DOMOversizeFactor_;

    /// Parameter: Ignore string numbers < 1 and OM numbers > 60. (AMANDA and IceTop)
    bool ignoreNonIceCubeOMNumbers_;

    /// Parameter: Name of the I3MCTree frame object. All particles except neutrinos will be read from this tree.
    std::string MCTreeName_;

    /// Parameter: Name of the I3CLSimFlasherPulseSeries frame object. Flasher pulses will be read from this object.
    /// Set this to the empty string to disable flashers.
    std::string flasherPulseSeriesName_;

    /// Parameter: Name of the I3CLSimPhotonSeriesMap frame object that will be written to the frame.
    std::string photonSeriesMapName_;

    /// Parameter: Name of a I3VectorOMKey with masked OMKeys. DOMs in this list will not record I3Photons.
    std::string omKeyMaskName_;
    
    /// Parameter: If set to True, muons will not be propagated.
    bool ignoreMuons_;
    
    /// Parameter: Geant4 physics list name. Examples are "QGSP_BERT_EMV" and "QGSP_BERT".
    std::string geant4PhysicsListName_;

    /// Parameter: Maximum change of beta=v/c per Geant4 step.
    double geant4MaxBetaChangePerStep_;

    /// Parameter: Approximate maximum number of Cherenkov photons generated per step by Geant4.
    uint32_t geant4MaxNumPhotonsPerStep_;

    /// Parameter: Collect statistics in this frame object (e.g. number of photons generated or reaching the DOMs)
    std::string statisticsName_;
    bool collectStatistics_;
    
    /// Parameter: Ignore all OMKeys with these string IDs
    std::vector<int> ignoreStrings_;

    /// Parameter: Ignore all OMKeys with these DOM IDs
    std::vector<unsigned int> ignoreDomIDs_;

    /// Parameter: Ignore all OMKeys with these subdetector names
    std::vector<std::string> ignoreSubdetectors_;

    /// Parmeter: If you have a geometry with multiple OMs per floor (e.g. Antares or KM3NeT-tower-like),
    ///   it will internally get split into subdetectors to save memory on the GPU. This is 100% transparent.
    ///   By default the split is according to the "floor index" of an OM on a floor. If you enable this
    ///   option, the split will also be done according to the x-y projected positions of the OMs per string.
    ///   This may be necessary for "tower" geometries.
    bool splitGeometryIntoPartsAcordingToPosition_;

    /// Parameter: Split off DeepCore as its own two subdetectors (upper and lower part).
    ///  This may save constant memory on your GPU.
    ///  Assumes that strings [79..86] are DeepCore strings with an upper part at z>-30m
    ///  and a lower part at z<-30m.
    bool useHardcodedDeepCoreSubdetector_;

    /// Parmeter: Disables or enables double-buffered GPU usage. Double buffering will use
    ///   two command queues and two sets of input and output buffers in order to transfer
    ///   data to the GPU while a kernel is executing on the other buffer.
    ///   This has been observed to yield empty results results on older drivers for the nVidia
    ///   architecture, so it is disabled by default.
    ///   
    ///   Before enabling this for a certain driver/hardware combination, make sure that both correct results
    ///   are returned. Most of the time the second buffer results are always empty, so this error should be
    ///   easy to observe.
    bool enableDoubleBuffering_;
    
    /// Parmeter: Enables double-precision support in the kernel. This slows down calculations and
    ///   requires more memory. The performance hit is minimal on CPUs but up to an order
    ///   of magnitude on GPUs.
    bool doublePrecision_;
    
    /// Parmeter: Configures behaviour for photons that hit a DOM. If this is true (the default)
    ///   photons will be stopped once they hit a DOM. If this is false, they continue to
    ///   propagate.
    bool stopDetectedPhotons_;

    /// Parmeter: Saves all photons, even if they don't hit a DOM. Cannot be used with "StopDetectedPhotons".
    bool saveAllPhotons_;

    /// Parmeter: Sets the prescale factor of photons being generated in "saveAllPhotons" mode.
    ///   Only this fraction of photons is actually generated.
    double saveAllPhotonsPrescale_;

    /// Parameter: Sets the number of absorption lengths each photon
    ///   should be propagated. If set to NaN (the default),
    ///   the number is sampled from an exponential distribution.
    ///   (This is what you want for "normal" propagation.)
    ///
    ///   Use this override for table-making.
    double fixedNumberOfAbsorptionLengths_;

    /// Parameter:  Sets the "pancake" factor for DOMs. For standard
    ///  oversized-DOM simulations, this should be the
    ///  radius oversizing factor. This will flatten the
    ///  DOM in the direction parallel to the photon.
    ///  The DOM will have a pancake-like shape, elongated
    ///  in the directions perpendicular to the photon direction.
    ///
    ///  The DOM radius (supplied by the geometry) must also include
    ///  the oversizing factor.
    double pancakeFactor_;
    
    /// Parmeter: Sets the number of scattering step positions that are saved for a photon hitting
    ///   a DOM. The last N photons are saved if there are more scattering points than available entries.
    uint32_t photonHistoryEntries_;

    /// Parmeter: Limits the OpenCL workgroup size (the number of bunches to be processed in parallel).
    ///   If set to zero (the default) the largest possible workgroup size will be chosen.
    uint32_t limitWorkgroupSize_;


private:
    // default, assignment, and copy constructor declared private
    I3CLSimModule();
    I3CLSimModule(const I3CLSimModule&);
    I3CLSimModule& operator=(const I3CLSimModule&);
    
    // thread stuff
    void StopThread();
    void StartThread();
    void Thread_starter();
    bool Thread(boost::this_thread::disable_interruption &di);
    boost::shared_ptr<boost::thread> threadObj_;
    boost::condition_variable_any threadStarted_cond_;
    boost::mutex threadStarted_mutex_;
    bool threadStarted_;
    bool threadFinishedOK_;
    std::vector<uint64_t> numBunchesSentToOpenCL_;

    
    // helper functions
    std::size_t FlushFrameCache();
    void ConvertMCTreeToLightSources(const I3MCTree &mcTree,
                                     std::deque<I3CLSimLightSource> &lightSources,
                                     std::deque<double> &timeOffsets);
    void ConvertFlasherPulsesToLightSources(const I3CLSimFlasherPulseSeries &flasherPulses,
                                            std::deque<I3CLSimLightSource> &lightSources,
                                            std::deque<double> &timeOffsets);

    
    // statistics will be collected here:
    std::map<uint32_t, uint64_t> photonNumGeneratedPerParticle_;
    std::map<uint32_t, double> photonWeightSumGeneratedPerParticle_;



    
    
    bool geometryIsConfigured_;
    uint32_t currentParticleCacheIndex_;
    double totalSimulatedEnergyForFlush_;
    uint64_t totalNumParticlesForFlush_;
    
    // this is calculated using wavelengthGenerationBias:
    std::vector<I3CLSimRandomValueConstPtr> wavelengthGenerators_;

    I3CLSimSimpleGeometryFromI3GeometryPtr geometry_;
    std::vector<I3CLSimStepToPhotonConverterOpenCLPtr> openCLStepsToPhotonsConverters_;
    I3CLSimLightSourceToStepConverterGeant4Ptr geant4ParticleToStepsConverter_;
    
    // list of all currently held frames, in order
    std::size_t frameListPhysicsFrameCounter_;
    std::vector<I3FramePtr> frameList_;
    std::deque<I3FramePtr> frameList2_;
    std::vector<I3PhotonSeriesMapPtr> photonsForFrameList_;
    std::vector<int32_t> currentPhotonIdForFrame_;
    std::vector<bool> frameIsBeingWorkedOn_;
    std::vector<std::set<ModuleKey> > maskedOMKeys_;
    
    struct particleCacheEntry
    {
        std::size_t frameListEntry; // pointer to the frame list by entry number
        uint64_t particleMajorID;
        int particleMinorID;
        double timeShift; // optional time that needs to be added to the final output photon
    };
    
    // list of all particles (with pointrs to their frames)
    // currently being simulated
    std::map<uint32_t, particleCacheEntry> particleCache_;
    
    static void AddPhotonsToFrames(const I3CLSimPhotonSeries &photons,
                                   I3CLSimPhotonHistorySeriesConstPtr photonHistories,
                                   const std::vector<I3PhotonSeriesMapPtr> &photonsForFrameList_,
                                   std::vector<int32_t> &currentPhotonIdForFrame_,
                                   const std::vector<I3FramePtr> &frameList_,
                                   const std::map<uint32_t, particleCacheEntry> &particleCache_,
                                   const std::vector<std::set<ModuleKey> > &maskedOMKeys_,
                                   bool collectStatistics_,
                                   std::map<uint32_t, uint64_t> &photonNumAtOMPerParticle,
                                   std::map<uint32_t, double> &photonWeightSumAtOMPerParticle
                                   );

    SET_LOGGER("I3CLSimModule");
};

#endif //I3CLSIMMODULE_H_INCLUDED
