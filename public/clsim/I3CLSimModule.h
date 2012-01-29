/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3CLSimModule.h
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

#ifndef I3CLSIMMODULE_H_INCLUDED
#define I3CLSIMMODULE_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "phys-services/I3RandomService.h"

#include "clsim/I3CLSimOpenCLDevice.h"

#include "clsim/I3CLSimRandomValue.h"
#include "clsim/I3CLSimWlenDependentValue.h"

#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"
#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

#include "clsim/I3CLSimLightSourceParameterization.h"

#include "clsim/I3Photon.h"
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
    ~I3CLSimModule();
    
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
    
    
private:
    /**
     * The module needs to process Physics frames
     */
    void DigestOtherFrame(I3FramePtr frame);
    
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
    
    /// Parameter: An instance of I3CLSimWlenDependentValue describing the reciprocal weight a photon gets assigned as a function of its wavelength.
    ///            You can set this to the wavelength depended acceptance of your DOM to pre-scale the number of generated photons.
    I3CLSimWlenDependentValueConstPtr wavelengthGenerationBias_;

    /// Parameter: An instance of I3CLSimMediumProperties describing the ice/water properties.
    I3CLSimMediumPropertiesPtr mediumProperties_;

    /// Parameter: Maximum number of events that will be processed by the GPU in parallel.
    unsigned int maxNumParallelEvents_;

    /// Parameter: A vector of I3CLSimOpenCLDevice objects, describing the devices to be used for simulation.
    I3CLSimOpenCLDeviceSeries openCLDeviceList_;

    /// Parameter: The DOM radius used during photon tracking.
    double DOMRadius_;

    /// Parameter: Ignore string numbers < 1 and OM numbers > 60. (AMANDA and IceTop)
    bool ignoreNonIceCubeOMNumbers_;

    /// Parameter: Name of the I3MCTree frame object. All particles except neutrinos will be read from this tree.
    std::string MCTreeName_;
    
    /// Parameter: Name of the I3CLSimPhotonSeriesMap frame object that will be written to the frame.
    std::string photonSeriesMapName_;
    
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
    void FlushFrameCache();
    void AddPhotonsToFrames(const I3CLSimPhotonSeries &photons);

    // statistics will be collected here:
    std::map<uint32_t, uint64_t> photonNumGeneratedPerParticle_;
    std::map<uint32_t, double> photonWeightSumGeneratedPerParticle_;

    std::map<uint32_t, uint64_t> photonNumAtOMPerParticle_;
    std::map<uint32_t, double> photonWeightSumAtOMPerParticle_;


    
    
    bool geometryIsConfigured_;
    uint32_t currentParticleCacheIndex_;
    double totalSimulatedEnergyForFlush_;
    uint64_t totalNumParticlesForFlush_;
    
    // this is calculated from wavelengthGenerationBias:
    I3CLSimRandomValueConstPtr wavelengthGenerator_;

    I3CLSimSimpleGeometryFromI3GeometryPtr geometry_;
    std::vector<I3CLSimStepToPhotonConverterOpenCLPtr> openCLStepsToPhotonsConverters_;
    I3CLSimLightSourceToStepConverterGeant4Ptr geant4ParticleToStepsConverter_;
    
    // list of all currently held frames, in order
    unsigned int frameListPhysicsFrameCounter_;
    std::vector<I3FramePtr> frameList_;
    std::vector<I3PhotonSeriesMapPtr> photonsForFrameList_;
    std::vector<int32_t> currentPhotonIdForFrame_;
    
    struct particleCacheEntry
    {
        std::size_t frameListEntry; // pointer to the frame list by entry number
        uint64_t particleMajorID;
        int particleMinorID;
    };
    
    // list of all particles (with pointrs to their frames)
    // currently being simulated
    std::map<uint32_t, particleCacheEntry> particleCache_;
    
    SET_LOGGER("I3CLSimModule");
};

#endif //I3CLSIMMODULE_H_INCLUDED
