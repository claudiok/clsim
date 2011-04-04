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

#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimPhotonSeriesToPhotonSeriesMapConverter.h"
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"
#include "clsim/I3CLSimParticleToStepConverterGeant4.h"

#include "clsim/I3CLSimParticleParameterization.h"

#include "clsim/I3Photon.h"

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include <vector>
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
     * The module needs to process Physics frames
     */
    virtual void Physics(I3FramePtr frame);

    /**
     * The module needs to process Geometry frames
     */
    virtual void Geometry(I3FramePtr frame);

    /**
     * Warn the user if the module is aborted prematurely
     */
    virtual void Finish();
    
    
private:
    // parameters
    
    /// Parameter: An instance I3CLSimParticleParameterizationSeries specifying the fast simulation parameterizations to be used.
    I3CLSimParticleParameterizationSeries parameterizationList_;
    
    /// Parameter: A random number generating service (derived from I3RandomService).
    I3RandomServicePtr randomService_;

    /// Parameter: An instance of I3CLSimMediumProperties describing the ice/water properties.
    I3CLSimMediumPropertiesPtr mediumProperties_;

    /// Parameter: Maximum number of events that will be processed by the GPU in parallel.
    unsigned int maxNumParallelEvents_;
    
    /// Parameter: Name of the I3MCTree frame object. All particles except neutrinos will be read from this tree.
    std::string MCTreeName_;
    
    /// Parameter: Name of the I3CLSimPhotonSeriesMap frame object that will be written to the frame.
    std::string photonSeriesMapName_;
    
    /// Parameter: If set to True, muons will not be propagated.
    bool ignoreMuons_;
    
    
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
    uint64_t numBunchesSentToOpenCL_;

    
    // helper functions
    void FlushFrameCache();
    void AddPhotonsToFrames(const I3CLSimPhotonSeries &photons);

    
    bool geometryIsConfigured_;
    uint32_t currentParticleCacheIndex_;
    double totalSimulatedEnergyForFlush_;
    uint64_t totalNumParticlesForFlush_;
    
    I3CLSimPhotonSeriesToPhotonSeriesMapConverterPtr photonSeriesMapConverter_;
    I3CLSimSimpleGeometryFromI3GeometryPtr geometry_;
    I3CLSimStepToPhotonConverterOpenCLPtr openCLStepsToPhotonsConverter_;
    I3CLSimParticleToStepConverterGeant4Ptr geant4ParticleToStepsConverter_;
    
    // list of all currently held frames, in order
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
