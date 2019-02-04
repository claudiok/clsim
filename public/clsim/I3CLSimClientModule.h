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
 * @file I3CLSimClientModule.h   
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMCLIENTMODULE_H_INCLUDED
#define I3CLSIMCLIENTMODULE_H_INCLUDED

#include <atomic>
#include <deque>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/thread.hpp>

#include "icetray/I3ConditionalModule.h"
#include "simclasses/I3CompressedPhoton.h"
#include "sim-services/MCPEMCPulseTools.hpp"
#include "sim-services/I3CosmicEventGenerator.h"
#include "clsim/I3CLSimFlasherPulse.h"

I3_FORWARD_DECLARATION(I3CLSimClient);
I3_FORWARD_DECLARATION(I3CLSimLightSource);
I3_FORWARD_DECLARATION(I3CLSimLightSourceToStepConverterAsync);
I3_FORWARD_DECLARATION(I3CLSimPhotonToMCPEConverter);

/**
 * @brief
 */
class I3CLSimClientModule : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3CLSimClientModule(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    virtual ~I3CLSimClientModule();
    
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
    bool DigestOtherFrame(I3FramePtr frame);
    
    /**
     * The module needs to process Geometry frames
     */
    void DigestGeometry(I3FramePtr frame);

    // parameters
    
    /// Parameter: work on MCTrees found in the stream types ("stops") specified in this list
    std::vector<I3Frame::Stream> workOnTheseStops_;
    std::set<I3Frame::Stream> workOnTheseStops_set_; // will be converted to a set internally

    /// Parameter: Name of the I3MCTree frame object. All particles except neutrinos will be read from this tree.
    std::string MCTreeName_;

    /// Parameter: Name of the I3CLSimFlasherPulseSeries frame object. Flasher pulses will be read from this object.
    /// Set this to the empty string to disable flashers.
    std::string flasherPulseSeriesName_;

    /// Parameter: Name of the I3CompressedPhotonSeriesMap frame object that will be written to the frame.
    std::string photonSeriesMapName_;
    
    /// Parameter: Name of the I3MCPESeriesMap frame object that will be written to the frame.
    std::string mcpeSeriesMapName_;

    /// Parameter: Name of a I3VectorOMKey with masked OMKeys. DOMs in this list will not record I3Photons.
    std::string omKeyMaskName_;
    
    /// Parameter: If set to True, muons will not be propagated.
    bool ignoreMuons_;

    /// Parameter: Collect statistics in this frame object (e.g. number of photons generated or reaching the DOMs)
    std::string statisticsName_;
    bool collectStatistics_;

    /// Parameter: Ignore all OMKeys with these subdetector names
    std::vector<std::string> ignoreSubdetectors_;

    /// Parmeter: Limits the OpenCL workgroup size (the number of bunches to be processed in parallel).
    ///   If set to zero (the default) the largest possible workgroup size will be chosen.
    uint32_t limitWorkgroupSize_;

    /// Parameter: do not even start light from sources that do not have any DOMs closer to
    ///   to them than this distance. (default is 300m)
    double closestDOMDistanceCutoff_;
    boost::shared_ptr<I3Surfaces::ExtrudedPolygon> detectorHull_;

    boost::shared_ptr<I3CosmicEventGenerator> cosmicGenerator_;    

    I3CLSimLightSourceToStepConverterAsyncPtr stepGenerator_;
    
    I3CLSimPhotonToMCPEConverterPtr mcpeGenerator_;

private:
    // default, assignment, and copy constructor declared private
    I3CLSimClientModule();
    I3CLSimClientModule(const I3CLSimClientModule&);
    I3CLSimClientModule& operator=(const I3CLSimClientModule&);
    
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
    
    // helper functions
    std::size_t FlushFrameCache();
    void ConvertMCTreeToLightSources(const I3MCTree &mcTree,
                                     std::deque<I3CLSimLightSource> &lightSources,
                                     std::deque<double> &timeOffsets);
    void ConvertFlasherPulsesToLightSources(const I3CLSimFlasherPulseSeries &flasherPulses,
                                            std::deque<I3CLSimLightSource> &lightSources,
                                            std::deque<double> &timeOffsets);

    void CosmicGeneraterToLightSources(I3MCTree &mcTree,
                                       std::deque<I3CLSimLightSource> &lightSources,
                                       std::deque<double> &timeOffsets);
  
    bool FrameContainsWork(const I3FramePtr&, I3MCTreeConstPtr&, I3CLSimFlasherPulseSeriesConstPtr&);


    boost::shared_ptr<I3CLSimClient> stepsToPhotonsConverter_;

public:
    struct frameCacheEntry 
    {
        uint32_t frameId;
        uint32_t currentPhotonId;
        uint32_t numPendingParticles;
        I3FramePtr frame;
        std::set<ModuleKey> ignoreModules;
        std::vector<uint32_t> particles;
        boost::shared_ptr<I3CompressedPhotonSeriesMap> photons;
        typedef std::map<OMKey, MCHitMerging::MCPEStream> MCPEStreamMap;
        boost::shared_ptr<MCPEStreamMap> hits;
        
        frameCacheEntry(uint32_t identifier, const I3FramePtr &f)
            : frameId(identifier), currentPhotonId(0), numPendingParticles(0), frame(f)
        {}
        
        uint32_t key() const { return frameId; }
    };
    
private:
    boost::mutex frameCache_mutex_;
    std::atomic<bool> newFramesAvailable_;
    uint32_t currentFrameId_;
    uint32_t currentParticleCacheIndex_;
    uint32_t currentBunchId_;
    std::deque<frameCacheEntry> frameCache_;
    decltype(frameCache_)::difference_type framesInKernel_;
    
    // left: frame id, right: bunch id
    typedef boost::bimaps::bimap< boost::bimaps::multiset_of<uint32_t>, boost::bimaps::multiset_of<uint32_t> > bimap;
    bimap framesForBunches_;
public:
    struct particleCacheEntry
    {
        uint32_t particleId;
        uint32_t frameId; // pointer to the frame list by entry number
        struct photonStatistics { 
            uint64_t count;
            double weightSum;
            photonStatistics() : count(0), weightSum(0) {}
        };
        photonStatistics generatedPhotons, photonsAtDOMs;
        double timeShift; // optional time that needs to be added to the final output photon
        uint64_t particleMajorID;
        int particleMinorID;
        
        particleCacheEntry(uint32_t id, const I3CLSimLightSource &lightSource, uint32_t frame, double dt);
        
        uint32_t key() const { return particleId; }
    };
private:
    // list of all particles (with pointrs to their frames)
    // currently being simulated
    std::deque<particleCacheEntry> particleCache_;

    SET_LOGGER("I3CLSimClientModule");
};

#endif //I3CLSIMCLIENTMODULE_H_INCLUDED
