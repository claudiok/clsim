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
 * @file I3CLSimClientModule
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include "dataclasses/geometry/I3ModuleGeo.h"
#include "phys-services/surfaces/ExtrudedPolygon.h"
#include "clsim/I3CLSimClientModule.h"
#include "clsim/I3CLSimServer.h"
#include "clsim/I3CLSimLightSourceToStepConverterAsync.h"
#include "clsim/I3CLSimEventStatistics.h"
#include "clsim/dom/I3CLSimPhotonToMCPEConverter.h"

I3_MODULE(I3CLSimClientModule);

I3CLSimClientModule::particleCacheEntry::particleCacheEntry(uint32_t id, const I3CLSimLightSource &lightSource, uint32_t frame, double dt)
    : particleId(id), frameId(frame), timeShift(dt)
{
    if (lightSource.GetType() == I3CLSimLightSource::Particle) {
        particleMajorID = lightSource.GetParticle().GetMajorID();
        particleMinorID = lightSource.GetParticle().GetMinorID();
    } else {
        particleMajorID = 0; // flashers, etc. do get ID 0,0
        particleMinorID = 0;
    }
}

I3CLSimClientModule::I3CLSimClientModule(const I3Context& context) 
: I3ConditionalModule(context)
{
    // define parameters
    workOnTheseStops_.clear();
    workOnTheseStops_.push_back(I3Frame::DAQ);
    AddParameter("WorkOnTheseStops",
                 "Work on MCTrees found in the stream types (\"stops\") specified in this list",
                 workOnTheseStops_);
    
    AddParameter("ServerAddress",
                 "Address of an I3CLSimServer instance",
                 std::string());
    
    AddParameter("StepGenerator",
                 "Instance of I3CLSimLightSourceToStepConverterAsync",
                 NULL);
    
    AddParameter("CosmicEventGenerator",
                 "Instance of I3CosmicEventGenerator to ",
                 NULL);
    
    MCTreeName_="I3MCTree";
    AddParameter("MCTreeName",
                 "Name of the I3MCTree frame object. All particles except neutrinos will be read from this tree.",
                 MCTreeName_);

    flasherPulseSeriesName_="";
    AddParameter("FlasherPulseSeriesName",
                 "Name of the I3CLSimFlasherPulseSeries frame object. Flasher pulses will be read from this object.\n"
                 "Set this to the empty string to disable flashers.",
                 flasherPulseSeriesName_);

    mcpeSeriesMapName_="";
    AddParameter("MCPESeriesMapName",
                 "Name of the I3MCPESeriesMap frame object that will be written to the frame.",
                 mcpeSeriesMapName_);

    AddParameter("MCPEGenerator",
                 "Instance of I3CLSimPhotonToMCPEConverter",
                 mcpeGenerator_);

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

    AddParameter("IgnoreSubdetectors",
                 "Ignore all OMKeys with these subdetector names",
                 ignoreSubdetectors_);

    statisticsName_="";
    AddParameter("StatisticsName",
                 "Collect statistics in this frame object (e.g. number of photons generated or reaching the DOMs)",
                 statisticsName_);

    closestDOMDistanceCutoff_=300.*I3Units::m;
    AddParameter("ClosestDOMDistanceCutoff",
                 "Do not even start light from sources that do not have any DOMs closer to\n"
                 "to them than this distance.",
                 closestDOMDistanceCutoff_);

    // add an outbox
    AddOutBox("OutBox");

    currentParticleCacheIndex_ = 1;
    currentFrameId_=0;
    currentBunchId_=0;
    framesInKernel_=0;
    newFramesAvailable_=false;
}

I3CLSimClientModule::~I3CLSimClientModule()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    StopThread();
    
}

void I3CLSimClientModule::StartThread()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    if (threadObj_) {
        log_debug("Thread is already running. Not starting a new one.");
        return;
    }
    
    log_trace("Thread not running. Starting a new one..");
    
    // re-set flags
    threadStarted_=false;
    threadFinishedOK_=false;
    
    threadObj_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&I3CLSimClientModule::Thread_starter, this)));
    
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

void I3CLSimClientModule::StopThread()
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    if (threadObj_)
    {
        if (threadObj_->joinable())
        {
            log_trace("Stopping the thread..");
            // threadObj_->interrupt();
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

void I3CLSimClientModule::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("WorkOnTheseStops", workOnTheseStops_);
    workOnTheseStops_set_ = std::set<I3Frame::Stream>(workOnTheseStops_.begin(), workOnTheseStops_.end());
    
    {
        std::string serverAddress;
        GetParameter("ServerAddress", serverAddress);
        stepsToPhotonsConverter_.reset(new I3CLSimClient(serverAddress));
        GetParameter("StepGenerator", stepGenerator_);
        if (!stepGenerator_)
            log_fatal("You need to provide a step generator");
        uint64_t granularity = stepsToPhotonsConverter_->GetWorkgroupSize();
        uint64_t maxBunchSize = stepsToPhotonsConverter_->GetMaxNumWorkitems();
        stepGenerator_->SetBunchSizeGranularity(granularity);
        stepGenerator_->SetMaxBunchSize(maxBunchSize);
        stepGenerator_->Initialize();
    }

    GetParameter("MCTreeName", MCTreeName_);
    GetParameter("FlasherPulseSeriesName", flasherPulseSeriesName_);
    GetParameter("PhotonSeriesMapName", photonSeriesMapName_);
    GetParameter("MCPESeriesMapName", mcpeSeriesMapName_);
    GetParameter("OMKeyMaskName", omKeyMaskName_);
    GetParameter("IgnoreMuons", ignoreMuons_);
        
    GetParameter("MCPEGenerator", mcpeGenerator_);
    if (mcpeGenerator_ && mcpeSeriesMapName_.empty()) {
        log_fatal("You set MCPEGenerator, but not MCPESeriesMapName");
    } else if (!mcpeGenerator_ && !mcpeSeriesMapName_.empty()) {
        log_fatal("You set MCPESeriesMapName, but not MCPEGenerator");
    }

    GetParameter("StatisticsName", statisticsName_);
    collectStatistics_ = (statisticsName_!="");
    
    GetParameter("IgnoreSubdetectors", ignoreSubdetectors_);
    GetParameter("ClosestDOMDistanceCutoff", closestDOMDistanceCutoff_);
    GetParameter("CosmicEventGenerator",cosmicGenerator_);
    
    if ((flasherPulseSeriesName_=="") && (MCTreeName_==""))
        log_fatal("You need to set at least one of the \"MCTreeName\" and \"FlasherPulseSeriesName\" parameters.");

    StartThread();
    
    log_info("Initialization complete.");
}

namespace {
    static inline ModuleKey ModuleKeyFromOpenCLSimIDs(int16_t stringID, uint16_t domID)
    {
        return ModuleKey(stringID, domID);
    }
}

namespace {

// Compare items in a queue with strictly increasing indices where the index
// may overflow, e.g. UINT32_MAX-2, UINT32_MAX-1, UINT32_MAX, 0, 1
template <typename T>
struct compare_with_overflow
{
    typedef typename std::result_of<decltype(&T::key)(T)>::type Key;
    
    compare_with_overflow(const T& front) : min_(front.key()) {}
    
    bool operator()(const T &a, Key b) { return a.key()-min_ < b-min_; }
    bool operator()(Key a, const T &b) { return a-min_ < b.key()-min_; }
    
    Key min_;
};

template <typename Iterator>
Iterator
lower_bound(Iterator begin, Iterator end, typename compare_with_overflow<typename std::iterator_traits<Iterator>::value_type>::Key value)
{
    typedef typename std::iterator_traits<Iterator>::value_type T;
    typedef typename compare_with_overflow<T>::Key Key;
    i3_assert( std::distance(begin, end) < std::numeric_limits<Key>::max() );
    
    return (begin == end) ? end
        : std::lower_bound(begin, end, value, compare_with_overflow<T>(*begin));
}

template <typename Iterator>
Iterator
upper_bound(Iterator begin, Iterator end, typename compare_with_overflow<typename std::iterator_traits<Iterator>::value_type>::Key value)
{
    typedef typename std::iterator_traits<Iterator>::value_type T;
    typedef typename compare_with_overflow<T>::Key Key;
    i3_assert( std::distance(begin, end) < std::numeric_limits<Key>::max() );
    
    return (begin == end) ? end
        : std::upper_bound(begin, end, value, compare_with_overflow<T>(*begin));
}

template <typename Container, typename Key>
typename Container::iterator
lower_bound(Container &container, Key value)
{
    return ::lower_bound(container.begin(), container.end(), value);
};

template <typename Iterator>
Iterator
find(Iterator begin, Iterator end, typename compare_with_overflow<typename std::iterator_traits<Iterator>::value_type>::Key value)
{
    typedef typename std::iterator_traits<Iterator>::value_type T;
    typedef typename compare_with_overflow<T>::Key Key;
    i3_assert( std::distance(begin, end) < std::numeric_limits<Key>::max() );
    
    if (begin == end)
        return end;
    
    compare_with_overflow<T> compare(*begin);
    Iterator iter = std::lower_bound(begin, end, value, compare);
    return (iter == end || compare(value,*iter)) ? end : iter;
};

template <typename Container, typename Key>
typename Container::iterator
find(Container &container, Key value)
{
    return ::find(container.begin(), container.end(), value);
};

template <typename T>
T Emit(const I3CLSimPhoton &photon, uint32_t currentPhotonId,
    double timeShift, uint32_t particleMinorID, uint64_t particleMajorID);

template <>
I3CompressedPhoton Emit<I3CompressedPhoton>(const I3CLSimPhoton &photon, uint32_t currentPhotonId,
    double timeShift, uint32_t particleMinorID, uint64_t particleMajorID)
{
    I3CompressedPhoton outputPhoton;

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
    
    return outputPhoton;
}

void AddHistoryEntries(const I3CLSimPhotonHistory &photonHistory, I3CompressedPhotonSeries &outputPhotonSeries)
{
    
}

typedef typename I3CLSimClientModule::particleCacheEntry particleCacheEntry;
typedef typename I3CLSimClientModule::frameCacheEntry frameCacheEntry;

void AddPhotonsToFrames(const I3CLSimPhotonSeries &photons,
                        I3CLSimPhotonHistorySeriesConstPtr photonHistories,
                        std::deque<frameCacheEntry> &frameCache,
                        std::deque<particleCacheEntry> &particleCache,
                        I3CLSimPhotonToMCPEConverterPtr &mcpeGenerator,
                        bool collectStatistics_
                        )
{
    if (photonHistories) {
        if (photonHistories->size() != photons.size())
        {
            log_fatal("Error: photon history vector size (%zu) != photon vector size (%zu)",
                      photonHistories->size(), photons.size());
        }
    }
    
    typedef I3CompressedPhotonSeriesMap PhotonSeriesMap;
    typedef typename PhotonSeriesMap::mapped_type PhotonSeries;
    typedef typename PhotonSeries::value_type Photon;
    typedef typename I3CLSimClientModule::particleCacheEntry particleCacheEntry;
    
    size_t nhits = 0;
    size_t nhits_unique = 0;
    for (std::size_t i=0;i<photons.size();++i)
    {
        const I3CLSimPhoton &photon = photons[i];
        
        // find identifier in particle cache
        auto particle = ::find(particleCache, photon.identifier);
        if (particle == particleCache.end())
            log_fatal("Internal error: unknown particle id from OpenCL: %" PRIu32,
                      photon.identifier);
        
        auto frame = ::find(frameCache, particle->frameId);
        if (frame == frameCache.end())
            log_fatal("Internal error: particle cache entry uses invalid frame cache id");
        
        // generate the OMKey
        const ModuleKey key = ModuleKeyFromOpenCLSimIDs(photon.stringID, photon.omID);

        if (frame->ignoreModules.count(key) > 0) continue; // ignore masked DOMs
        
        if (collectStatistics_) {
            // collect statistics
            particle->photonsAtDOMs.count++;
            particle->photonsAtDOMs.weightSum+=photon.GetWeight();
        }
        
        Photon outputPhoton = Emit<Photon>(photon,
                frame->currentPhotonId, particle->timeShift,
                particle->particleMinorID, particle->particleMajorID);
        
        if (frame->photons) {
            // this either inserts a new vector or retrieves an existing one
            PhotonSeries &outputPhotonSeries = (*(frame->photons))[key];
        
            outputPhotonSeries.push_back(outputPhoton);
        
            if (photonHistories) {
                const I3CLSimPhotonHistory &photonHistory = (*photonHistories)[i];
                AddHistoryEntries(photonHistory, outputPhotonSeries);
            }
        
        }
        
        if (frame->hits) {
            OMKey omkey;
            I3MCPE hit;
            bool survived;
            std::tie(omkey, hit, survived) = mcpeGenerator->Convert(key, outputPhoton);
            if (survived) {
                if ((*frame->hits)[omkey].insert(hit))
                    ++nhits_unique;
                ++nhits;
            }
        }
        
        frame->currentPhotonId++;
    }
    log_debug_stream(photons.size()<<" photons -> "<<nhits<<" mcpe ("<<nhits_unique<<" unique)");
}

} // namespace

/**
 * This thread takes care of passing steps from Geant4 to OpenCL
 */

template <typename T>
std::ostream& operator<<(std::ostream &s, const std::deque<T> &vec)
{
    s << "[ ";
    for(const T& v : vec)
        s << v << " ";
    s << "]";
    return s;
}


bool I3CLSimClientModule::Thread(boost::this_thread::disable_interruption &di)
{
    // do some setup while the main thread waits..
    
    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(threadStarted_mutex_);
        threadStarted_=true;
    }
    threadStarted_cond_.notify_all();

    // the main thread is running again
    
    for (;;)
    {
        // retrieve steps from Geant4
        I3CLSimStepSeriesConstPtr steps;
        boost::shared_ptr<std::vector<uint32_t>> finished;
        bool barrierWasJustReset=false;
        
        {
            boost::this_thread::restore_interruption ri(di);
            try {
                std::tie(steps, finished) = stepGenerator_->GetConversionResultWithBarrierInfoAndMarkers(barrierWasJustReset);
            } catch(boost::thread_interrupted &i) {
                return false;
            }
        }
        
        if (!steps) 
        {
            log_debug("Got NULL I3CLSimStepSeriesConstPtr from Geant4.");
            continue;
        }
        else if (steps->empty())
        {
            log_error("Got 0 steps from Geant4, nothing to do for OpenCL.");
            continue;
        }
        else
        {
            log_debug("Got %zu steps from Geant4, sending them to OpenCL",
                     steps->size());
            
            {
                boost::unique_lock<boost::mutex> guard(frameCache_mutex_);
                
                assert( framesInKernel_ >= 0 && framesInKernel_ <= std::distance(frameCache_.begin(), frameCache_.end()) );
                
                // note the frames that have work in this bunch
                {
                    std::set<uint32_t> framesInThisBunch;
                    for (const I3CLSimStep &step : *steps) {
                        if (step.GetNumPhotons() == 0)
                            continue;
                        auto particle = ::find(particleCache_, step.GetID());
                        if (particle != particleCache_.end()) {
                            framesInThisBunch.insert(particle->frameId);
                            // collect statistics if requested
                            if (collectStatistics_) {
                                particle->generatedPhotons.count += step.numPhotons;
                                particle->generatedPhotons.weightSum += static_cast<double>(step.numPhotons)*step.weight;
                            }
                        } else
                            log_fatal("Unknown particle ID %u!", step.GetID());
                    }
                    for (uint32_t frameId : framesInThisBunch) {
                        framesForBunches_.insert(bimap::value_type(frameId, currentBunchId_));
                    }
                }
                
                // find the first frame that could still receive steps
                auto lastFinishedFrame = frameCache_.begin() + framesInKernel_;
                
                // note which light sources were finalized in this bunch
                for (uint32_t identifier : *finished) {
                    auto particle = ::find(particleCache_, identifier);
                    if (particle != particleCache_.end()) {
                        uint32_t frameId = particle->frameId;
                        auto frame = ::find(lastFinishedFrame, frameCache_.end(), frameId);
                        if (frame == frameCache_.end())
                            log_fatal("Unknown frame ID %u!", frameId);
                        assert( frame->numPendingParticles > 0 );
                        // Advance end marker if the frame is final
                        if (--frame->numPendingParticles == 0)
                            lastFinishedFrame = frame;
                    } else
                        log_fatal("Unknown particle ID %u!", identifier);
                }
                
                // At this point, lastFinishedFrame points to the last frame
                // that was finalized in this bunch. All previous frames are
                // finished with step generation.
                {
                    // Frames that contained no work are never associated with
                    // any pending particles. Skip over any that might appear.
                    auto firstPendingFrame = std::find_if(lastFinishedFrame, frameCache_.end(),
                        [](const frameCacheEntry &entry) { return entry.numPendingParticles > 0; });
                    log_debug("%zu frames done with Geant", std::distance(frameCache_.begin() + framesInKernel_, firstPendingFrame));
                    framesInKernel_ = std::distance(frameCache_.begin(), firstPendingFrame);
                }
                
                assert( framesInKernel_ >= 0 && framesInKernel_ <= std::distance(frameCache_.begin(), frameCache_.end()) );
            }


            // send to OpenCL
            {
                boost::this_thread::restore_interruption ri(di);
                try {
                    stepsToPhotonsConverter_->EnqueueSteps(steps, currentBunchId_);
                } catch(boost::thread_interrupted &i) {
                    return false;
                }
            }
            
            ++currentBunchId_; // this may overflow, but only when there are more than 4 billion workitems in flight
        }
        
        // now wait for OpenCL to finish; retrieve results
        
        I3CLSimStepToPhotonConverter::ConversionResult_t res;
        {
            boost::this_thread::restore_interruption ri(di);
            try {
                res = stepsToPhotonsConverter_->GetConversionResult();
            } catch(boost::thread_interrupted &i) {
                return false;
            }
        }
        
        if (!res.photons) log_fatal("Internal error: received NULL photon series from OpenCL.");
        
        {
            boost::unique_lock<boost::mutex> guard(frameCache_mutex_);
            
            // convert to I3Photons and add to their respective frames
            AddPhotonsToFrames(*(res.photons), res.photonHistories,
                               frameCache_,
                               particleCache_,
                               mcpeGenerator_,
                               collectStatistics_);
            
            // This bunch is done. Remove all frame-bunch pairs with this bunch id.
            auto bounds = framesForBunches_.right.equal_range(res.identifier);
            log_debug("bunch %u had work from %zu frames", res.identifier, std::distance(bounds.first, bounds.second));
            
            assert( bounds.first != bounds.second );
            framesForBunches_.right.erase(bounds.first, bounds.second);
            
            newFramesAvailable_ = true;
        }

        if (barrierWasJustReset) {
            log_trace("Geant4 barrier has been reached. Exiting thread.");
            break;
        }

    }
    
    return true;
}

void I3CLSimClientModule::Thread_starter()
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


void I3CLSimClientModule::DigestGeometry(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    if (!(closestDOMDistanceCutoff_ >= 0) || !std::isfinite(closestDOMDistanceCutoff_))
        return;
    
    std::set<std::string> ignoreSubdetectorsSet(ignoreSubdetectors_.begin(), ignoreSubdetectors_.end());
    I3ModuleGeoMapConstPtr moduleGeoMap = frame->Get<I3ModuleGeoMapConstPtr>("I3ModuleGeoMap");
    
    if (!moduleGeoMap) {
        if (frame->Has("I3Geometry")) {
            log_fatal("No I3ModuleGeoMap found in frame. There *is* an I3Geometry object. Please run the \"I3GeometryDecomposer\" module before this.");
        } else {
            log_fatal("No I3ModuleGeoMap found in frame. There does not seem to be a geometry.");
        }
    }

    I3MapModuleKeyStringConstPtr subdetectors = frame->Get<I3MapModuleKeyStringConstPtr>("Subdetectors");
    if (!subdetectors) log_error("No subdetector configuration in frame. Missing a \"Subdetectors\" object. Assuming all modules are on the same detector.");
    
    std::vector<I3Position> doms;
    for(auto &i : *moduleGeoMap) {
        const ModuleKey &key = i.first;
        const I3ModuleGeo &geo = i.second;
        
        std::string subdetectorName = "Unknown"; // use this if there is no Subdetectors object
        if (subdetectors) {
            I3MapModuleKeyString::const_iterator subdetector_it =
            subdetectors->find(key);
            if (subdetector_it == subdetectors->end()) {
                log_fatal("ModuleKey(%i/%u) not found in \"Subdetectors\".",
                          key.GetString(), key.GetOM());
            }
            subdetectorName = subdetector_it->second;
        }
        
        if (ignoreSubdetectorsSet.find(subdetectorName) == ignoreSubdetectorsSet.end())
            doms.push_back(geo.GetPos());
    }
    
    detectorHull_.reset(new I3Surfaces::ExtrudedPolygon(doms, closestDOMDistanceCutoff_));
}

std::size_t I3CLSimClientModule::FlushFrameCache()
{
    // If there are no frames in photon propagation, there can be nothing to flush.
    if (framesInKernel_ == 0 || !newFramesAvailable_)
        return 0;
    
    boost::unique_lock<boost::mutex> guard(frameCache_mutex_);
    
    assert( framesInKernel_ >= 0 && framesInKernel_ <= std::distance(frameCache_.begin(), frameCache_.end()) );
    
    // Any frame in the first framesInKernel_ items but not in framesForBunches_
    // is finished and can be pushed. Conversely, the id of the first frame
    // that is still pending is the first entry in framesForBunches_.
    auto lastPendingFrame = framesForBunches_.empty() ? frameCache_.begin() + framesInKernel_
        : ::lower_bound(frameCache_.begin(), frameCache_.begin() + framesInKernel_, framesForBunches_.left.begin()->first);
    
    // Push frames to outbox
    std::size_t numFrames = std::distance(frameCache_.begin(), lastPendingFrame);
    for (auto frameCacheEntry = frameCache_.begin(); frameCacheEntry != lastPendingFrame; frameCacheEntry++) {
        
        // Commit photons to frame
        if (frameCacheEntry->photons) {
            frameCacheEntry->frame->Put(photonSeriesMapName_, frameCacheEntry->photons);
        }
        // Commit hits
        if (frameCacheEntry->hits) {
            auto hits = boost::make_shared<I3MCPESeriesMap>();
            auto pidmap = boost::make_shared<I3ParticleIDMap>();
            for (auto &pair : *frameCacheEntry->hits) {
                std::tie((*hits)[pair.first], (*pidmap)[pair.first]) = pair.second.extractMCPEsWithPIDInfo();
            }
            frameCacheEntry->frame->Put(mcpeSeriesMapName_, hits);
            frameCacheEntry->frame->Put(mcpeSeriesMapName_+"ParticleIDMap", pidmap);
        }
        
        // Clean up particle cache and collect statistics
        auto &particles = frameCacheEntry->particles;
        if (!particles.empty()) {
            auto begin = ::lower_bound(particleCache_.begin(), particleCache_.end(), particles.front());
            auto end = ::upper_bound(begin, particleCache_.end(), particles.back());
            if (begin == particleCache_.end() || begin == end)
                log_fatal("Could not find particle IDs in the cache");
            if (collectStatistics_) {
                auto statistics = boost::make_shared<I3CLSimEventStatistics>();
                for (auto particle = begin; particle != end; particle++) {
                    statistics->AddNumPhotonsGeneratedWithWeights(particle->generatedPhotons.count,
                                                                  particle->generatedPhotons.weightSum,
                                                                  particle->particleMajorID,
                                                                  particle->particleMinorID);
                    statistics->AddNumPhotonsAtDOMsWithWeights(particle->photonsAtDOMs.count,
                                                               particle->photonsAtDOMs.weightSum,
                                                               particle->particleMajorID,
                                                               particle->particleMinorID);
                }
                frameCacheEntry->frame->Put(statisticsName_, statistics);
            }
            particleCache_.erase(begin, end);
        }
        
        PushFrame(frameCacheEntry->frame);
    }
    frameCache_.erase(frameCache_.begin(), lastPendingFrame);
    framesInKernel_ -= numFrames;
    
    if (numFrames > 0)
        log_debug_stream("Pushed " << numFrames << " frames, " << frameCache_.size() << " remaining");
    
    assert( framesInKernel_ >= 0 && framesInKernel_ <= std::distance(frameCache_.begin(), frameCache_.end()) );
    
    newFramesAvailable_ = false;
    
    return numFrames;
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

bool I3CLSimClientModule::ShouldDoProcess(I3FramePtr frame)
{
    return true;
}

void I3CLSimClientModule::Process()
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
    if ((frameCache_.empty()) && (workOnTheseStops_set_.count(frame->GetStop()) == 0) )
    {
        PushFrame(frame);
        return;
    }
    
    FlushFrameCache();
    DigestOtherFrame(frame);
}

bool I3CLSimClientModule::FrameContainsWork(const I3FramePtr &frame, I3MCTreeConstPtr &MCTree, I3CLSimFlasherPulseSeriesConstPtr &flasherPulses)
{
    // a few cases where we don't work with the frame:
    
    //// not our designated Stop
    if (workOnTheseStops_set_.count(frame->GetStop()) == 0 || !I3ConditionalModule::ShouldDoProcess(frame)) {
        return false;
    }
    
    if (MCTreeName_ != "")
        MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);
    if (flasherPulseSeriesName_ != "")
        flasherPulses = frame->Get<I3CLSimFlasherPulseSeriesConstPtr>(flasherPulseSeriesName_);

    if ((!MCTree) && (!flasherPulses)) {
        // ignore frames without any MCTree and/or Flashers
        return false;
    }
    
    return true;
}

bool I3CLSimClientModule::DigestOtherFrame(I3FramePtr frame)
{
    frameCacheEntry frameCacheEntry(currentFrameId_, frame);
    
    if (!threadObj_)
        log_fatal("No connector thread running");
    
    I3MCTreeConstPtr MCTree;
    I3CLSimFlasherPulseSeriesConstPtr flasherPulses;
    std::deque<I3CLSimLightSource> lightSources;
    std::deque<double> timeOffsets;
    
    if (FrameContainsWork(frame, MCTree, flasherPulses)) {
        I3VectorOMKeyConstPtr omKeyMask;
        I3VectorModuleKeyConstPtr moduleKeyMask;
        if (omKeyMaskName_ != "") {
            omKeyMask = frame->Get<I3VectorOMKeyConstPtr>(omKeyMaskName_);
        
            if (!omKeyMask) {
                moduleKeyMask = frame->Get<I3VectorModuleKeyConstPtr>(omKeyMaskName_);
            }
        }

        if (MCTree){
          I3MCTreePtr newtree = boost::make_shared<I3MCTree>(*MCTree);

          if (cosmicGenerator_){
            CosmicGeneraterToLightSources(*newtree, lightSources, timeOffsets);
          }else{        
            ConvertMCTreeToLightSources(*newtree, lightSources, timeOffsets);
          }

          frame->Delete(MCTreeName_);
          frame->Put(MCTreeName_,newtree);
        }
        if (flasherPulses){
          ConvertFlasherPulsesToLightSources(*flasherPulses, lightSources, timeOffsets);
        }

        // support both vectors of OMKeys and vectors of ModuleKeys
        if (omKeyMask) {
            // assign the current OMKey mask if there is one
            BOOST_FOREACH(const OMKey &key, *omKeyMask) {
                frameCacheEntry.ignoreModules.insert(ModuleKey(key.GetString(), key.GetOM()));
            }
        }
    
        if (moduleKeyMask) {
            // assign the current ModuleKey mask if there is one
            BOOST_FOREACH(const ModuleKey &key, *moduleKeyMask) {
                frameCacheEntry.ignoreModules.insert(key);
            }
        }        
    }

    std::vector<uint32_t> particleIndices;
    {
        boost::unique_lock<boost::mutex> guard(frameCache_mutex_);
        for (std::size_t i=0;i<lightSources.size();++i)
        {
            const I3CLSimLightSource &lightSource = lightSources[i];
            const double timeOffset = timeOffsets[i];
            
            particleCache_.emplace_back(currentParticleCacheIndex_, lightSource, frameCacheEntry.frameId, timeOffset);
            frameCacheEntry.particles.push_back(currentParticleCacheIndex_);
    
            // make a new index. This will eventually overflow,
            // but at that time, index 0 should be unused again.
            ++currentParticleCacheIndex_;
            if (currentParticleCacheIndex_==0) ++currentParticleCacheIndex_; // never use index==0
        }
        
        // make a copy that we can safely reference without holding the mutex
        particleIndices = frameCacheEntry.particles;

        frameCacheEntry.numPendingParticles = frameCacheEntry.particles.size();
        if (!photonSeriesMapName_.empty())
            frameCacheEntry.photons = boost::make_shared<I3CompressedPhotonSeriesMap>();
        if (mcpeGenerator_ && !mcpeSeriesMapName_.empty())
            frameCacheEntry.hits = boost::make_shared<frameCacheEntry::MCPEStreamMap>();
        frameCache_.push_back(std::move(frameCacheEntry));
    
        currentFrameId_++;
    }
    
    assert( particleIndices.size() == lightSources.size() );
    // Do the potentially blocking enqueue outside of the critical section
    for (std::size_t i=0;i<lightSources.size();++i) {
        stepGenerator_->EnqueueLightSource(lightSources[i], particleIndices[i]);
        log_debug_stream("Enqueued "<<i+1<<" of "<<lightSources.size()<<" light sources");
    }

    return true;
}

void I3CLSimClientModule::Finish()
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    log_info("Flushing I3Tray..");
    Flush();
    
    // Flush step generator and wait for collector thread to stop
    stepGenerator_->EnqueueBarrier();
    StopThread();
    
    // Emit any frames still in flight
    FlushFrameCache();
    if (!frameCache_.empty()) {
        log_error("frame cache still has %zu entries", frameCache_.size());
        for(const frameCacheEntry &entry : frameCache_)
            log_error_stream(entry.frameId << ": " << entry.numPendingParticles << " pending of "<<entry.particles);
    }
    if (!particleCache_.empty())
        log_error("particle cache still has %zu entries", particleCache_.size());
    
    // Send them down the chain
    log_info("Flushing I3Tray (again)..");
    Flush();

    log_info("I3CLSimModule is done.");
}




void I3CLSimClientModule::CosmicGeneraterToLightSources(I3MCTree &mcTree,
                                                        std::deque<I3CLSimLightSource> &lightSources,
                                                        std::deque<double> &timeOffsets)
  
{
  const I3Surfaces::ExtrudedPolygon& hull=*detectorHull_;
  auto emitParticle = [&hull,&lightSources,&timeOffsets](I3Particle &p)
                      {
                        log_debug_stream("Got Particle from CosmicGenerator:  "<<p.GetTypeString());
                        const bool isMuon = ((p.GetType() == I3Particle::MuMinus) ||
                                             (p.GetType() == I3Particle::MuPlus));                                            
                        if (isMuon){                                                   
                          double a=0.212/1.2, b=0.251e-3/1.2;
                          double max_range = std::log(1 + p.GetEnergy()*b/a)/b;                          
                          double range =  hull.GetIntersection(p.GetPos(), p.GetDir()).first;
                          if (range < max_range){                            
                            const double particleTime = p.GetTime();
                            p.SetLocationType(I3Particle::InIce);
                            lightSources.push_back(I3CLSimLightSource(p));
                            timeOffsets.push_back(0);
                          }
                        }
                      };
  I3Frame dummy_frame;          
  cosmicGenerator_->Generate(mcTree, dummy_frame, emitParticle);
}

//////////////

void I3CLSimClientModule::ConvertMCTreeToLightSources(const I3MCTree &mcTree,
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


void I3CLSimClientModule::ConvertFlasherPulsesToLightSources(const I3CLSimFlasherPulseSeries &flasherPulses,
                                                       std::deque<I3CLSimLightSource> &lightSources,
                                                       std::deque<double> &timeOffsets)
{
    BOOST_FOREACH(const I3CLSimFlasherPulse &flasherPulse, flasherPulses)
    {
        lightSources.push_back(I3CLSimLightSource(flasherPulse));
        timeOffsets.push_back(0.);
    }
}
