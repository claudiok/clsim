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
 * @file I3CLSimStepToPhotonConverterOpenCL.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED
#define I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED

#include "clsim/I3CLSimStepToPhotonConverter.h"

#include "phys-services/I3RandomService.h"

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "clsim/I3CLSimQueue.h"

#include "clsim/I3CLSimOpenCLDevice.h"

#include <vector>
#include <map>
#include <string>
#include <stdexcept>

// forward declarations
namespace cl {
    class Device;
    class Platform;
    class CommandQueue;
    class Context;
    class Kernel;
    class Buffer;
    class Event;
};

/**
 * @brief Creates photons from a given list of steps and propagates
 * them to a DOM using an OpenCL-enabled algorithm
 */
struct I3CLSimStepToPhotonConverterOpenCL : public I3CLSimStepToPhotonConverter
{
public:
    static const bool default_useNativeMath;
    
    I3CLSimStepToPhotonConverterOpenCL(I3RandomServicePtr randomService,
                                       bool useNativeMath=default_useNativeMath);
    virtual ~I3CLSimStepToPhotonConverterOpenCL();

    /**
     * Sets the workgroup size. A value of 0 
     * uses the maximum possible workgroup size
     * for the kernel.
     *
     * Will throw if already initialized.
     */
    void SetWorkgroupSize(std::size_t val);

    /**
     * Sets the number of parallel work items.
     *
     * Will throw if already initialized.
     */
    void SetMaxNumWorkitems(std::size_t val);

    /**
     * Gets the current workgroup size.
     */
    std::size_t GetWorkgroupSize() const;
    
    /**
     * Gets the number of parallel work items.
     */
    std::size_t GetMaxNumWorkitems() const;
    
    /**
     * Sets the OpenCL device.
     *
     * Will throw if already initialized.
     */
    void SetDevice(const I3CLSimOpenCLDevice &device);
    
    /**
     * Returns the maximum workgroup size for the
     * current kernel.
     *
     * Will throw if compilation fails.
     * Will throw if initialized.
     */
    uint64_t GetMaxWorkgroupSize() const;

    /**
     * Disables or enables double-buffered
     * GPU usage. Double buffering will use
     * two command queues and two sets of input
     * and output buffers in order to transfer
     * data to the GPU while a kernel is executing
     * on the other buffer.
     *
     * This has been observed to yield empty results
     * results on older drivers for the nVidia
     * architecture, so it is disabled by default.
     *
     * Before enabling this for a certain driver/hardware
     * combination, make sure that both correct results
     * are returned. Most of the time the second buffer
     * results are always empty, so this error should be
     * easy to observe.
     *
     * Will throw if already initialized.
     */
    void SetEnableDoubleBuffering(bool value);

    /**
     * Returns true if double buffering is enabled.
     */
    bool GetEnableDoubleBuffering() const;

    /**
     * Enables double-precision support in the
     * kernel. This slows down calculations and
     * requires more memory.
     *
     * The performance hit is minimal on CPUs
     * but up to an order of magnitude on GPUs.
     *
     * Will throw if already initialized.
     */
    void SetDoublePrecision(bool value);
        
    /**
     * Returns true if double precision is enabled.
     */
    bool GetDoublePrecision() const;

    /**
     * Configures behaviour for photons that
     * hit a DOM. If this is true (the default)
     * photons will be stopped once they hit a
     * DOM. If this is false, they continue to
     * propagate.
     *
     * Will throw if already initialized.
     */
    void SetStopDetectedPhotons(bool value);
    
    /**
     * Returns true if detected photons are stopped.
     */
    bool GetStopDetectedPhotons() const;

    /**
     * Tell the propagator to save all photons,
     * regardless of detection. All photons will
     * be assigned to DOM(0,0).
     *
     * Will throw if already initialized.
     */
    void SetSaveAllPhotons(bool value);
    
    /**
     * Returns true if all photons are saved,
     * regardless of detection.
     */
    bool GetSaveAllPhotons() const;

    /**
     * Sets the prescale factor of photons
     * being generated in "saveAllPhotons" mode.
     * Only this fraction of photons is actually
     * generated.
     *
     * Will throw if already initialized.
     */
    void SetSaveAllPhotonsPrescale(double value);
    
    /**
     * Returns true if all photons are saved,
     * regardless of detection.
     */
    double GetSaveAllPhotonsPrescale() const;

    /**
     * Sets the maximum number of entries in the photon
     * history table. Each point in the table
     * will store the position of the photon
     * at each point of scatter. (Only the most
     * recent points are stored if there are
     * more scatters than available entries.)
     *
     * Will throw if already initialized.
     */
    void SetPhotonHistoryEntries(uint32_t value);
    
    /**
     * Returns the maximum number of photon
     * history entries.
     */
    uint32_t GetPhotonHistoryEntries() const;

    /**
     * Sets the number of absorption lengths each photon
     * should be propagated. If set to NaN (the default),
     * the number is sampled from an exponential distribution.
     * Use this override for table-making.
     *
     * Will throw if already initialized.
     */
    void SetFixedNumberOfAbsorptionLengths(double value);
    
    /**
     * Returns number of absorption lengths each photon
     * should be propagated. If set to NaN (the default),
     * the number is sampled from an exponential distribution.
     * Use this override for table-making.
     */
    double GetFixedNumberOfAbsorptionLengths() const;

    /**
     * Sets the "pancake" factor for DOMs. For standard
     * oversized-DOM simulations, this should be the
     * radius oversizing factor. This will flatten the
     * DOM in the direction parallel to the photon.
     * The DOM will have a pancake-like shape, elongated
     * in the directions perpendicular to the photon direction.
     *
     * The DOM radius (supplied by the geometry) must also include
     * the oversizing factor.
     *
     * Will throw if already initialized.
     */
    void SetDOMPancakeFactor(double value);
    
    /**
     * Returns the "pancake" factor for DOMs.
     */
    double GetDOMPancakeFactor() const;

    /**
     * Sets the wavelength generators. 
     * The first generator (index 0) is assumed to return a Cherenkov
     * spectrum that may have a bias applied to it. This bias factor
     * needs to be set using SetWlenBias().
     * All other generator indices are assumed to be for flasher/laser
     * light generation. During generation, no Cherenkov angle
     * rotation will be applied to those photons with indices >= 1.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenGenerators(const std::vector<I3CLSimRandomValueConstPtr> &wlenGenerators);
    
    /**
     * Sets the wavelength weights. Set this to a constant value
     * of 1 if you do not need biased photon generation.
     * The wavelength spectrum set with SetWlenGenerator()
     * is assumed to have a biassing factor already applied to it.
     * This call sets this factor in order to be able to assign
     * correct weights.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias);

    /**
     * Sets the medium properties.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties);
    
    /**
     * Sets the geometry.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetGeometry(I3CLSimSimpleGeometryConstPtr geometry);

    /**
     * Compiles the kernel. Can only be used
     * after medium properties, geometry and device have
     * been set.
     * 
     * Will throw if already initialized.
     */
    virtual void Compile();

    /**
     * Gets the full generated OpenCL source code. 
     *
     * Will throw if not yet compiled.
     */
    virtual std::string GetFullSource();
    
    /**
     * Provide individual bits of the OpenCL source code. 
     */
    virtual std::string GetPreambleSource();
    std::string GetMediumPropertiesSource();
    std::string GetWlenGeneratorSource();
    std::string GetWlenBiasSource();
    virtual std::string GetGeometrySource();
    virtual std::string GetCollisionDetectionSource(bool header=true);
    
    /**
     * Initializes the simulation.
     * Will throw if already initialized.
     */
    virtual void Initialize();

    /**
     * Returns true if initialized.
     * Never throws.
     */
    virtual bool IsInitialized() const;
    
    /**
     * Adds a new I3CLSimStepSeries to the queue.
     * The resulting I3CLSimPhotonSeries can be retrieved from the
     * I3CLSimStepToPhotonConverter after some processing time.
     *
     * Enqueuing a vector after calling EnqueueBarrier 
     * will throw if not all photons have been retrieved.
     *
     * Will throw if not initialized.
     */
    virtual void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier);

    /**
     * Reports the current queue size. The queue works asynchronously,
     * so this value will probably have changed once you use it.
     *
     * Will throw if not initialized.
     */
    virtual std::size_t QueueSize() const; 

    /**
     * Returns true if more photons are available.
     * If the return value is false, the current simulation is finished
     * and a new step vector may be set.
     * 
     * Will throw if not initialized.
     */
    virtual bool MorePhotonsAvailable() const;

    /**
     * Returns a bunch of photons stored in a vector<I3CLSimPhoton>.
     *
     * The return value is a pair<uint, vector<I3CLSimPhoton> >.
     * The integer is the same identifier as specified in the call
     * to EnqueueSteps().
     *
     * Might block if no photons are available.
     * 
     * Will throw if not initialized.
     */
    virtual I3CLSimStepToPhotonConverter::ConversionResult_t GetConversionResult();
    
    virtual std::map<std::string, double> GetStatistics() const;
    
    inline double GetTotalDeviceTime() const {boost::unique_lock<boost::mutex> guard(statistics_mutex_); return static_cast<double>(statistics_total_device_duration_in_nanoseconds_);}
    inline double GetTotalHostTime() const {boost::unique_lock<boost::mutex> guard(statistics_mutex_); return static_cast<double>(statistics_total_host_duration_in_nanoseconds_);}
    inline uint64_t GetNumKernelCalls() const {boost::unique_lock<boost::mutex> guard(statistics_mutex_); return statistics_total_kernel_calls_;}
    inline uint64_t GetTotalNumPhotonsGenerated() const {boost::unique_lock<boost::mutex> guard(statistics_mutex_); return statistics_total_num_photons_generated_;}
    inline uint64_t GetTotalNumPhotonsAtDOMs() const {boost::unique_lock<boost::mutex> guard(statistics_mutex_); return statistics_total_num_photons_atDOMs_;}
    
private:
    typedef std::pair<uint32_t, I3CLSimStepSeriesConstPtr> ToOpenCLPair_t;

    // sets up OpenCL
    void SetupQueueAndKernel(const cl::Platform& platform, const cl::Device &device);

    
    void OpenCLThread();
    void OpenCLThread_impl(boost::this_thread::disable_interruption &di);
    bool OpenCLThread_impl_uploadSteps(boost::this_thread::disable_interruption &di,
                                       bool &shouldBreak,
                                       unsigned int bufferIndex,
                                       uint32_t &out_stepsIdentifier,
                                       uint64_t &out_totalNumberOfPhotons,
                                       std::size_t &out_numberOfInputSteps,
                                       bool blocking=true
                                       );
    void OpenCLThread_impl_downloadPhotons(boost::this_thread::disable_interruption &di,
                                           bool &shouldBreak,
                                           unsigned int bufferIndex,
                                           uint32_t stepsIdentifier);
    void OpenCLThread_impl_runKernel(unsigned int bufferIndex,
                                     cl::Event &kernelFinishEvent,
                                     std::size_t numberOfInputSteps);

    boost::posix_time::ptime DumpStatistics(const cl::Event &kernelFinishEvent,
                                            const boost::posix_time::ptime &last_timestamp,
                                            uint64_t totalNumberOfPhotons,
                                            bool starving,
                                            const std::string &platformName,
                                            const std::string &deviceName,
                                            uint64_t deviceProfilingResolution);

    mutable boost::mutex statistics_mutex_;
    uint64_t statistics_total_device_duration_in_nanoseconds_;
    uint64_t statistics_total_host_duration_in_nanoseconds_;
    uint64_t statistics_total_kernel_calls_;
    uint64_t statistics_total_num_photons_generated_;
    uint64_t statistics_total_num_photons_atDOMs_;

    
    boost::shared_ptr<boost::thread> openCLThreadObj_;
    boost::condition_variable_any openCLStarted_cond_;
    boost::mutex openCLStarted_mutex_;
    bool openCLStarted_;
    
    boost::shared_ptr<I3CLSimQueue<ToOpenCLPair_t> > queueToOpenCL_;
    boost::shared_ptr<I3CLSimQueue<I3CLSimStepToPhotonConverter::ConversionResult_t> > queueFromOpenCL_;

    I3RandomServicePtr randomService_;
    
    bool initialized_;
    bool compiled_;
    std::vector<I3CLSimRandomValueConstPtr> wlenGenerators_;
    I3CLSimFunctionConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    I3CLSimSimpleGeometryConstPtr geometry_;
    
    I3CLSimOpenCLDevicePtr device_;
    bool useNativeMath_;
    bool deviceIsSelected_;
    
    bool disableDoubleBuffering_;
    bool doublePrecision_;
    bool stopDetectedPhotons_;
    bool saveAllPhotons_;
    double saveAllPhotonsPrescale_;
    double fixedNumberOfAbsorptionLengths_;
    double pancakeFactor_;
    
    uint32_t photonHistoryEntries_;
    
    // some kernel sources loaded on construction
    std::string prependSource_;
    std::string mwcrngKernelSource_;
    std::string wlenGeneratorSource_;
    std::string wlenBiasSource_;
    std::string mediumPropertiesSource_;
    std::string geometrySource_;
    std::string propagationKernelSource_;
    
    // this is extra geometry information, we upload it to global memory
    std::vector<unsigned short> geoLayerToOMNumIndexPerStringSetInfo_;
    
    // this allows us to convert the string index back to the string ID (which may be negative and non-contiguous)
    std::vector<int> stringIndexToStringIDBuffer_;

    // this allows us to convert the DOM index back to the DOM ID (which may be non-contiguous)
    std::vector<std::vector<unsigned int> > domIndexToDomIDBuffer_perStringIndex_;
    
    // OpenCL command queue and kernel
    std::vector<boost::shared_ptr<cl::CommandQueue> > queue_;
    std::vector<boost::shared_ptr<cl::Kernel> > kernel_;
    boost::shared_ptr<cl::Context> context_;
    
    // maximum workgroup size for current kernel
    uint64_t maxWorkgroupSize_;
    
    // configured workgroup size and maximum number of work items
    std::size_t workgroupSize_;
    std::size_t maxNumWorkitems_;
    
    // rng state per workitem
    std::vector<uint64_t> MWC_RNG_x;
    std::vector<uint32_t> MWC_RNG_a;
    boost::shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_x;
    boost::shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_a;
    
    // Memory buffers on the device
    std::vector<boost::shared_ptr<cl::Buffer> > deviceBuffer_InputSteps;
    std::vector<boost::shared_ptr<cl::Buffer> > deviceBuffer_OutputPhotons;
    std::vector<boost::shared_ptr<cl::Buffer> > deviceBuffer_CurrentNumOutputPhotons;
    std::vector<boost::shared_ptr<cl::Buffer> > deviceBuffer_PhotonHistory;
    
    // this one is constant, so we only need one
    boost::shared_ptr<cl::Buffer> deviceBuffer_GeoLayerToOMNumIndexPerStringSet;
    
    // Size of output photon storage (maximum amount of photons per step bunch)
    uint32_t maxNumOutputPhotons_;
    
    SET_LOGGER("I3CLSimStepToPhotonConverterOpenCL");
};

I3_POINTER_TYPEDEFS(I3CLSimStepToPhotonConverterOpenCL);

#endif //I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED
