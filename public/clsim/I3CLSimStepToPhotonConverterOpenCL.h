#ifndef I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED
#define I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED

#include "clsim/I3CLSimStepToPhotonConverter.h"

#include "phys-services/I3RandomService.h"

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "clsim/I3CLSimQueue.h"

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
};

/**
 * @brief Creates photons from a given list of steps and propagates
 * them to a DOM using an OpenCL-enabled algorithm
 */
struct I3CLSimStepToPhotonConverterOpenCL : public I3CLSimStepToPhotonConverter
{
public:
    static const bool default_useNativeMath;
    static const bool default_cpuOnly;
    static const bool default_gpuOnly;
    
    I3CLSimStepToPhotonConverterOpenCL(I3RandomServicePtr randomService,
                                       bool useNativeMath=default_useNativeMath,
                                       bool cpuOnly=default_cpuOnly,
                                       bool gpuOnly=default_gpuOnly);
    virtual ~I3CLSimStepToPhotonConverterOpenCL();

    /**
     * Gets a list of all available OpenCL devices as
     * <platform, device> pairs.
     *
     * Will throw if already initialized.
     */
    shared_ptr<const std::vector<std::pair<std::string, std::string> > > GetDeviceList() const;

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
     * Sets the device by index.
     * The index should be chosen from the list
     * returned by GetDeviceList().
     *
     * Will throw if already initialized.
     */
    void SetDeviceIndex(std::size_t selectedDeviceIndex);

    /**
     * Sets the device by platform/device name.
     * The name should be chosen from the list
     * returned by GetDeviceList().
     *
     * Will throw if already initialized.
     */
    void SetDeviceName(const std::string &platform, const std::string &device);
    
    /**
     * Returns the maximum workgroup size for the
     * current kernel.
     *
     * Will throw if compilation fails.
     * Will throw if initialized.
     */
    uint64_t GetMaxWorkgroupSize() const;
    
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
    
private:
    typedef std::pair<uint32_t, I3CLSimStepSeriesConstPtr> ToOpenCLPair_t;

    // sets up OpenCL
    void SetupQueueAndKernel(const cl::Platform& platform, const cl::Device &device);

    
    void OpenCLThread();
    void OpenCLThread_impl(boost::this_thread::disable_interruption &di);
    boost::shared_ptr<boost::thread> openCLThreadObj_;
    boost::condition_variable_any openCLStarted_cond_;
    boost::mutex openCLStarted_mutex_;
    bool openCLStarted_;
    
    boost::shared_ptr<I3CLSimQueue<ToOpenCLPair_t> > queueToOpenCL_;
    boost::shared_ptr<I3CLSimQueue<I3CLSimStepToPhotonConverter::ConversionResult_t> > queueFromOpenCL_;

    I3RandomServicePtr randomService_;
    
    bool initialized_;
    bool compiled_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    I3CLSimSimpleGeometryConstPtr geometry_;
    
    shared_ptr<const std::vector<std::pair<std::string, std::string> > > deviceNameList_;
    std::vector<std::pair<shared_ptr<cl::Platform>, shared_ptr<cl::Device> > > clPlatformDeviceList_;
    bool cpuOnly_;
    bool gpuOnly_;
    bool useNativeMath_;
    std::size_t selectedDeviceIndex_;
    bool deviceIndexIsSelected_;
    
    // some kernel sources loaded on construction
    std::string mwcrngKernelSource_;
    std::string mediumPropertiesSource_;
    std::string geometrySource_;
    std::string propagationKernelSource_;
    
    // this is extra geometry information, we upload it to global memory
    std::vector<unsigned char> geoLayerToOMNumIndexPerStringSetInfo_;
    
    // this allows us to convert the string index back to the string ID (which may be negative and non-contiguous)
    std::vector<int> stringIndexToStringIDBuffer_;
    
    // OpenCL command queue and kernel
    shared_ptr<cl::CommandQueue> queue_;
    shared_ptr<cl::Kernel> kernel_;
    shared_ptr<cl::Context> context_;
    
    // maximum workgroup size for current kernel
    uint64_t maxWorkgroupSize_;
    
    // configured workgroup size and maximum number of work items
    std::size_t workgroupSize_;
    std::size_t maxNumWorkitems_;
    
    // rng state per workitem
    std::vector<uint64_t> MWC_RNG_x;
    std::vector<uint32_t> MWC_RNG_a;
    shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_x;
    shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_a;
    
    // Memory buffers on the device
    shared_ptr<cl::Buffer> deviceBuffer_InputSteps;
    shared_ptr<cl::Buffer> deviceBuffer_OutputPhotons;
    shared_ptr<cl::Buffer> deviceBuffer_CurrentNumOutputPhotons;
    shared_ptr<cl::Buffer> deviceBuffer_GeoLayerToOMNumIndexPerStringSet;
    
    // Size of output photon storage (maximum amount of photons per step bunch)
    uint32_t maxNumOutputPhotons_;
    
    
};

I3_POINTER_TYPEDEFS(I3CLSimStepToPhotonConverter);

#endif //I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED
