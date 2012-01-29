#ifndef I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED
#define I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED

#include "clsim/I3CLSimStepToPhotonConverter.h"

#include "phys-services/I3RandomService.h"

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "clsim/I3CLSimQueue.h"

#include "clsim/I3CLSimOpenCLDevice.h"

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
    void SetDisableDoubleBuffering(bool value);

    /**
     * Returns true if double buffering is disabled.
     */
    bool GetDisableDoubleBuffering() const;

    /**
     * Sets the wavelength generator. By default it should
     * return wavelengths according to a Cherenkov spectrum.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenGenerator(I3CLSimRandomValueConstPtr wlenGenerator);
    
    /**
     * Sets the wavelength weights. Set this to a constant value
     * of 1 if you do not need biased photon generation.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias);

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

    
    boost::shared_ptr<boost::thread> openCLThreadObj_;
    boost::condition_variable_any openCLStarted_cond_;
    boost::mutex openCLStarted_mutex_;
    bool openCLStarted_;
    
    boost::shared_ptr<I3CLSimQueue<ToOpenCLPair_t> > queueToOpenCL_;
    boost::shared_ptr<I3CLSimQueue<I3CLSimStepToPhotonConverter::ConversionResult_t> > queueFromOpenCL_;

    I3RandomServicePtr randomService_;
    
    bool initialized_;
    bool compiled_;
    I3CLSimRandomValueConstPtr wlenGenerator_;
    I3CLSimWlenDependentValueConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    I3CLSimSimpleGeometryConstPtr geometry_;
    
    I3CLSimOpenCLDevicePtr device_;
    bool useNativeMath_;
    std::size_t selectedDeviceIndex_;
    bool deviceIsSelected_;
    
    bool disableDoubleBuffering_;
    
    // some kernel sources loaded on construction
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
    std::vector<shared_ptr<cl::CommandQueue> > queue_;
    std::vector<shared_ptr<cl::Kernel> > kernel_;
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
    std::vector<shared_ptr<cl::Buffer> > deviceBuffer_InputSteps;
    std::vector<shared_ptr<cl::Buffer> > deviceBuffer_OutputPhotons;
    std::vector<shared_ptr<cl::Buffer> > deviceBuffer_CurrentNumOutputPhotons;
    
    // this one is constant, so we only need one
    shared_ptr<cl::Buffer> deviceBuffer_GeoLayerToOMNumIndexPerStringSet;
    
    // Size of output photon storage (maximum amount of photons per step bunch)
    uint32_t maxNumOutputPhotons_;
    
    
};

I3_POINTER_TYPEDEFS(I3CLSimStepToPhotonConverterOpenCL);

#endif //I3CLSIMSTEPTOPHOTONCONVERTEROPENCL_H_INCLUDED
