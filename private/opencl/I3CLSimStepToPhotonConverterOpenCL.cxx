#define __STDC_FORMAT_MACROS 
#include <inttypes.h>

#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"

// debugging: show GPUtime/photon
#define DUMP_STATISTICS


#include <string>
#include <sstream>
#include <algorithm>
#include <limits>

#include <stdlib.h>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"
#include "opencl/I3CLSimHelperGenerateGeometrySource.h"

#include "opencl/mwcrng_init.h"

#define __CL_ENABLE_EXCEPTIONS
#include "clsim/cl.hpp"

const bool I3CLSimStepToPhotonConverterOpenCL::default_useNativeMath=true;


I3CLSimStepToPhotonConverterOpenCL::I3CLSimStepToPhotonConverterOpenCL(I3RandomServicePtr randomService,
                                                                       bool useNativeMath)
:
openCLStarted_(false),
queueToOpenCL_(new I3CLSimQueue<ToOpenCLPair_t>(5)),
queueFromOpenCL_(new I3CLSimQueue<I3CLSimStepToPhotonConverter::ConversionResult_t>(0)),
randomService_(randomService),
initialized_(false),
compiled_(false),
useNativeMath_(useNativeMath),
selectedDeviceIndex_(0),
deviceIsSelected_(false),
maxWorkgroupSize_(0),
workgroupSize_(0),
maxNumWorkitems_(10240)
{
    if (!randomService_) log_fatal("You need to supply a I3RandomService.");
    
    // load program source from files
    const std::string I3_SRC(getenv("I3_SRC"));
    const std::string kernelBaseDir = I3_SRC+"/clsim/resources/kernels";
    
    try {
        mwcrngKernelSource_ = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/mwcrng_kernel.cl");
        propagationKernelSource_ = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/propagation_kernel.cl");
    } catch (std::runtime_error &e) {
        throw I3CLSimStepToPhotonConverter_exception((std::string("Could not load kernel: ") + e.what()).c_str());
    }
    
}

I3CLSimStepToPhotonConverterOpenCL::~I3CLSimStepToPhotonConverterOpenCL()
{
    if (openCLThreadObj_)
    {
        if (openCLThreadObj_->joinable())
        {
            log_debug("Stopping the OpenCL worker thread..");
            
            openCLThreadObj_->interrupt();
            
            openCLThreadObj_->join(); // wait for it indefinitely
            
            log_debug("OpenCL worker thread stopped.");
        }
        
        openCLThreadObj_.reset();
    }
    
    // reset buffers
    deviceBuffer_MWC_RNG_x.reset();
    deviceBuffer_MWC_RNG_a.reset();
    
    deviceBuffer_InputSteps.clear();
    deviceBuffer_OutputPhotons.clear();
    deviceBuffer_CurrentNumOutputPhotons.clear();

    deviceBuffer_GeoLayerToOMNumIndexPerStringSet.reset();
    
    // reset pointers
    compiled_=false;
    context_.reset();
    kernel_.clear();
    queue_.clear();
    
}

uint64_t I3CLSimStepToPhotonConverterOpenCL::GetMaxWorkgroupSize() const
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    if (!compiled_)
        throw I3CLSimStepToPhotonConverter_exception("You need to compile the kernel first. Call Compile().");
    
    return maxWorkgroupSize_;
}

void I3CLSimStepToPhotonConverterOpenCL::SetDevice(const I3CLSimOpenCLDevice &device)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    if (device_)
    {
        if (!(*device_ == device))
        {
            compiled_=false;
            kernel_.clear();
            queue_.clear();
            device_.reset();
        }
    }
    
    device_ = I3CLSimOpenCLDevicePtr(new I3CLSimOpenCLDevice(device)); // make a copy of the device
    
    deviceIsSelected_=true;
}

void I3CLSimStepToPhotonConverterOpenCL::SetWorkgroupSize(std::size_t val)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    workgroupSize_=val;
}

void I3CLSimStepToPhotonConverterOpenCL::SetMaxNumWorkitems(std::size_t val)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    if (val <= 0)
        throw I3CLSimStepToPhotonConverter_exception("Invalid maximum number of work items!");
    
    maxNumWorkitems_=val;
}

std::size_t I3CLSimStepToPhotonConverterOpenCL::GetWorkgroupSize() const
{
    if (workgroupSize_==0)
    {
        if (!compiled_)
            throw I3CLSimStepToPhotonConverter_exception("Automatic workgroup size cannot be returned before Compile() has been called!");
        
        return maxWorkgroupSize_;
    }
    
    return workgroupSize_;
}

std::size_t I3CLSimStepToPhotonConverterOpenCL::GetMaxNumWorkitems() const
{
    return maxNumWorkitems_;
}


void I3CLSimStepToPhotonConverterOpenCL::Initialize()
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    log_info("Setting up OpenCL..");
    
    Compile();
    
    if (workgroupSize_>maxWorkgroupSize_)
        throw I3CLSimStepToPhotonConverter_exception("Workgroup size too large!");
    
    log_debug("basic OpenCL setup done.");
    
    // start with a maximum number of output photons of the same size as the number of
    // input steps. Should be plenty..
    maxNumOutputPhotons_ = static_cast<uint32_t>(std::min(maxNumWorkitems_, static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())));
    
    // set up rng
    log_debug("Setting up RNG for %zu workitems.", maxNumWorkitems_);
    
    MWC_RNG_x.resize(maxNumWorkitems_);
    MWC_RNG_a.resize(maxNumWorkitems_);
    
    const std::string I3_SRC(getenv("I3_SRC"));
    if (init_MWC_RNG(&(MWC_RNG_x[0]), &(MWC_RNG_a[0]), maxNumWorkitems_, (I3_SRC+"/clsim/resources/safeprimes_base32.txt").c_str(), randomService_)!=0) 
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    log_debug("RNG is set up..");
    
    log_debug("Setting up device buffers..");
    
    // reset all buffers first
    deviceBuffer_MWC_RNG_x.reset();
    deviceBuffer_MWC_RNG_a.reset();
    deviceBuffer_InputSteps.clear();
    deviceBuffer_OutputPhotons.clear();
    deviceBuffer_CurrentNumOutputPhotons.clear();
    deviceBuffer_GeoLayerToOMNumIndexPerStringSet.reset();
    
    
    // set up device buffers from existing host buffers
    deviceBuffer_MWC_RNG_x = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_x.size() * sizeof(uint64_t), &(MWC_RNG_x[0])));
    
    deviceBuffer_MWC_RNG_a = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_a.size() * sizeof(uint32_t), &(MWC_RNG_a[0])));
    
    deviceBuffer_GeoLayerToOMNumIndexPerStringSet = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, geoLayerToOMNumIndexPerStringSetInfo_.size() * sizeof(unsigned short), &(geoLayerToOMNumIndexPerStringSetInfo_[0])));
    
    
    // allocate empty buffers on the device
    for (unsigned int i=0;i<2;++i)
    {
        deviceBuffer_InputSteps.push_back(shared_ptr<cl::Buffer>
        (new cl::Buffer(*context_, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, maxNumWorkitems_*sizeof(I3CLSimStep), NULL)));
        
        deviceBuffer_OutputPhotons.push_back(shared_ptr<cl::Buffer>
        (new cl::Buffer(*context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, maxNumOutputPhotons_*sizeof(I3CLSimPhoton), NULL)));
        
        deviceBuffer_CurrentNumOutputPhotons.push_back(shared_ptr<cl::Buffer>
        (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(uint32_t), NULL)));
    }
    
    log_debug("Device buffers are set up.");
    
    log_debug("Configuring kernel.");
    for (unsigned int i=0;i<2;++i)
    {
        kernel_[i]->setArg(0, *(deviceBuffer_CurrentNumOutputPhotons[i]));     // hit counter
        kernel_[i]->setArg(1, maxNumOutputPhotons_);                           // maximum number of possible hits
        
        kernel_[i]->setArg(2, *deviceBuffer_GeoLayerToOMNumIndexPerStringSet); // additional geometry information (did not fit into constant memory)
        
        kernel_[i]->setArg(3, *(deviceBuffer_InputSteps[i]));                  // the input steps
        kernel_[i]->setArg(4, *(deviceBuffer_OutputPhotons[i]));               // the output photons
        
        kernel_[i]->setArg(5, *deviceBuffer_MWC_RNG_x);                        // rng state
        kernel_[i]->setArg(6, *deviceBuffer_MWC_RNG_a);                        // rng state
    }
    log_debug("Kernel configured.");
    
    log_debug("Starting the OpenCL worker thread..");
    openCLStarted_=false;
    
    openCLThreadObj_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&I3CLSimStepToPhotonConverterOpenCL::OpenCLThread, this)));
    
    // wait for startup
    {
        boost::unique_lock<boost::mutex> guard(openCLStarted_mutex_);
        for (;;)
        {
            if (openCLStarted_) break;
            openCLStarted_cond_.wait(guard);
        }
    }        
    
    log_debug("OpenCL worker thread started.");
    
    log_info("OpenCL setup complete.");
    
    initialized_=true;
}


void I3CLSimStepToPhotonConverterOpenCL::Compile()
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    if (compiled_) return; // silently
    
    if (!wlenGenerator_)
        throw I3CLSimStepToPhotonConverter_exception("WlenGenerator not set!");
    
    if (!wlenBias_)
        throw I3CLSimStepToPhotonConverter_exception("WlenBias not set!");
    
    if (!mediumProperties_)
        throw I3CLSimStepToPhotonConverter_exception("MediumProperties not set!");
    
    if (!geometry_)
        throw I3CLSimStepToPhotonConverter_exception("Geometry not set!");
    
    if (!deviceIsSelected_)
        throw I3CLSimStepToPhotonConverter_exception("Device not selected!");
    
    wlenGeneratorSource_ = wlenGenerator_->GetOpenCLFunction("generateWavelength", // name
                                                             // these are all defined as macros by the rng code:
                                                             "RNG_ARGS",               // function arguments for rng
                                                             "RNG_ARGS_TO_CALL",       // if we call anothor function, this is how we pass on the rng state
                                                             "RNG_CALL_UNIFORM_CO",    // the call to the rng for creating a uniform number [0;1[
                                                             "RNG_CALL_UNIFORM_OC"     // the call to the rng for creating a uniform number ]0;1]
                                                             );
    wlenBiasSource_ = wlenBias_->GetOpenCLFunction("getWavelengthBias"); // name
    
    mediumPropertiesSource_ = I3CLSimHelper::GenerateMediumPropertiesSource(*mediumProperties_);
    geometrySource_ = I3CLSimHelper::GenerateGeometrySource(*geometry_,
                                                            geoLayerToOMNumIndexPerStringSetInfo_,
                                                            stringIndexToStringIDBuffer_,
                                                            domIndexToDomIDBuffer_perStringIndex_);
    
    SetupQueueAndKernel(*(device_->GetPlatformHandle()),
                        *(device_->GetDeviceHandle()));
    
    compiled_=true;
}

std::string I3CLSimStepToPhotonConverterOpenCL::GetFullSource()
{
    std::ostringstream code;
    
    code << mwcrngKernelSource_;
    code << wlenGeneratorSource_;
    code << wlenBiasSource_;
    code << mediumPropertiesSource_;
    code << geometrySource_;
    code << propagationKernelSource_;
    
    return code.str();
}

void I3CLSimStepToPhotonConverterOpenCL::SetupQueueAndKernel(const cl::Platform &platform,
                                                             const cl::Device &device)
{
    VECTOR_CLASS<cl::Device> devices(1, device);
    
    // prepare a device vector (containing a single device)
    cl_context_properties properties[] = 
    { CL_CONTEXT_PLATFORM, (cl_context_properties)(platform)(), 0};
    
    unsigned int createContextRetriesLeft = 20;
    unsigned long retryDelayMilliseconds = 500;
    bool hadToRetry=false;
    
    context_.reset(); // make sure the pointer is NULL
    
    // the newer NVIDIA drivers sometimes fail to create a context
    // with a "CL_OUT_OF_RESOURCES" error, but work just fine if you
    // try again after a short time.
    for(;;)
    {
        try {
            // create a context
            context_ = shared_ptr<cl::Context>(new cl::Context(devices, properties));
        } catch (cl::Error &err) {
            if ((err.err() == CL_OUT_OF_RESOURCES) && (createContextRetriesLeft>0)) {
                --createContextRetriesLeft;
                
                log_warn("OpenCL ERROR while creating conetxt: CL_OUT_OF_RESOURCES. Trying again in %fs.. (%u tries left)",
                         static_cast<double>(retryDelayMilliseconds)/1000., createContextRetriesLeft);
                
                boost::this_thread::sleep(boost::posix_time::milliseconds(retryDelayMilliseconds));
                
                log_warn("Re-trying to create OpenCL context..");
                hadToRetry=true;
                context_.reset(); // make sure the pointer is NULL
            } else {
                // an OpenCL error here most probably means that there are no devices of the
                // requested type. So just continue quietly.
                log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
                throw I3CLSimStepToPhotonConverter_exception("OpenCL error: could not set up context!");
            }
        }
        
        if (context_) break;
    }
    
    if (hadToRetry)
        log_warn("OpenCL context created successfully!");
    
    {
        log_info("Running on: ");
        std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
        log_info("  * \"%s\"", deviceName.c_str());
        log_info("      ->                      CL_DEVICE_TYPE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_TYPE>()).c_str());
        log_info("      ->         CL_DEVICE_MAX_COMPUTE_UNITS: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()).c_str());
        log_info("      ->    CL_DEVICE_MAX_WORK_ITEM_SIZES[0]: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0]).c_str());
        log_info("      ->       CL_DEVICE_MAX_WORK_GROUP_SIZE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()).c_str());
        log_info("      ->       CL_DEVICE_MAX_CLOCK_FREQUENCY: %sMHz", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>()).c_str());
        log_info("      ->           CL_DEVICE_GLOBAL_MEM_SIZE: %sMiB", boost::lexical_cast<std::string>(static_cast<double>(device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>())/1024./1024.).c_str());
        log_info("      ->  CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE: %sKiB", boost::lexical_cast<std::string>(static_cast<double>(device.getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>())/1024.).c_str());
        log_info("      ->            CL_DEVICE_LOCAL_MEM_TYPE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_LOCAL_MEM_TYPE>()).c_str());
        log_info("      ->            CL_DEVICE_LOCAL_MEM_SIZE: %sKiB", boost::lexical_cast<std::string>(static_cast<double>(device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>())/1024.).c_str());
        log_info("      ->  CL_DEVICE_ERROR_CORRECTION_SUPPORT: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_ERROR_CORRECTION_SUPPORT>()).c_str());
        log_info("      ->             CL_DEVICE_ENDIAN_LITTLE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_ENDIAN_LITTLE>()).c_str());
        log_info("      ->                 CL_DEVICE_AVAILABLE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_AVAILABLE>()).c_str());
        log_info("      ->                    CL_DEVICE_VENDOR: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_VENDOR>()).c_str());
        log_info("      ->                   CL_DEVICE_VERSION: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_VERSION>()).c_str());
        log_info("      ->                CL_DEVICE_EXTENSIONS: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_EXTENSIONS>()).c_str());
    }    
    
    log_debug("Compiling..");
    // accumulate the build options
    std::string BuildOptions;
    
    //BuildOptions += "-Werror ";
    BuildOptions += "-cl-mad-enable ";
    //BuildOptions += "-cl-opt-disable ";
    //BuildOptions += "-cl-no-signed-zeros ";
    //BuildOptions += "-cl-unsafe-math-optimizations ";
    BuildOptions += "-cl-fast-relaxed-math ";
    
    // only valid if extension "cl_nv_compiler_options" is present
    //BuildOptions += "-cl-nv-verbose ";          // Passed on to ptxas as --verbose
    //BuildOptions += "-cl-nv-maxrregcount=60 ";  // Passed on to ptxas as --maxrregcount <N>
    //BuildOptions += "-cl-nv-opt-level=3 ";     // Passed on to ptxas as --opt-level <N>
    
    if (useNativeMath_) {BuildOptions += "-DUSE_NATIVE_MATH ";}
    
    cl::Program program;
    try {
        // build the program
        cl::Program::Sources source;
        
        source.push_back(std::make_pair(mwcrngKernelSource_.c_str(),mwcrngKernelSource_.size()));
        source.push_back(std::make_pair(wlenGeneratorSource_.c_str(),wlenGeneratorSource_.size()));
        source.push_back(std::make_pair(wlenBiasSource_.c_str(),wlenBiasSource_.size()));
        source.push_back(std::make_pair(mediumPropertiesSource_.c_str(),mediumPropertiesSource_.size()));
        source.push_back(std::make_pair(geometrySource_.c_str(),geometrySource_.size()));
        source.push_back(std::make_pair(propagationKernelSource_.c_str(),propagationKernelSource_.size()));
        
        program = cl::Program(*context_, source);
        log_debug("building...");
        program.build(devices, BuildOptions.c_str());
        log_debug("...building finished.");
    } catch (cl::Error &err) {
        log_error("OpenCL ERROR (compile): %s (%i)", err.what(), err.err());
        
        std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
        log_error("  * build status on %s\"", deviceName.c_str());
        log_error("==============================");
        log_error("Build Status: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(device)).c_str());
        log_error("Build Options: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device)).c_str());
        log_error("Build Log: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)).c_str());
        log_error("==============================");
        
        throw I3CLSimStepToPhotonConverter_exception("OpenCL error: could build the OpenCL program!");;
    }
    log_debug("code compiled.");
    
    // instantiate the command queue
    log_debug("Initializing..");
    try {
        for (unsigned int i=0;i<2;++i)
        {
#ifdef DUMP_STATISTICS
            queue_.push_back(shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context_, device, CL_QUEUE_PROFILING_ENABLE)));
#else
            queue_.push_back(shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context_, device, 0)));
#endif
        }
    } catch (cl::Error err) {
        queue_.clear(); // throw away command queue.
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw I3CLSimStepToPhotonConverter_exception("OpenCL error: could not set up command queue!");
    }
    log_debug("initialized.");
    
    // create the kernel
    log_debug("Creating kernel..");
    try {
        // instantiate the kernel object
        for (unsigned int i=0;i<2;++i)
        {
            kernel_.push_back(shared_ptr<cl::Kernel>(new cl::Kernel(program, "propKernel")));
        }

        maxWorkgroupSize_ = kernel_[0]->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
        if (kernel_[0]->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device) != maxWorkgroupSize_) {
            log_fatal("created two identical kernels and got different maximum work group sizes.");
        }
        
        log_debug("Maximum workgroup sizes for the kernel is %" PRIu64, maxWorkgroupSize_);
    } catch (cl::Error err) {
        kernel_.clear(); // throw away command queue.
        queue_.clear(); // throw away command queue.
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw I3CLSimStepToPhotonConverter_exception("OpenCL error: could not create kernel!");
    }
    log_debug("created.");
    
    
    
    
}


void I3CLSimStepToPhotonConverterOpenCL::OpenCLThread()
{
    // do not interrupt this thread by default
    boost::this_thread::disable_interruption di;
    
    try {
        OpenCLThread_impl(di);
    } catch(...) { // any exceptions?
        std::cerr << "OpenCL worker thread died unexpectedly.." << std::endl;
        exit(0); // get out as quickly as possible, we probably just had a FATAL error anyway..
        throw; // will never be reached
    }
}

bool I3CLSimStepToPhotonConverterOpenCL::OpenCLThread_impl_uploadSteps(boost::this_thread::disable_interruption &di,
                                                                       bool &shouldBreak,
                                                                       unsigned int bufferIndex,
                                                                       uint32_t &out_stepsIdentifier,
                                                                       uint64_t &out_totalNumberOfPhotons,
                                                                       std::size_t &out_numberOfInputSteps,
                                                                       bool blocking
                                                                       )
{
    shouldBreak=false;
    
    uint32_t stepsIdentifier;
    I3CLSimStepSeriesConstPtr steps;
    
    const uint32_t zeroCounterBufferSource=0;
    VECTOR_CLASS<cl::Event> bufferWriteEvents(2);

    while (!steps)
    {
        // we need to fetch new steps
        
        boost::this_thread::restore_interruption ri(di);
        try {
            if (blocking) {
                // this can block until there is something on the queue:
                log_trace("[%u] waiting for input queue..", bufferIndex);
                ToOpenCLPair_t val = queueToOpenCL_->Get();
                log_trace("[%u] returned value from input queue..", bufferIndex);
                stepsIdentifier = val.first;
                steps = val.second;
            } else {
                ToOpenCLPair_t val;
                // this will never block:
                log_trace("[%u] waiting for queue.. (non-blocking)", bufferIndex);
                const bool ret = queueToOpenCL_->GetNonBlocking(val);
                
                if (!ret) {
                    log_trace("[%u] returned value from queue (empty), size==%zu/%zu!", bufferIndex, queueToOpenCL_->size(), queueToOpenCL_->max_size());
                    // queue is empty
                    return false;
                }

                log_trace("[%u] returned value from queue (non-empty), size==%zu/%zu!", bufferIndex, queueToOpenCL_->size(), queueToOpenCL_->max_size());

                stepsIdentifier = val.first;
                steps = val.second;
            }
        }
        catch(boost::thread_interrupted &i)
        {
            log_trace("[%u] OpenCL worker thread was interrupted. closing.", bufferIndex);
            shouldBreak=true;
            return true;
        }
    }
    
    log_trace("[%u] OpenCL thread got steps with id %zu", bufferIndex, static_cast<std::size_t>(stepsIdentifier));
    out_stepsIdentifier = stepsIdentifier;
    
#ifdef DUMP_STATISTICS
    uint64_t totalNumberOfPhotons=0;
    BOOST_FOREACH(const I3CLSimStep &step, *steps)
    {
        totalNumberOfPhotons+=step.numPhotons;
    }
    out_totalNumberOfPhotons = totalNumberOfPhotons;
#else
    out_totalNumberOfPhotons = 0;
#endif //DUMP_STATISTICS
    
    log_trace("[%u] copy steps to device", bufferIndex);
    // copy steps to device
    try {
        queue_[bufferIndex]->enqueueWriteBuffer(*deviceBuffer_CurrentNumOutputPhotons[bufferIndex], CL_FALSE, 0, sizeof(uint32_t), &zeroCounterBufferSource, NULL, &(bufferWriteEvents[0]));
        queue_[bufferIndex]->enqueueWriteBuffer(*deviceBuffer_InputSteps[bufferIndex], CL_FALSE, 0, steps->size()*sizeof(I3CLSimStep), &((*steps)[0]), NULL, &(bufferWriteEvents[1]));
        queue_[bufferIndex]->flush(); // make sure it starts executing on the device
        
        log_trace("[%u] waiting for copy to finish", bufferIndex);
        cl::Event::waitForEvents(bufferWriteEvents);
    } catch (cl::Error &err) {
        log_fatal("[%u] OpenCL ERROR (memcpy to device): %s (%i)", bufferIndex, err.what(), err.err());
    }
    log_trace("[%u] copied steps to device", bufferIndex);
    
    out_numberOfInputSteps = steps->size();
    
    return true;
}

void I3CLSimStepToPhotonConverterOpenCL::OpenCLThread_impl_runKernel(unsigned int bufferIndex,
                                                                     cl::Event &kernelFinishEvent,
                                                                     std::size_t numberOfInputSteps)
{
    // run the kernel
    log_trace("[%u] enqueuing kernel..", bufferIndex);

    try {
        // configure which input buffers to use
        queue_[bufferIndex]->enqueueNDRangeKernel(*(kernel_[bufferIndex]), 
                                                  cl::NullRange,	// current implementations force this to be NULL
                                                  cl::NDRange(numberOfInputSteps),	// number of work items
                                                  cl::NDRange(workgroupSize_),
                                                  NULL, //&(bufferWriteEvents),  // wait for buffers to be filled
                                                  &kernelFinishEvent); // signal when finished
        queue_[bufferIndex]->flush(); // make sure it begins executing on the device
    } catch (cl::Error &err) {
        log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
    }

    log_trace("[%u] kernel in queue..", bufferIndex);
}

void I3CLSimStepToPhotonConverterOpenCL::OpenCLThread_impl_downloadPhotons(boost::this_thread::disable_interruption &di,
                                                                           bool &shouldBreak,
                                                                           unsigned int bufferIndex,
                                                                           uint32_t stepsIdentifier)
{
    shouldBreak=false;
   
    I3CLSimPhotonSeriesPtr photons;
    
    try {
        uint32_t numberOfGeneratedPhotons;
        {
            cl::Event copyComplete;
            queue_[bufferIndex]->enqueueReadBuffer(*deviceBuffer_CurrentNumOutputPhotons[bufferIndex], CL_FALSE, 0, sizeof(uint32_t), &numberOfGeneratedPhotons, NULL, &copyComplete);
            queue_[bufferIndex]->flush(); // make sure it starts executing on the device
            copyComplete.wait();
        }
        
        log_trace("Num photons to copy: %" PRIu32, numberOfGeneratedPhotons);
        
        if (numberOfGeneratedPhotons > maxNumOutputPhotons_)
        {
            log_error("Maximum number of photons exceeded, only receiving %" PRIu32 " of %" PRIu32 " photons",
                      maxNumOutputPhotons_, numberOfGeneratedPhotons);
            numberOfGeneratedPhotons = maxNumOutputPhotons_;
        }
        
        if (numberOfGeneratedPhotons>0)
        {
            cl::Event copyComplete;
            
            // allocate the result vector while waiting for the mapping operation to complete
            photons = I3CLSimPhotonSeriesPtr(new I3CLSimPhotonSeries(numberOfGeneratedPhotons));
            
            queue_[bufferIndex]->enqueueReadBuffer(*deviceBuffer_OutputPhotons[bufferIndex], CL_FALSE, 0, numberOfGeneratedPhotons*sizeof(I3CLSimPhoton), &((*photons)[0]), NULL, &copyComplete);
            queue_[bufferIndex]->flush(); // make sure it starts executing on the device
            copyComplete.wait(); // wait for the buffer to be copied
        }
        else
        {
            // empty vector
            photons = I3CLSimPhotonSeriesPtr(new I3CLSimPhotonSeries());
        }
        
    } catch (cl::Error &err) {
        log_fatal("OpenCL ERROR (memcpy from device): %s (%i)", err.what(), err.err());
    }
    
    // we finished simulating.
    // signal the caller by putting it's id on the 
    // output queue.
    {
        boost::this_thread::restore_interruption ri(di);
        try {
            queueFromOpenCL_->Put(ConversionResult_t(stepsIdentifier, photons));
        } catch(boost::thread_interrupted &i) {
            log_debug("OpenCL thread was interrupted. closing.");
            shouldBreak=true;
            return;
        }
    }
    
    
}

#ifdef DUMP_STATISTICS
namespace {
    inline boost::posix_time::ptime DumpStatistics(const cl::Event &kernelFinishEvent,
                                                   const boost::posix_time::ptime &last_timestamp,
                                                   uint64_t totalNumberOfPhotons,
                                                   const std::string &platformName,
                                                   const std::string &deviceName)
    {
        // calculate time since last kernel execution
        boost::posix_time::ptime this_timestamp(boost::posix_time::microsec_clock::universal_time());
        boost::posix_time::time_duration posix_duration = this_timestamp - last_timestamp;
        const double host_duration_in_nanoseconds = static_cast<double>(posix_duration.total_nanoseconds());
        
        uint64_t timeStart, timeEnd;
        kernelFinishEvent.getProfilingInfo(CL_PROFILING_COMMAND_START, &timeStart);
        kernelFinishEvent.getProfilingInfo(CL_PROFILING_COMMAND_END, &timeEnd);
        
        const double kernel_duration_in_nanoseconds = static_cast<double>(timeEnd-timeStart);
        
        const double utilization = kernel_duration_in_nanoseconds/host_duration_in_nanoseconds;
        
#ifdef I3_OPTIMIZE
        std::cout << "kernel statistics: " 
        << kernel_duration_in_nanoseconds/static_cast<double>(totalNumberOfPhotons) 
        << " nanoseconds/photon (util: " << utilization*100. << "%) "
        << "(" << platformName << " " << deviceName << ")"
        << std::endl;
#else
        log_info("kernel statistics: %g nanoseconds/photon (util: %.0f%%) (%s %s)",
                 kernel_duration_in_nanoseconds/static_cast<double>(totalNumberOfPhotons),
                 utilization*100.,
                 platformName.c_str(), deviceName.c_str());
#endif
        
        return this_timestamp;
    }
}
#endif

void I3CLSimStepToPhotonConverterOpenCL::OpenCLThread_impl(boost::this_thread::disable_interruption &di)
{
    // set things up here
    if (!context_) log_fatal("Internal error: context is (null)");

    if (queue_.size() != 2) log_fatal("Internal error: queue_.size() != 2!");
    if (kernel_.size() != 2) log_fatal("Internal error: kernel_.size() != 2!");

    BOOST_FOREACH(shared_ptr<cl::CommandQueue> &ptr, queue_) {
        if (!ptr) log_fatal("Internal error: queue_[] is (null)");
    }
    BOOST_FOREACH(shared_ptr<cl::Kernel> &ptr, kernel_) {
        if (!ptr) log_fatal("Internal error: kernel_[] is (null)");
    }

    if (deviceBuffer_InputSteps.size() != 2) log_fatal("Internal error: deviceBuffer_InputSteps.size() != 2!");
    if (deviceBuffer_OutputPhotons.size() != 2) log_fatal("Internal error: deviceBuffer_OutputPhotons.size() != 2!");
    if (deviceBuffer_CurrentNumOutputPhotons.size() != 2) log_fatal("Internal error: deviceBuffer_CurrentNumOutputPhotons.size() != 2!");
    
    BOOST_FOREACH(shared_ptr<cl::Buffer> &ptr, deviceBuffer_InputSteps) {
        if (!ptr) log_fatal("Internal error: deviceBuffer_InputSteps[] is (null)");
    }
    BOOST_FOREACH(shared_ptr<cl::Buffer> &ptr, deviceBuffer_OutputPhotons) {
        if (!ptr) log_fatal("Internal error: deviceBuffer_OutputPhotons[] is (null)");
    }
    BOOST_FOREACH(shared_ptr<cl::Buffer> &ptr, deviceBuffer_CurrentNumOutputPhotons) {
        if (!ptr) log_fatal("Internal error: deviceBuffer_CurrentNumOutputPhotons[] is (null)");
    }

    if (!deviceBuffer_GeoLayerToOMNumIndexPerStringSet) log_fatal("Internal error: deviceBuffer_GeoLayerToOMNumIndexPerStringSet is (null)");
    if (!deviceBuffer_MWC_RNG_x) log_fatal("Internal error: deviceBuffer_MWC_RNG_x is (null)");
    if (!deviceBuffer_MWC_RNG_a) log_fatal("Internal error: deviceBuffer_MWC_RNG_a is (null)");
    
    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(openCLStarted_mutex_);
        openCLStarted_=true;
    }
    openCLStarted_cond_.notify_all();
    
    std::vector<uint32_t> stepsIdentifier(2, 0);
    std::vector<uint64_t> totalNumberOfPhotons(2, 0);
    std::vector<std::size_t> numberOfSteps(2, 0);
    
#ifdef DUMP_STATISTICS
    boost::posix_time::ptime last_timestamp(boost::posix_time::microsec_clock::universal_time());
#endif
    
    unsigned int thisBuffer=1;
    unsigned int otherBuffer=0;
    bool otherBufferHasBeenCopied=false;
    
    // start the main loop
    for (;;)
    {
        // swap buffers
        std::swap(thisBuffer, otherBuffer);
        
        log_trace("buffers indices now: this==%u, other==%u", thisBuffer, otherBuffer);
        
        if (!otherBufferHasBeenCopied) {
            log_trace("[%u] starting \"this\" buffer copy (need to block)..", thisBuffer);
            {
                bool shouldBreak=false; // shouldBreak is true if this thread has been signalled to terminate
                OpenCLThread_impl_uploadSteps(di, shouldBreak, thisBuffer, stepsIdentifier[thisBuffer], totalNumberOfPhotons[thisBuffer], numberOfSteps[thisBuffer]);
                if (shouldBreak) break; // is thread termination being requested?
            }
            log_trace("[%u] this buffer has been copied..", thisBuffer);
        } else {
            // else: this buffer is already there!
            log_trace("[%u] buffer is already there!", thisBuffer);
        }
        
        // reset the "has-been-copied" flag
        otherBufferHasBeenCopied=false;
        
        
        // start the kernel
        cl::Event kernelFinishEvent;
        OpenCLThread_impl_runKernel(thisBuffer, kernelFinishEvent, numberOfSteps[thisBuffer]);

        // if there already is a new buffer available, copy it now, while the kernel is running
        log_trace("[%u] Starting copy (other buffer)..", otherBuffer);
        {
            bool shouldBreak=false; // shouldBreak is true if this thread has been signalled to terminate
            bool gotSomething = OpenCLThread_impl_uploadSteps(di, shouldBreak, otherBuffer, stepsIdentifier[otherBuffer], totalNumberOfPhotons[otherBuffer], numberOfSteps[otherBuffer], false);
            if (shouldBreak) break;
            
            if (!gotSomething) {
                log_trace("[%u] copy (other buffer): queue empty!", otherBuffer);
                // nothing on the queue
                otherBufferHasBeenCopied=false;
            } else {
                log_trace("[%u] copy (other buffer):  done!", otherBuffer);
                otherBufferHasBeenCopied=true;
            }
        }
        
        
        log_trace("[%u] waiting for kernel..", thisBuffer);

        try {
            // wait for the kernel to finish
            kernelFinishEvent.wait();
        } catch (cl::Error &err) {
            log_fatal("[%u] OpenCL ERROR (running kernel): %s (%i)", thisBuffer, err.what(), err.err());
        }

        log_trace("[%u] kernel finished..", thisBuffer);

#ifdef DUMP_STATISTICS
        log_trace("[%u] dumping statistics..", thisBuffer);

        last_timestamp = DumpStatistics(kernelFinishEvent,
                                        last_timestamp,
                                        totalNumberOfPhotons[thisBuffer],
                                        device_->GetPlatformName(),
                                        device_->GetDeviceName());
#endif

        log_trace("[%u] waiting for queue..", thisBuffer);

        try {
            // wait for the queue to really finish (just to make sure)
            queue_[thisBuffer]->finish();
        } catch (cl::Error &err) {
            log_fatal("[%u] OpenCL ERROR (running kernel): %s (%i)", thisBuffer, err.what(), err.err());
        }
        
        log_trace("[%u] queue finished!", thisBuffer);
        
        // receive results
        log_trace("[%u] receiving results..!", thisBuffer);
        {
            bool shouldBreak;
            OpenCLThread_impl_downloadPhotons(di, shouldBreak, thisBuffer, stepsIdentifier[thisBuffer]);
            if (shouldBreak) break; // is thread termination being requested?
        }
        log_trace("[%u] results received.", thisBuffer);

    }
    
    log_debug("OpenCL thread terminating...");
    
    // shut down
    
    log_debug("OpenCL thread terminated.");
}

bool I3CLSimStepToPhotonConverterOpenCL::IsInitialized() const
{
    return initialized_;
}

void I3CLSimStepToPhotonConverterOpenCL::SetWlenGenerator(I3CLSimRandomValueConstPtr wlenGenerator)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    wlenGenerator_=wlenGenerator;
}

void I3CLSimStepToPhotonConverterOpenCL::SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    wlenBias_=wlenBias;
}

void I3CLSimStepToPhotonConverterOpenCL::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    mediumProperties_=mediumProperties;
}

void I3CLSimStepToPhotonConverterOpenCL::SetGeometry(I3CLSimSimpleGeometryConstPtr geometry)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    geometry_=geometry;
}

void I3CLSimStepToPhotonConverterOpenCL::EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier)
{
    if (!initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL is not initialized!");
    
    if (!steps)
        throw I3CLSimStepToPhotonConverter_exception("Steps pointer is (null)!");
    
    if (steps->empty())
        throw I3CLSimStepToPhotonConverter_exception("Steps are empty!");
    
    if (steps->size() > maxNumWorkitems_)
        throw I3CLSimStepToPhotonConverter_exception("Number of steps is greater than maximum number of work items!");
    
    if (steps->size() % workgroupSize_ != 0)
        throw I3CLSimStepToPhotonConverter_exception("The number of steps is not a multiple of the workgroup size!");
    
    
    queueToOpenCL_->Put(make_pair(identifier, steps));
}

std::size_t I3CLSimStepToPhotonConverterOpenCL::QueueSize() const
{
    if (!initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL is not initialized!");
    
    return queueToOpenCL_->size();
}


bool I3CLSimStepToPhotonConverterOpenCL::MorePhotonsAvailable() const
{
    if (!initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL is not initialized!");
    
    return (!queueFromOpenCL_->empty());
}

// helper
namespace {
    inline void ReplaceStringDOMIndexWithStringDOMIDs(I3CLSimPhotonSeries &photons,
                                                      const std::vector<int> &stringIndexToStringIDBuffer,
                                                      const std::vector<std::vector<unsigned int> > &domIndexToDomIDBuffer_perStringIndex)
    {
        BOOST_FOREACH(I3CLSimPhoton &photon, photons)
        {
            const int16_t stringIndex = photon.stringID;
            const uint16_t DOMIndex = photon.omID;
            
            const int stringID = stringIndexToStringIDBuffer.at(stringIndex);
            const unsigned int domID = domIndexToDomIDBuffer_perStringIndex.at(stringIndex).at(DOMIndex);
            
            if ((stringID < std::numeric_limits<int16_t>::min()) ||
                (stringID > std::numeric_limits<int16_t>::max()))
                log_fatal("Your detector I3Geometry uses a string ID \"%i\". Large IDs like that are currently not supported by clsim.",
                          stringID);
            
            if (domID > std::numeric_limits<uint16_t>::max())
                log_fatal("Your detector I3Geometry uses a OM ID \"%u\". Large IDs like that are currently not supported by clsim.",
                          domID);
            
            photon.stringID = static_cast<int16_t>(stringID);
            photon.omID = static_cast<uint16_t>(domID);
            
            log_trace("Replaced ID (%" PRIi16 "/%" PRIu16 ") with ID (%" PRIi16 "/%" PRIu16 ") (photon @ pos=(%g,%g,%g))",
                      stringIndex,
                      DOMIndex,
                      photon.stringID,
                      photon.omID,
                      photon.GetPosX(),
                      photon.GetPosY(),
                      photon.GetPosZ()
                      );
            
        }
    }
    
}

I3CLSimStepToPhotonConverter::ConversionResult_t I3CLSimStepToPhotonConverterOpenCL::GetConversionResult()
{
    if (!initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL is not initialized!");
    
    ConversionResult_t result = queueFromOpenCL_->Get();
    
    if (result.second) {
        ReplaceStringDOMIndexWithStringDOMIDs(*result.second,
                                              stringIndexToStringIDBuffer_,
                                              domIndexToDomIDBuffer_perStringIndex_);
        
    }
    
    return result;
}
