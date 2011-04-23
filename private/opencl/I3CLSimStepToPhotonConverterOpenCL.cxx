#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"

// debugging: show GPUtime/photon
#define DUMP_STATISTICS

#define __STDC_FORMAT_MACROS 
#include <inttypes.h>

#include <string>
#include <sstream>
#include <algorithm>
#include <limits>

#include <stdlib.h>
#include <boost/foreach.hpp>

#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"
#include "opencl/I3CLSimHelperGenerateGeometrySource.h"

#include "opencl/mwcrng_init.h"

#define __CL_ENABLE_EXCEPTIONS
#include "clsim/cl.hpp"

const bool I3CLSimStepToPhotonConverterOpenCL::default_useNativeMath=true;
const bool I3CLSimStepToPhotonConverterOpenCL::default_cpuOnly=false;
const bool I3CLSimStepToPhotonConverterOpenCL::default_gpuOnly=false;


I3CLSimStepToPhotonConverterOpenCL::I3CLSimStepToPhotonConverterOpenCL(I3RandomServicePtr randomService,
                                                                       bool useNativeMath,
                                                                       bool cpuOnly,
                                                                       bool gpuOnly)
:
openCLStarted_(false),
queueToOpenCL_(new I3CLSimQueue<ToOpenCLPair_t>(5)),
queueFromOpenCL_(new I3CLSimQueue<I3CLSimStepToPhotonConverter::ConversionResult_t>(0)),
randomService_(randomService),
initialized_(false),
compiled_(false),
cpuOnly_(cpuOnly),
gpuOnly_(gpuOnly),
useNativeMath_(useNativeMath),
selectedDeviceIndex_(0),
deviceIndexIsSelected_(false),
maxWorkgroupSize_(0),
workgroupSize_(0),
maxNumWorkitems_(10240)
{
    if (!randomService_) log_fatal("You need to supply a I3RandomService.");
    
    // enumerate platforms and devices
    shared_ptr<std::vector<std::pair<std::string, std::string> > > deviceNameList(new std::vector<std::pair<std::string, std::string> >());
    
    std::vector<cl::Platform> platforms;
    
    try {
		// get a list of available platforms
		cl::Platform::get(&platforms);
	} catch (cl::Error &err) {
		log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw I3CLSimStepToPhotonConverter_exception("OpenCL error.");
	}
    
    BOOST_FOREACH(cl::Platform &platform, platforms)
    {
        const std::string platformName = platform.getInfo<CL_PLATFORM_NAME>();
        
        std::vector<cl::Device> devices;
        
        try {
            if ((!cpuOnly) && (!gpuOnly))
                platform.getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);
            else if (cpuOnly)
                platform.getDevices(CL_DEVICE_TYPE_CPU, &devices);
            else if (gpuOnly)
                platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
            else
                throw I3CLSimStepToPhotonConverter_exception("Cannot set both CPUOnly and GPUOnly!.");
        } catch (cl::Error &err) {
            // an OpenCL error here most probably means that there are no devices of the
            // requested type. So just continue quietly.
            log_debug("OpenCL ERROR: %s (%i)", err.what(), err.err());
            continue;
        }
        
        BOOST_FOREACH(cl::Device &device, devices)
        {
            const std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
            
            deviceNameList->push_back(std::make_pair(platformName, deviceName));
            clPlatformDeviceList_.push_back(std::make_pair(shared_ptr<cl::Platform>(new cl::Platform(platform)), shared_ptr<cl::Device>(new cl::Device(device))));
            
            log_trace("PLATFORM: \"%s\" -> DEVICE: \"%s\"",
                      platformName.c_str(),
                      deviceName.c_str());
        }
    }
    
    deviceNameList_=deviceNameList;
    
    
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
    deviceBuffer_InputSteps.reset();
    deviceBuffer_OutputPhotons.reset();
    deviceBuffer_CurrentNumOutputPhotons.reset();
    deviceBuffer_GeoLayerToOMNumIndexPerStringSet.reset();
    
    // reset pointers
    compiled_=false;
    context_.reset();
    kernel_.reset();
    queue_.reset();
    
}

shared_ptr<const std::vector<std::pair<std::string, std::string> > >
I3CLSimStepToPhotonConverterOpenCL::GetDeviceList() const
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");

    return deviceNameList_;
}

uint64_t I3CLSimStepToPhotonConverterOpenCL::GetMaxWorkgroupSize() const
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");

    if (!compiled_)
        throw I3CLSimStepToPhotonConverter_exception("You need to compile the kernel first. Call Compile().");

    return maxWorkgroupSize_;
}

void I3CLSimStepToPhotonConverterOpenCL::SetDeviceIndex(std::size_t selectedDeviceIndex)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    if (selectedDeviceIndex_ >= clPlatformDeviceList_.size())
        throw I3CLSimStepToPhotonConverter_exception("Invalid device index!");
    
    if (selectedDeviceIndex!=selectedDeviceIndex_)
    {
        compiled_=false;
        kernel_.reset();
        queue_.reset();
    }
    
    selectedDeviceIndex_ = selectedDeviceIndex;
    deviceIndexIsSelected_=true;
}

void I3CLSimStepToPhotonConverterOpenCL::SetDeviceName(const std::string &platform, const std::string &device)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");

    for (std::size_t i=0; i<deviceNameList_->size(); ++i)
    {
        if (((*deviceNameList_)[i].first == platform) &&
            ((*deviceNameList_)[i].second == device))
            return SetDeviceIndex(i);
    }
    
    throw I3CLSimStepToPhotonConverter_exception("Device does not exist!");
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
    log_debug("Setting up RNG for %" PRIu64 " workitems.", maxNumWorkitems_);

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
    deviceBuffer_InputSteps.reset();
    deviceBuffer_OutputPhotons.reset();
    deviceBuffer_CurrentNumOutputPhotons.reset();
    deviceBuffer_GeoLayerToOMNumIndexPerStringSet.reset();
    
    
    // set up device buffers from existing host buffers
    deviceBuffer_MWC_RNG_x = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_x.size() * sizeof(uint64_t), &(MWC_RNG_x[0])));

    deviceBuffer_MWC_RNG_a = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_a.size() * sizeof(uint32_t), &(MWC_RNG_a[0])));

    deviceBuffer_GeoLayerToOMNumIndexPerStringSet = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, geoLayerToOMNumIndexPerStringSetInfo_.size() * sizeof(unsigned char), &(geoLayerToOMNumIndexPerStringSetInfo_[0])));

    
    // allocate empty buffers on the device
    deviceBuffer_InputSteps = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, maxNumWorkitems_*sizeof(I3CLSimStep), NULL));
    
    deviceBuffer_OutputPhotons = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, maxNumOutputPhotons_*sizeof(I3CLSimPhoton), NULL));

    deviceBuffer_CurrentNumOutputPhotons = shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(uint32_t), NULL));

    log_debug("Device buffers are set up.");

    log_debug("Configuring kernel.");
    {
        kernel_->setArg(0, *deviceBuffer_CurrentNumOutputPhotons);          // hit counter
        kernel_->setArg(1, maxNumOutputPhotons_);                           // maximum number of possible hits
        
        kernel_->setArg(2, *deviceBuffer_GeoLayerToOMNumIndexPerStringSet); // additional geometry information (did not fit into constant memory)
        
        kernel_->setArg(3, *deviceBuffer_InputSteps);                       // the input steps
        kernel_->setArg(4, *deviceBuffer_OutputPhotons);                    // the output photons
        
        kernel_->setArg(5, *deviceBuffer_MWC_RNG_x);                        // rng state
        kernel_->setArg(6, *deviceBuffer_MWC_RNG_a);                        // rng state
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

    initialized_=true;
}


void I3CLSimStepToPhotonConverterOpenCL::Compile()
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");

    if (compiled_) return; // silently
    
    if (clPlatformDeviceList_.empty())
        throw I3CLSimStepToPhotonConverter_exception("No OpenCL devices found!");

    if (!wlenGenerator_)
        throw I3CLSimStepToPhotonConverter_exception("WlenGenerator not set!");

    if (!wlenBias_)
        throw I3CLSimStepToPhotonConverter_exception("WlenBias not set!");

    if (!mediumProperties_)
        throw I3CLSimStepToPhotonConverter_exception("MediumProperties not set!");
    
    if (!geometry_)
        throw I3CLSimStepToPhotonConverter_exception("Geometry not set!");
    
    if (!deviceIndexIsSelected_)
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
    
    SetupQueueAndKernel(*(clPlatformDeviceList_[selectedDeviceIndex_].first),
                        *(clPlatformDeviceList_[selectedDeviceIndex_].second));
    
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
    std::vector<cl::Device> devices(1, device);

    try {
        // prepare a device vector (containing a single device)
        cl_context_properties properties[] = 
        { CL_CONTEXT_PLATFORM, (cl_context_properties)(platform)(), 0};
        
        // create a context
        context_ = shared_ptr<cl::Context>(new cl::Context(devices, properties));
    } catch (cl::Error &err) {
        // an OpenCL error here most probably means that there are no devices of the
        // requested type. So just continue quietly.
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw I3CLSimStepToPhotonConverter_exception("OpenCL error: could not set up context!");
    }
    
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
#ifdef DUMP_STATISTICS
		queue_ = shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context_, device, CL_QUEUE_PROFILING_ENABLE));
#else
		queue_ = shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context_, device, 0));
#endif
	} catch (cl::Error err) {
        queue_.reset(); // throw away command queue.
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw I3CLSimStepToPhotonConverter_exception("OpenCL error: could not set up command queue!");
	}
    log_debug("initialized.");
    
    // create the kernel
	log_debug("Creating kernel..");
    try {
		// instantiate the kernel object
		kernel_ = shared_ptr<cl::Kernel>(new cl::Kernel(program, "propKernel"));
        maxWorkgroupSize_ = kernel_->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
        
		log_debug("Maximum workgroup sizes for the kernel is %" PRIu64, maxWorkgroupSize_);
	} catch (cl::Error err) {
        kernel_.reset(); // throw away command queue.
        queue_.reset(); // throw away command queue.
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
        log_fatal("OpenCL worker thread died unexpectedly..");
        throw;
    }
}

void I3CLSimStepToPhotonConverterOpenCL::OpenCLThread_impl(boost::this_thread::disable_interruption &di)
{
    // set things up here
    if (!context_) log_fatal("Internal error: context is (null)");
    if (!queue_) log_fatal("Internal error: queue is (null)");
    if (!kernel_) log_fatal("Internal error: kernel is (null)");

    if (!deviceBuffer_InputSteps) log_fatal("Internal error: deviceBuffer_InputSteps is (null)");
    if (!deviceBuffer_OutputPhotons) log_fatal("Internal error: deviceBuffer_OutputPhotons is (null)");
    if (!deviceBuffer_CurrentNumOutputPhotons) log_fatal("Internal error: deviceBuffer_CurrentNumOutputPhotons is (null)");
    if (!deviceBuffer_GeoLayerToOMNumIndexPerStringSet) log_fatal("Internal error: deviceBuffer_GeoLayerToOMNumIndexPerStringSet is (null)");
    if (!deviceBuffer_MWC_RNG_x) log_fatal("Internal error: deviceBuffer_MWC_RNG_x is (null)");
    if (!deviceBuffer_MWC_RNG_a) log_fatal("Internal error: deviceBuffer_MWC_RNG_a is (null)");
    
    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(openCLStarted_mutex_);
        openCLStarted_=true;
    }
    openCLStarted_cond_.notify_all();


    I3CLSimStepSeriesConstPtr steps;
    uint32_t stepsIdentifier=0;

    // start the main loop
    for (;;)
    {
        
        if (!steps)
        {
            // we need to fetch new steps
            
            boost::this_thread::restore_interruption ri(di);
            try {
                ToOpenCLPair_t val = queueToOpenCL_->Get();
                stepsIdentifier = val.first;
                steps = val.second;
            }
            catch(boost::thread_interrupted &i)
            {
                log_debug("OpenCL worker thread was interrupted. closing.");
                break;
            }
        }
        
        if (!steps) {
            log_warn("OpenCL thread got NULL! ignoring.");
            continue;
        }

        log_trace("OpenCL thread got steps with id %zu", static_cast<std::size_t>(stepsIdentifier));
        
#ifdef DUMP_STATISTICS
        uint64_t totalNumberOfPhotons=0;
        BOOST_FOREACH(const I3CLSimStep &step, *steps)
        {
            totalNumberOfPhotons+=step.numPhotons;
        }
        
#endif //DUMP_STATISTICS
        
        // copy steps to device
        try {
            // blocking!
            uint32_t *mapped_CurrentNumOutputPhotons = (uint32_t *)queue_->enqueueMapBuffer(*deviceBuffer_CurrentNumOutputPhotons, CL_TRUE, CL_MAP_WRITE, 0, sizeof(uint32_t));
            I3CLSimStep *mapped_InputSteps = (I3CLSimStep *)queue_->enqueueMapBuffer(*deviceBuffer_InputSteps, CL_TRUE, CL_MAP_WRITE, 0, steps->size()*sizeof(I3CLSimStep));
            //queue_->finish();

            *mapped_CurrentNumOutputPhotons=0; // reset the photon count (== result count)
            memcpy(mapped_InputSteps, &((*steps)[0]), steps->size()*sizeof(I3CLSimStep));
            
            queue_->enqueueUnmapMemObject(*deviceBuffer_InputSteps, mapped_InputSteps);
            queue_->enqueueUnmapMemObject(*deviceBuffer_CurrentNumOutputPhotons, mapped_CurrentNumOutputPhotons);
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (memcpy to device): %s (%i)", err.what(), err.err());
        }

        const std::size_t thisNumberOfInputSteps = steps->size();
        
        log_trace("Starting kernel..");
        
        cl::Event kernelFinishEvent;
        
        // run the kernel
        try {
            queue_->enqueueNDRangeKernel(*kernel_, 
                                         cl::NullRange,	// current implementations force this to be NULL
                                         cl::NDRange(thisNumberOfInputSteps),	// number of work items
                                         cl::NDRange(workgroupSize_),
                                         NULL,
                                         &kernelFinishEvent);
            queue_->flush(); // make sure it begins executing on the device
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
        }

        // we don't need the current steps anymore
        steps.reset();

        try {
            // wait for the kernel to finish
            kernelFinishEvent.wait();

#ifdef DUMP_STATISTICS
            uint64_t timeStart, timeEnd;
            kernelFinishEvent.getProfilingInfo(CL_PROFILING_COMMAND_START, &timeStart);
            kernelFinishEvent.getProfilingInfo(CL_PROFILING_COMMAND_END, &timeEnd);

            log_warn("kernel statistics: %g nanoseconds/photon",
                     static_cast<double>(timeEnd-timeStart)/static_cast<double>(totalNumberOfPhotons));
#endif

            // wait for the queue to really finish (just to make sure)
            queue_->finish();
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
        }
        
        log_debug("Kernel finished!");

        // receive results
        I3CLSimPhotonSeriesPtr photons;

        try {
            uint32_t numberOfGeneratedPhotons;
            {
                cl::Event mappingComplete;
                uint32_t *mapped_CurrentNumOutputPhotons = (uint32_t *)queue_->enqueueMapBuffer(*deviceBuffer_CurrentNumOutputPhotons, CL_FALSE, CL_MAP_READ, 0, sizeof(uint32_t), NULL, &mappingComplete);
                queue_->flush(); // make sure it begins executing on the device
                mappingComplete.wait();
                numberOfGeneratedPhotons = *mapped_CurrentNumOutputPhotons;
                queue_->enqueueUnmapMemObject(*deviceBuffer_CurrentNumOutputPhotons, mapped_CurrentNumOutputPhotons);
            }
            
            log_trace("Num photons to copy: %" PRIu32, numberOfGeneratedPhotons);
            
            if (numberOfGeneratedPhotons > maxNumOutputPhotons_)
            {
                log_warn("Maximum number of photons exceeded, only receiving %" PRIu32 " of %" PRIu32 " photons",
                         maxNumOutputPhotons_, numberOfGeneratedPhotons);
                numberOfGeneratedPhotons = maxNumOutputPhotons_;
            }
            
            if (numberOfGeneratedPhotons>0)
            {
                cl::Event mappingComplete;
                I3CLSimPhoton *mapped_OutputPhotons = (I3CLSimPhoton *)queue_->enqueueMapBuffer(*deviceBuffer_OutputPhotons, CL_FALSE, CL_MAP_READ, 0, numberOfGeneratedPhotons*sizeof(I3CLSimPhoton), NULL, &mappingComplete);
                queue_->flush(); // make sure it begins executing on the device
                
                // allocate the result vector while waiting for the mapping operation to complete
                photons = I3CLSimPhotonSeriesPtr(new I3CLSimPhotonSeries(numberOfGeneratedPhotons));
                
                // wait for the buffer to be mapped
                mappingComplete.wait();
                
                // copy it to the host
                memcpy(&((*photons)[0]), mapped_OutputPhotons, numberOfGeneratedPhotons*sizeof(I3CLSimPhoton));
                
                queue_->enqueueUnmapMemObject(*deviceBuffer_OutputPhotons, mapped_OutputPhotons);
            }
            else
            {
                // empty vector
                photons = I3CLSimPhotonSeriesPtr(new I3CLSimPhotonSeries());
            }
            
            
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (memcpy to device): %s (%i)", err.what(), err.err());
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
                break;
            }
        }
        
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
    kernel_.reset();
    queue_.reset();
    
    wlenGenerator_=wlenGenerator;
}

void I3CLSimStepToPhotonConverterOpenCL::SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.reset();
    queue_.reset();
    
    wlenBias_=wlenBias;
}

void I3CLSimStepToPhotonConverterOpenCL::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");

    compiled_=false;
    kernel_.reset();
    queue_.reset();
    
    mediumProperties_=mediumProperties;
}

void I3CLSimStepToPhotonConverterOpenCL::SetGeometry(I3CLSimSimpleGeometryConstPtr geometry)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.reset();
    queue_.reset();
    
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
            const uint32_t stringIndexWithDOMIndex = *reinterpret_cast<cl_uint *>(&photon.dummy);
            
            const uint32_t stringIndex = stringIndexWithDOMIndex/1000;
            const uint32_t DOMIndex = stringIndexWithDOMIndex%1000;
            
            const int stringID = stringIndexToStringIDBuffer.at(stringIndex);

            const unsigned int domID = domIndexToDomIDBuffer_perStringIndex.at(stringIndex).at(DOMIndex);

            if (stringID >= 0)
                photon.dummy = stringID*1000 + domID;
            else
                photon.dummy = -((-stringID)*1000 + domID);
            
            log_trace("Replaced dummy %u with %i (photon @ pos=(%g,%g,%g))",
                      stringIndexWithDOMIndex,
                      photon.dummy,
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
