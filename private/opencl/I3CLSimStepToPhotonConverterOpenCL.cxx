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
 * @file I3CLSimStepToPhotonConverterOpenCL.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <cmath>

#include <boost/date_time/posix_time/posix_time.hpp>
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

#include <icetray/I3Units.h>

#include "clsim/I3CLSimHelperToFloatString.h"
#include "opencl/I3CLSimHelperMath.h"

#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"
#include "opencl/I3CLSimHelperGenerateGeometrySource.h"

#include "opencl/mwcrng_init.h"

#define __CL_ENABLE_EXCEPTIONS
#include "clsim/cl.hpp"

using namespace I3CLSimHelper;

const bool I3CLSimStepToPhotonConverterOpenCL::default_useNativeMath=true;


I3CLSimStepToPhotonConverterOpenCL::I3CLSimStepToPhotonConverterOpenCL(I3RandomServicePtr randomService,
                                                                       bool useNativeMath)
:
statistics_total_device_duration_in_nanoseconds_(0),
statistics_total_host_duration_in_nanoseconds_(0),
statistics_total_kernel_calls_(0),
statistics_total_num_photons_generated_(0),
statistics_total_num_photons_atDOMs_(0),
openCLStarted_(false),
queueToOpenCL_(new I3CLSimQueue<ToOpenCLPair_t>(5)),
queueFromOpenCL_(new I3CLSimQueue<I3CLSimStepToPhotonConverter::ConversionResult_t>(0)),
randomService_(randomService),
initialized_(false),
compiled_(false),
useNativeMath_(useNativeMath),
deviceIsSelected_(false),
disableDoubleBuffering_(true),
doublePrecision_(false),
stopDetectedPhotons_(false),
saveAllPhotons_(false),
saveAllPhotonsPrescale_(0.001), // only save .1% of all photons when in "AllPhotons" mode
fixedNumberOfAbsorptionLengths_(NAN),
pancakeFactor_(1.),
photonHistoryEntries_(0),
maxWorkgroupSize_(0),
workgroupSize_(0),
maxNumWorkitems_(10240)
{
    if (!randomService_) log_fatal("You need to supply a I3RandomService.");
    
    // load program source from files
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels";
    
    try {
        mwcrngKernelSource_ = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/mwcrng_kernel.cl");
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
    deviceBuffer_PhotonHistory.clear();

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
    
    log_debug("Setting up OpenCL..");
    
    Compile();
    
    // use the maximum workgroup size (==granularity) if none has been configured
    if (workgroupSize_==0) workgroupSize_=maxWorkgroupSize_;
    
    if (workgroupSize_>maxWorkgroupSize_)
        throw I3CLSimStepToPhotonConverter_exception("Workgroup size too large!");
    
    const unsigned int numBuffers = disableDoubleBuffering_?1:2;
    
    // For GPUs, choose the largest number of work items such that they still
    // fit in device memory
    if (device_->IsGPU()) {
      size_t sizePerWorkitem = numBuffers*(
        sizeof(I3CLSimStep)  // the input step
        + 2*sizeof(uint64_t) // MWC multipliers
        + std::max(10u,
          unsigned(10000.*(saveAllPhotons_ ? saveAllPhotonsPrescale_ : 0)))
          *sizeof(I3CLSimPhoton) // the output buffer
        )
        + kernel_[0]->getWorkGroupInfo<CL_KERNEL_LOCAL_MEM_SIZE>(
          *device_->GetDeviceHandle()) // the kernel itself
        ;
      maxNumWorkitems_ = (device_->GetGlobalMemSize()
        - geoLayerToOMNumIndexPerStringSetInfo_.size()*sizeof(unsigned short)
        - numBuffers*sizeof(uint32_t))/sizePerWorkitem;
      size_t numMultipliers = 6139850;
      if (maxNumWorkitems_ > numMultipliers) {
        log_info_stream("Limiting number of work items to "<<numMultipliers<<" (maximum number of prime multipliers)");
        maxNumWorkitems_ = numMultipliers;
      }
      // Choose a bunch size that is a multiple of both the number of cores
      // and the workgroup size
      size_t granularity = device_->GetMaxComputeUnits()*workgroupSize_;
      maxNumWorkitems_ = (maxNumWorkitems_/granularity)*granularity;
    }
    
    if (maxNumWorkitems_%workgroupSize_ != 0)
        throw I3CLSimStepToPhotonConverter_exception("The maximum number of work items (" + boost::lexical_cast<std::string>(maxNumWorkitems_) + ") must be a multiple of the workgroup size (" + boost::lexical_cast<std::string>(workgroupSize_) + ").");
    
    log_debug("basic OpenCL setup done.");
    
    if (!saveAllPhotons_) {
        // start with a maximum number of output photons of the same size as the number of
        // input steps. Should be plenty..
        maxNumOutputPhotons_ = static_cast<uint32_t>(std::min(maxNumWorkitems_*10, static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())));
        if (maxNumOutputPhotons_ < 1000) maxNumOutputPhotons_=1000; // use a sane minimum output buffer size
    } else {
        // we need a lot more space for photon storage in case all photons are to be saved
        std::size_t sizeIncreaseFactor = 10000.*saveAllPhotonsPrescale_;
        if (sizeIncreaseFactor < 1) sizeIncreaseFactor=1;
        
        maxNumOutputPhotons_ = static_cast<uint32_t>(std::min(maxNumWorkitems_*sizeIncreaseFactor, static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())));
    }
    
    // set up rng
    log_debug("Setting up RNG for %zu workitems.", maxNumWorkitems_);
    
    MWC_RNG_x.resize(maxNumWorkitems_);
    MWC_RNG_a.resize(maxNumWorkitems_);
    
    if (init_MWC_RNG(&(MWC_RNG_x[0]), &(MWC_RNG_a[0]), maxNumWorkitems_, randomService_)!=0) 
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    log_debug("RNG is set up..");
    
    log_debug("Setting up device buffers..");
    
    // reset all buffers first
    deviceBuffer_MWC_RNG_x.reset();
    deviceBuffer_MWC_RNG_a.reset();
    deviceBuffer_InputSteps.clear();
    deviceBuffer_OutputPhotons.clear();
    deviceBuffer_PhotonHistory.clear();
    deviceBuffer_CurrentNumOutputPhotons.clear();
    deviceBuffer_GeoLayerToOMNumIndexPerStringSet.reset();
    
    
    // set up device buffers from existing host buffers
    deviceBuffer_MWC_RNG_x = boost::shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_x.size() * sizeof(uint64_t), &(MWC_RNG_x[0])));
    
    deviceBuffer_MWC_RNG_a = boost::shared_ptr<cl::Buffer>
    (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_a.size() * sizeof(uint32_t), &(MWC_RNG_a[0])));
    
    if (!saveAllPhotons_) {
        // no need for a geometry buffer if all photons are saved and no
        // geometry is necessary.
        deviceBuffer_GeoLayerToOMNumIndexPerStringSet = boost::shared_ptr<cl::Buffer>
        (new cl::Buffer(*context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, geoLayerToOMNumIndexPerStringSetInfo_.size() * sizeof(unsigned short), &(geoLayerToOMNumIndexPerStringSetInfo_[0])));
    }
    
    // allocate empty buffers on the device
    for (unsigned int i=0;i<numBuffers;++i)
    {
        deviceBuffer_InputSteps.push_back(boost::shared_ptr<cl::Buffer>
        (new cl::Buffer(*context_, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, maxNumWorkitems_*sizeof(I3CLSimStep), NULL)));
        
        deviceBuffer_OutputPhotons.push_back(boost::shared_ptr<cl::Buffer>
        (new cl::Buffer(*context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, maxNumOutputPhotons_*sizeof(I3CLSimPhoton), NULL)));
        
        deviceBuffer_CurrentNumOutputPhotons.push_back(boost::shared_ptr<cl::Buffer>
        (new cl::Buffer(*context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(uint32_t), NULL)));

        if (photonHistoryEntries_>0) {
            deviceBuffer_PhotonHistory.push_back
            (boost::shared_ptr<cl::Buffer>
             (new cl::Buffer(*context_,
                             CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR,
                             static_cast<std::size_t>(maxNumOutputPhotons_)*static_cast<std::size_t>(photonHistoryEntries_)*sizeof(cl_float4),
                             NULL
                            )
             )
            );
        }
    }
    
    log_debug("Device buffers are set up.");
    
    log_debug("Configuring kernel.");
    for (unsigned int i=0;i<numBuffers;++i)
    {
        unsigned argN=0;
        
        kernel_[i]->setArg(argN++, *(deviceBuffer_CurrentNumOutputPhotons[i]));     // hit counter
        kernel_[i]->setArg(argN++, maxNumOutputPhotons_);                           // maximum number of possible hits
        
        if (!saveAllPhotons_) {
            kernel_[i]->setArg(argN++, *deviceBuffer_GeoLayerToOMNumIndexPerStringSet); // additional geometry information (did not fit into constant memory)
        }
        
        kernel_[i]->setArg(argN++, *(deviceBuffer_InputSteps[i]));                  // the input steps
        kernel_[i]->setArg(argN++, *(deviceBuffer_OutputPhotons[i]));               // the output photons

        if (photonHistoryEntries_>0) {
            kernel_[i]->setArg(argN++, *(deviceBuffer_PhotonHistory[i]));           // the photon history (the last N points where the photon scattered)
        }

        kernel_[i]->setArg(argN++, *deviceBuffer_MWC_RNG_x);                    // rng state
        kernel_[i]->setArg(argN++, *deviceBuffer_MWC_RNG_a);                    // rng state

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
    
    log_debug("OpenCL setup complete.");
    
    initialized_=true;
}

std::string I3CLSimStepToPhotonConverterOpenCL::GetPreambleSource()
{
    std::string preamble = I3CLSimHelper::GetMathPreamble(*device_, doublePrecision_);

    // tell the kernel if photons should be stopped once they are detected
    if (stopDetectedPhotons_) {
        preamble = preamble + "#define STOP_PHOTONS_ON_DETECTION\n";
        log_debug("detected photons are stopped");
    } else {
        log_debug("detected photons continue to propagate");
    }

    // Tell the kernel that all photons should be saved.
    // This also turns off all collision detection in the kernel.
    if (saveAllPhotons_) {
        preamble = preamble + "#define SAVE_ALL_PHOTONS\n";
        if (doublePrecision_) {
            preamble = preamble + "#define SAVE_ALL_PHOTONS_PRESCALE " + ToDoubleString(saveAllPhotonsPrescale_) + "\n";
        } else {
            preamble = preamble + "#define SAVE_ALL_PHOTONS_PRESCALE " + ToFloatString(saveAllPhotonsPrescale_) + "\n";
        }
        
    }
    
    
    // should the photon history be saved?
    if (photonHistoryEntries_>0) {
        preamble = preamble + "#define SAVE_PHOTON_HISTORY\n";
        preamble = preamble + "#define NUM_PHOTONS_IN_HISTORY " + boost::lexical_cast<std::string>(photonHistoryEntries_) + "\n";
    }
    
    // Instead of sampling the number of absorption lengths from an
    // exponential distribution with mean 1, use a fixed defined number
    // of absorption lengths for table-making. Photonics uses a weight
    // of 1e-20 corresponding to about 46 absorption lengths.
    if (!std::isnan(fixedNumberOfAbsorptionLengths_)) {
        if (doublePrecision_) {
            preamble = preamble + "#define PROPAGATE_FOR_FIXED_NUMBER_OF_ABSORPTION_LENGTHS " + ToDoubleString(fixedNumberOfAbsorptionLengths_) + "\n";
        } else {
            preamble = preamble + "#define PROPAGATE_FOR_FIXED_NUMBER_OF_ABSORPTION_LENGTHS " + ToFloatString(fixedNumberOfAbsorptionLengths_) + "\n";
        }
    }

    if (pancakeFactor_ != 1.) {
        if (doublePrecision_) {
            preamble = preamble + "#define PANCAKE_FACTOR " + ToDoubleString(pancakeFactor_) + "\n";
        } else {
            preamble = preamble + "#define PANCAKE_FACTOR " + ToFloatString(pancakeFactor_) + "\n";
        }
    }
    
    return preamble;
}

std::string I3CLSimStepToPhotonConverterOpenCL::GetWlenGeneratorSource()
{
    return I3CLSimHelper::GenerateWavelengthGeneratorSource(wlenGenerators_);
}

std::string I3CLSimStepToPhotonConverterOpenCL::GetWlenBiasSource()
{
    return wlenBias_->GetOpenCLFunction("getWavelengthBias"); // name
}

std::string I3CLSimStepToPhotonConverterOpenCL::GetMediumPropertiesSource()
{
    return I3CLSimHelper::GenerateMediumPropertiesSource(*mediumProperties_);
}

std::string I3CLSimStepToPhotonConverterOpenCL::GetGeometrySource()
{
    if (!saveAllPhotons_) {
        return I3CLSimHelper::GenerateGeometrySource(*geometry_,
                                                      geoLayerToOMNumIndexPerStringSetInfo_,
                                                      stringIndexToStringIDBuffer_,
                                                      domIndexToDomIDBuffer_perStringIndex_);
    } else {
        return std::string("");
    }
}

static std::string 
loadKernel(const std::string& name, bool header)
{
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels/";
    const std::string ext = header ? ".h.cl" : ".c.cl";
    return I3CLSimHelper::LoadProgramSource(kernelBaseDir+name+ext);
}

std::string I3CLSimStepToPhotonConverterOpenCL::GetCollisionDetectionSource(bool header)
{
    return loadKernel("sparse_collision_kernel", header);
}

void I3CLSimStepToPhotonConverterOpenCL::Compile()
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    if (compiled_) return; // silently
    
    if (wlenGenerators_.empty())
        throw I3CLSimStepToPhotonConverter_exception("WlenGenerators not set!");
    
    if (!wlenBias_)
        throw I3CLSimStepToPhotonConverter_exception("WlenBias not set!");
    
    if (!mediumProperties_)
        throw I3CLSimStepToPhotonConverter_exception("MediumProperties not set!");
    
    if (!geometry_)
        throw I3CLSimStepToPhotonConverter_exception("Geometry not set!");
    
    if (!deviceIsSelected_)
        throw I3CLSimStepToPhotonConverter_exception("Device not selected!");
    
    if ((saveAllPhotons_) && (stopDetectedPhotons_))
        throw I3CLSimStepToPhotonConverter_exception("Internal error: both the saveAllPhotons and stopDetectedPhotons options are set at the same time.");
    
    prependSource_ = this->GetPreambleSource();
    wlenGeneratorSource_ = this->GetWlenGeneratorSource();
    wlenBiasSource_ = this->GetWlenBiasSource();
    
    mediumPropertiesSource_ = this->GetMediumPropertiesSource();
    
    if (!saveAllPhotons_) {
        geometrySource_ = this->GetGeometrySource();
    } else {
        geometrySource_ = "";
    }
    
    propagationKernelSource_  = loadKernel("propagation_kernel", true);
    if (!saveAllPhotons_) {
        propagationKernelSource_ += this->GetCollisionDetectionSource(true);
        propagationKernelSource_ += this->GetCollisionDetectionSource(false);
    }
    propagationKernelSource_ += loadKernel("propagation_kernel", false);
    
    SetupQueueAndKernel(*(device_->GetPlatformHandle()),
                        *(device_->GetDeviceHandle()));
    
    compiled_=true;
}

std::string I3CLSimStepToPhotonConverterOpenCL::GetFullSource()
{
    std::ostringstream code;
    
    code << prependSource_;
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
            context_ = boost::shared_ptr<cl::Context>(new cl::Context(devices, properties));
        } catch (cl::Error &err) {
            if ((err.err() == CL_OUT_OF_RESOURCES) && (createContextRetriesLeft>0)) {
                --createContextRetriesLeft;
                
                log_warn("Could not create OpenCL context: CL_OUT_OF_RESOURCES. Some drivers will work if we just try again. Waiting %fs.. (%u tries left)",
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
        std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
        log_info("Running on \"%s\"", deviceName.c_str());
        log_debug("      ->                      CL_DEVICE_TYPE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_TYPE>()).c_str());
        log_debug("      ->         CL_DEVICE_MAX_COMPUTE_UNITS: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()).c_str());
        log_debug("      ->    CL_DEVICE_MAX_WORK_ITEM_SIZES[0]: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0]).c_str());
        log_debug("      ->       CL_DEVICE_MAX_WORK_GROUP_SIZE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()).c_str());
        log_debug("      ->       CL_DEVICE_MAX_CLOCK_FREQUENCY: %sMHz", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>()).c_str());
        log_debug("      ->           CL_DEVICE_GLOBAL_MEM_SIZE: %sMiB", boost::lexical_cast<std::string>(static_cast<double>(device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>())/1024./1024.).c_str());
        log_debug("      ->  CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE: %sKiB", boost::lexical_cast<std::string>(static_cast<double>(device.getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>())/1024.).c_str());
        log_debug("      ->            CL_DEVICE_LOCAL_MEM_TYPE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_LOCAL_MEM_TYPE>()).c_str());
        log_debug("      ->            CL_DEVICE_LOCAL_MEM_SIZE: %sKiB", boost::lexical_cast<std::string>(static_cast<double>(device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>())/1024.).c_str());
        log_debug("      ->  CL_DEVICE_ERROR_CORRECTION_SUPPORT: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_ERROR_CORRECTION_SUPPORT>()).c_str());
        log_debug("      ->             CL_DEVICE_ENDIAN_LITTLE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_ENDIAN_LITTLE>()).c_str());
        log_debug("      ->                 CL_DEVICE_AVAILABLE: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_AVAILABLE>()).c_str());
        log_debug("      ->                    CL_DEVICE_VENDOR: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_VENDOR>()).c_str());
        log_debug("      ->                   CL_DEVICE_VERSION: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_VERSION>()).c_str());
        log_debug("      ->                CL_DEVICE_EXTENSIONS: %s", boost::lexical_cast<std::string>(device.getInfo<CL_DEVICE_EXTENSIONS>()).c_str());
    }    
    
    log_debug("Compiling..");
    // accumulate the build options
    std::string BuildOptions;
    
    //BuildOptions += "-w "; // no warnings
    //BuildOptions += "-Werror "; // warnings will become errors
    //BuildOptions += "-cl-opt-disable ";
    //BuildOptions += "-cl-no-signed-zeros ";
    //BuildOptions += "-cl-unsafe-math-optimizations ";
    BuildOptions += "-cl-mad-enable ";

    const bool nvidiaVerboseCompile=false;
    if (nvidiaVerboseCompile)
    {
        // only valid if extension "cl_nv_compiler_options" is present
        BuildOptions += "-cl-nv-verbose ";          // Passed on to ptxas as --verbose
    }
    
    // only valid if extension "cl_nv_compiler_options" is present
    //BuildOptions += "-cl-nv-maxrregcount=60 ";  // Passed on to ptxas as --maxrregcount <N>
    //BuildOptions += "-cl-nv-opt-level=3 ";     // Passed on to ptxas as --opt-level <N>
    
    if (useNativeMath_) {
        BuildOptions += "-cl-fast-relaxed-math ";
        BuildOptions += "-DUSE_NATIVE_MATH ";
    }

    // let the kernel code know that no flasher spectra will be used
    // (might be useful for some optimizations, i.e. less branches)
    if (wlenGenerators_.size() <= 1) {
        BuildOptions += "-DNO_FLASHER ";
    }

    cl::Program program;
    try {
        // build the program

        // combine into a single string first to work around Intel OpenCL
        // compiler issues (as found on OSX 10.11 for example)
        std::string combined_source;
        combined_source += prependSource_ + "\n";
        combined_source += mwcrngKernelSource_ + "\n";
        combined_source += wlenGeneratorSource_ + "\n";
        combined_source += wlenBiasSource_ + "\n";
        combined_source += mediumPropertiesSource_ + "\n";
        if (!saveAllPhotons_) {
            combined_source += geometrySource_ + "\n";
        }
        combined_source += propagationKernelSource_ + "\n";
        
        cl::Program::Sources source;
        source.push_back(std::make_pair(combined_source.c_str(),combined_source.size()));
        
        program = cl::Program(*context_, source);
        log_debug("building...");
        program.build(devices, BuildOptions.c_str());
        log_debug("...building finished.");
        
        if (nvidiaVerboseCompile) {
            std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
#ifdef I3_LOG4CPLUS_LOGGING
            // using LOG_IMPL will make this work even in Release build mode:
            LOG_IMPL(INFO, "  * build status on %s\"", deviceName.c_str());
            LOG_IMPL(INFO, "==============================");
            LOG_IMPL(INFO, "Build Status: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(device)).c_str());
            LOG_IMPL(INFO, "Build Options: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device)).c_str());
            LOG_IMPL(INFO, "Build Log: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)).c_str());
            LOG_IMPL(INFO, "==============================");
#else
            log_info("  * build status on %s\"", deviceName.c_str());
            log_info("==============================");
            log_info("Build Status: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(device)).c_str());
            log_info("Build Options: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device)).c_str());
            log_info("Build Log: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)).c_str());
            log_info("==============================");
#endif
        }
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
    
    const unsigned int numBuffers = disableDoubleBuffering_?1:2;

    // instantiate the command queue
    log_debug("Initializing..");
    try {
        for (unsigned int i=0;i<numBuffers;++i)
        {
#ifdef DUMP_STATISTICS
            queue_.push_back(boost::shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context_, device, CL_QUEUE_PROFILING_ENABLE)));
#else
            queue_.push_back(boost::shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context_, device, 0)));
#endif
        }
    } catch (cl::Error &err) {
        queue_.clear(); // throw away command queue.
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw I3CLSimStepToPhotonConverter_exception("OpenCL error: could not set up command queue!");
    }
    log_debug("initialized.");
    
    // create the kernel
    log_debug("Creating kernel..");
    try {
        // instantiate the kernel object
        for (unsigned int i=0;i<numBuffers;++i)
        {
            kernel_.push_back(boost::shared_ptr<cl::Kernel>(new cl::Kernel(program, "propKernel")));
        }

        maxWorkgroupSize_ = kernel_[0]->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
        
        if (!disableDoubleBuffering_) {
            if (kernel_[1]->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device) != maxWorkgroupSize_) {
                log_fatal("created two identical kernels and got different maximum work group sizes.");
            }
        }
        
        log_debug("Maximum workgroup sizes for the kernel is %" PRIu64, maxWorkgroupSize_);
    } catch (cl::Error &err) {
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

namespace {
#define YIELD_TIME_MICROSECONDS 1000
    
    inline void waitForOpenCLEventsYield(std::vector<cl::Event> &events)
    {
        for (;;)
        {
            bool allDone=true;
            BOOST_FOREACH(cl::Event &event, events)
            {
                if (event.getInfo<CL_EVENT_COMMAND_EXECUTION_STATUS>() != CL_COMPLETE)
                {
                    allDone=false;
                    break;
                }
            }

            if (allDone) break;
            
            // yield
            boost::this_thread::sleep(boost::posix_time::microseconds(YIELD_TIME_MICROSECONDS));
        }
        
        // to be sure, wait for all of them (-> this is a proper synchronization point)
        cl::Event::waitForEvents(events);
    }

    inline void waitForOpenCLEventYield(cl::Event &event)
    {
        for (;;)
        {
            if (event.getInfo<CL_EVENT_COMMAND_EXECUTION_STATUS>() == CL_COMPLETE)
            {
                break;
            }

            // yield
            boost::this_thread::sleep(boost::posix_time::microseconds(YIELD_TIME_MICROSECONDS));
        }
        
        // to be sure, wait for the event again (-> this is a proper synchronization point)
        event.wait();
    }

#undef YIELD_TIME_MILLISECONDS
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
    
    uint32_t stepsIdentifier=0;
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
        waitForOpenCLEventsYield(bufferWriteEvents);
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
                                                  cl::NullRange,    // current implementations force this to be NULL
                                                  cl::NDRange(numberOfInputSteps),  // number of work items
                                                  cl::NDRange(workgroupSize_),
                                                  NULL, //&(bufferWriteEvents),  // wait for buffers to be filled
                                                  &kernelFinishEvent); // signal when finished
        queue_[bufferIndex]->flush(); // make sure it begins executing on the device
    } catch (cl::Error &err) {
        log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
    }

    log_trace("[%u] kernel in queue..", bufferIndex);
}

namespace {
    // converts from the internal photon history fromat (flat array of float4)
    // to a vector of I3CLSimPhotonHistory objects. The output stores photons
    // in forward order (i.e. the most recent scatter listed last)
    I3CLSimPhotonHistorySeriesPtr ConvertPhotonHistories(const std::vector<cl_float4> &rawData,
                                                         const I3CLSimPhotonSeries &photons,
                                                         std::size_t photonHistoryEntries)
    {
        if (rawData.size() % photonHistoryEntries != 0)
            log_fatal("Internal logic error: rawData.size() (==%zu) is not a multiple of photonHistoryEntries (==%zu)",
                     rawData.size(), photonHistoryEntries);
        
        if (rawData.size()/photonHistoryEntries != photons.size())
            log_fatal("internal logic error: rawData.size()/photonHistoryEntries [==%zu/%zu] != photons.size() [==%zu]",
                      rawData.size(),photonHistoryEntries,photons.size());
        
        I3CLSimPhotonHistorySeriesPtr output(new I3CLSimPhotonHistorySeries());
        
        for (std::size_t i=0;i<rawData.size()/photonHistoryEntries;++i)
        {
            // insert a new history for the current photon
            output->push_back(I3CLSimPhotonHistory());
            I3CLSimPhotonHistory &current_history = output->back();
            
            const I3CLSimPhoton &current_photon = photons[i];
            
            const uint32_t numScatters = current_photon.GetNumScatters();
            if (numScatters==0) continue; // nothing to record for this photon
            if (photonHistoryEntries==0) continue; // no entries => nothing to record
            
            const uint32_t numRecordedScatters = static_cast<uint32_t>(std::min(static_cast<std::size_t>(numScatters), photonHistoryEntries));
            // start with the most recent index
            uint32_t currentScatterIndex;
            
            if (numScatters<=photonHistoryEntries) {
                currentScatterIndex=0; // start with index 0
            } else {
                currentScatterIndex=numScatters%photonHistoryEntries;
            }
            
            for (uint32_t j=0;j<numRecordedScatters;++j)
            {
                // [0], [1] and [2] are x,y,z
                // [3] is the distance the photon traveled in units of absorption lengths
                const cl_float4 &rawDataEntry = rawData[i*photonHistoryEntries + currentScatterIndex];
                current_history.push_back(((const cl_float *)&rawDataEntry)[0], ((const cl_float *)&rawDataEntry)[1], ((const cl_float *)&rawDataEntry)[2], ((const cl_float *)&rawDataEntry)[3]);

                ++currentScatterIndex;
                if (currentScatterIndex>=static_cast<uint32_t>(photonHistoryEntries)) currentScatterIndex=0;
            }
        }

        return output;
    }
    
}


void I3CLSimStepToPhotonConverterOpenCL::OpenCLThread_impl_downloadPhotons(boost::this_thread::disable_interruption &di,
                                                                           bool &shouldBreak,
                                                                           unsigned int bufferIndex,
                                                                           uint32_t stepsIdentifier)
{
    shouldBreak=false;
   
    I3CLSimPhotonSeriesPtr photons;
    I3CLSimPhotonHistorySeriesPtr photonHistories;
    boost::shared_ptr<std::vector<cl_float4> > photonHistoriesRaw;
    
    try {
        uint32_t numberOfGeneratedPhotons;
        {
            cl::Event copyComplete;
            queue_[bufferIndex]->enqueueReadBuffer(*deviceBuffer_CurrentNumOutputPhotons[bufferIndex], CL_FALSE, 0, sizeof(uint32_t), &numberOfGeneratedPhotons, NULL, &copyComplete);
            queue_[bufferIndex]->flush(); // make sure it starts executing on the device
            waitForOpenCLEventYield(copyComplete);
        }
        
#ifdef I3_LOG4CPLUS_LOGGING
        LOG_IMPL(INFO, "Num photons to copy (buffer %u): %" PRIu32, bufferIndex, numberOfGeneratedPhotons);
#else
        log_info("Num photons to copy (buffer %u): %" PRIu32, bufferIndex, numberOfGeneratedPhotons);
#endif

#ifdef DUMP_STATISTICS
        {
            boost::unique_lock<boost::mutex> guard(statistics_mutex_);
            statistics_total_num_photons_atDOMs_ += numberOfGeneratedPhotons;
        }
#endif
        
        if (numberOfGeneratedPhotons > maxNumOutputPhotons_)
        {
            log_error("Maximum number of photons exceeded, only receiving %" PRIu32 " of %" PRIu32 " photons",
                      maxNumOutputPhotons_, numberOfGeneratedPhotons);
            numberOfGeneratedPhotons = maxNumOutputPhotons_;
        }
        
        if (numberOfGeneratedPhotons>0)
        {
            VECTOR_CLASS<cl::Event> copyComplete((photonHistoryEntries_>0)?2:1);
            
            // allocate the result vector while waiting for the mapping operation to complete
            photons = I3CLSimPhotonSeriesPtr(new I3CLSimPhotonSeries(numberOfGeneratedPhotons));
            if (photonHistoryEntries_>0) {
                photonHistoriesRaw = boost::shared_ptr<std::vector<cl_float4> >(new std::vector<cl_float4>(numberOfGeneratedPhotons*static_cast<std::size_t>(photonHistoryEntries_)));
            }
            
            queue_[bufferIndex]->enqueueReadBuffer(*deviceBuffer_OutputPhotons[bufferIndex], CL_FALSE, 0, numberOfGeneratedPhotons*sizeof(I3CLSimPhoton), &((*photons)[0]), NULL, &copyComplete[0]);
            
            if (photonHistoryEntries_>0) {
                queue_[bufferIndex]->enqueueReadBuffer(*deviceBuffer_PhotonHistory[bufferIndex], CL_FALSE, 0, numberOfGeneratedPhotons*static_cast<std::size_t>(photonHistoryEntries_)*sizeof(cl_float4), &((*photonHistoriesRaw)[0]), NULL, &copyComplete[1]);
            }
            
            queue_[bufferIndex]->flush(); // make sure it starts executing on the device
            waitForOpenCLEventsYield(copyComplete); // wait for the buffer(s) to be copied

            // convert the histories to the external representation
            if (photonHistoriesRaw) {
                photonHistories = ConvertPhotonHistories(*photonHistoriesRaw, *photons, photonHistoryEntries_);
            }
        }
        else
        {
            // empty vector(s)
            photons = I3CLSimPhotonSeriesPtr(new I3CLSimPhotonSeries());
            if (photonHistoryEntries_>0) {
                photonHistories = I3CLSimPhotonHistorySeriesPtr(new I3CLSimPhotonHistorySeries());
            }
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
            queueFromOpenCL_->Put(ConversionResult_t(stepsIdentifier, photons, photonHistories));
        } catch(boost::thread_interrupted &i) {
            log_debug("OpenCL thread was interrupted. closing.");
            shouldBreak=true;
            return;
        }
    }
    
    
}

boost::posix_time::ptime 
I3CLSimStepToPhotonConverterOpenCL::DumpStatistics(const cl::Event &kernelFinishEvent,
                                                   const boost::posix_time::ptime &last_timestamp,
                                                   uint64_t totalNumberOfPhotons,
                                                   bool starving,
                                                   const std::string &platformName,
                                                   const std::string &deviceName,
                                                   uint64_t deviceProfilingResolution)
{
    // calculate time since last kernel execution
    boost::posix_time::ptime this_timestamp(boost::posix_time::microsec_clock::universal_time());

#ifdef DUMP_STATISTICS
    boost::posix_time::time_duration posix_duration = this_timestamp - last_timestamp;
    const uint64_t host_duration_in_nanoseconds = posix_duration.total_nanoseconds();
    
    uint64_t timeStart, timeEnd;
    kernelFinishEvent.getProfilingInfo(CL_PROFILING_COMMAND_START, &timeStart);
    kernelFinishEvent.getProfilingInfo(CL_PROFILING_COMMAND_END, &timeEnd);
    
    const uint64_t kernel_duration_in_nanoseconds = (timeStart==timeEnd)?deviceProfilingResolution:(timeEnd-timeStart);
    
    {
        boost::unique_lock<boost::mutex> guard(statistics_mutex_);
        
        statistics_total_device_duration_in_nanoseconds_ += kernel_duration_in_nanoseconds;
        statistics_total_host_duration_in_nanoseconds_ += host_duration_in_nanoseconds;
        statistics_total_kernel_calls_++;
        statistics_total_num_photons_generated_ += totalNumberOfPhotons;
    }
    
    const double utilization = static_cast<double>(kernel_duration_in_nanoseconds)/static_cast<double>(host_duration_in_nanoseconds);
    
#ifdef I3_LOG4CPLUS_LOGGING
    // use LOG_IMPL here to make it log this even when in Release build mode.
    LOG_IMPL(INFO, "kernel statistics: %s%g nanoseconds/photon (util: %.0f%%) (%s %s) %s",
             (timeStart==timeEnd)?"<=":"",
             static_cast<double>(kernel_duration_in_nanoseconds)/static_cast<double>(totalNumberOfPhotons),
             utilization*100.,
             platformName.c_str(), deviceName.c_str(),
             (starving?"[starving]":""));
#else
    log_info("kernel statistics: %s%g nanoseconds/photon (util: %.0f%%) (%s %s) %s",
             (timeStart==timeEnd)?"<=":"",
             static_cast<double>(kernel_duration_in_nanoseconds)/static_cast<double>(totalNumberOfPhotons),
             utilization*100.,
             platformName.c_str(), deviceName.c_str(),
             (starving?"[starving]":""));
#endif
#endif
    
    return this_timestamp;
}

void I3CLSimStepToPhotonConverterOpenCL::OpenCLThread_impl(boost::this_thread::disable_interruption &di)
{
    // set things up here
    if (!context_) log_fatal("Internal error: context is (null)");

    const std::size_t numBuffers = disableDoubleBuffering_?1:2;
    
    if (queue_.size() != numBuffers) log_fatal("Internal error: queue_.size() != 2!");
    if (kernel_.size() != numBuffers) log_fatal("Internal error: kernel_.size() != 2!");

    BOOST_FOREACH(boost::shared_ptr<cl::CommandQueue> &ptr, queue_) {
        if (!ptr) log_fatal("Internal error: queue_[] is (null)");
    }
    BOOST_FOREACH(boost::shared_ptr<cl::Kernel> &ptr, kernel_) {
        if (!ptr) log_fatal("Internal error: kernel_[] is (null)");
    }

    if (deviceBuffer_InputSteps.size() != numBuffers) log_fatal("Internal error: deviceBuffer_InputSteps.size() != 2!");
    if (deviceBuffer_OutputPhotons.size() != numBuffers) log_fatal("Internal error: deviceBuffer_OutputPhotons.size() != 2!");
    if (deviceBuffer_CurrentNumOutputPhotons.size() != numBuffers) log_fatal("Internal error: deviceBuffer_CurrentNumOutputPhotons.size() != 2!");
    if (photonHistoryEntries_ > 0) {
        if (deviceBuffer_PhotonHistory.size() != numBuffers) log_fatal("Internal error: deviceBuffer_PhotonHistory.size() != 2!");
    }
    
    BOOST_FOREACH(boost::shared_ptr<cl::Buffer> &ptr, deviceBuffer_InputSteps) {
        if (!ptr) log_fatal("Internal error: deviceBuffer_InputSteps[] is (null)");
    }
    BOOST_FOREACH(boost::shared_ptr<cl::Buffer> &ptr, deviceBuffer_OutputPhotons) {
        if (!ptr) log_fatal("Internal error: deviceBuffer_OutputPhotons[] is (null)");
    }
    BOOST_FOREACH(boost::shared_ptr<cl::Buffer> &ptr, deviceBuffer_CurrentNumOutputPhotons) {
        if (!ptr) log_fatal("Internal error: deviceBuffer_CurrentNumOutputPhotons[] is (null)");
    }

    if (photonHistoryEntries_ > 0) {
        BOOST_FOREACH(boost::shared_ptr<cl::Buffer> &ptr, deviceBuffer_PhotonHistory) {
            if (!ptr) log_fatal("Internal error: deviceBuffer_PhotonHistory[] is (null)");
        }
    }

    if (!saveAllPhotons_) {
        if (!deviceBuffer_GeoLayerToOMNumIndexPerStringSet) log_fatal("Internal error: deviceBuffer_GeoLayerToOMNumIndexPerStringSet is (null)");
    }
    if (!deviceBuffer_MWC_RNG_x) log_fatal("Internal error: deviceBuffer_MWC_RNG_x is (null)");
    if (!deviceBuffer_MWC_RNG_a) log_fatal("Internal error: deviceBuffer_MWC_RNG_a is (null)");
    
    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(openCLStarted_mutex_);
        openCLStarted_=true;
    }
    openCLStarted_cond_.notify_all();
    
    std::vector<uint32_t> stepsIdentifier(numBuffers, 0);
    std::vector<uint64_t> totalNumberOfPhotons(numBuffers, 0);
    std::vector<std::size_t> numberOfSteps(numBuffers, 0);
    
#ifdef DUMP_STATISTICS
    boost::posix_time::ptime last_timestamp(boost::posix_time::microsec_clock::universal_time());
#endif
    
    unsigned int thisBuffer=0;
    unsigned int otherBuffer=1;
    bool otherBufferHasBeenCopied=false;
    
    if (!disableDoubleBuffering_) {
        // swap buffers once to swap them back just a few lines later
        std::swap(thisBuffer, otherBuffer);
    }
    
    // start the main loop
    for (;;)
    {
        if (!disableDoubleBuffering_) {
            // swap buffers
            std::swap(thisBuffer, otherBuffer);
        }
        
        log_trace("buffers indices now: this==%u, other==%u", thisBuffer, otherBuffer);
        
        bool starving=false;
        if ((!otherBufferHasBeenCopied) || (disableDoubleBuffering_)) {
            if (!disableDoubleBuffering_) starving=true;
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

        if (!disableDoubleBuffering_)
        {
            // if there already is a new buffer available, copy it now, while the kernel is running
            log_trace("[%u] Starting copy (other buffer)..", otherBuffer);

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
            waitForOpenCLEventYield(kernelFinishEvent);
        } catch (cl::Error &err) {
            log_fatal("[%u] OpenCL ERROR (running kernel): %s (%i)", thisBuffer, err.what(), err.err());
        }

        log_trace("[%u] kernel finished..", thisBuffer);

#ifdef DUMP_STATISTICS
        log_trace("[%u] dumping statistics..", thisBuffer);

        last_timestamp = DumpStatistics(kernelFinishEvent,
                                        last_timestamp,
                                        totalNumberOfPhotons[thisBuffer],
                                        starving,
                                        device_->GetPlatformName(),
                                        device_->GetDeviceName(),
                                        (device_->GetDeviceHandle())->getInfo<CL_DEVICE_PROFILING_TIMER_RESOLUTION>() );
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

void I3CLSimStepToPhotonConverterOpenCL::SetEnableDoubleBuffering(bool value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");

    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    disableDoubleBuffering_=(!value);
}

bool I3CLSimStepToPhotonConverterOpenCL::GetEnableDoubleBuffering() const
{
    return (!disableDoubleBuffering_);
}



void I3CLSimStepToPhotonConverterOpenCL::SetDoublePrecision(bool value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    doublePrecision_=value;
}

bool I3CLSimStepToPhotonConverterOpenCL::GetDoublePrecision() const
{
    return doublePrecision_;
}


void I3CLSimStepToPhotonConverterOpenCL::SetStopDetectedPhotons(bool value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    if ((value) && (saveAllPhotons_))
        throw I3CLSimStepToPhotonConverter_exception("You cannot set stopDetectedPhotons, because saveAllPhotons is set. The options are mutually exclusive.");

    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    stopDetectedPhotons_=value;
}

bool I3CLSimStepToPhotonConverterOpenCL::GetStopDetectedPhotons() const
{
    return stopDetectedPhotons_;
}



void I3CLSimStepToPhotonConverterOpenCL::SetSaveAllPhotons(bool value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    if ((value) && (stopDetectedPhotons_))
        throw I3CLSimStepToPhotonConverter_exception("You cannot set saveAllPhotons, because stopDetectedPhotons is set. The options are mutually exclusive.");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    saveAllPhotons_=value;
}

bool I3CLSimStepToPhotonConverterOpenCL::GetSaveAllPhotons() const
{
    return saveAllPhotons_;
}



void I3CLSimStepToPhotonConverterOpenCL::SetSaveAllPhotonsPrescale(double value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    saveAllPhotonsPrescale_=value;
}

double I3CLSimStepToPhotonConverterOpenCL::GetSaveAllPhotonsPrescale() const
{
    return saveAllPhotonsPrescale_;
}



void I3CLSimStepToPhotonConverterOpenCL::SetPhotonHistoryEntries(uint32_t value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    photonHistoryEntries_=value;
}

uint32_t I3CLSimStepToPhotonConverterOpenCL::GetPhotonHistoryEntries() const
{
    return photonHistoryEntries_;
}


void I3CLSimStepToPhotonConverterOpenCL::SetFixedNumberOfAbsorptionLengths(double value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    fixedNumberOfAbsorptionLengths_=value;
}

double I3CLSimStepToPhotonConverterOpenCL::GetFixedNumberOfAbsorptionLengths() const
{
    return fixedNumberOfAbsorptionLengths_;
}


void I3CLSimStepToPhotonConverterOpenCL::SetDOMPancakeFactor(double value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    pancakeFactor_=value;
}

double I3CLSimStepToPhotonConverterOpenCL::GetDOMPancakeFactor() const
{
    return pancakeFactor_;
}



void I3CLSimStepToPhotonConverterOpenCL::SetWlenGenerators(const std::vector<I3CLSimRandomValueConstPtr> &wlenGenerators)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    
    compiled_=false;
    kernel_.clear();
    queue_.clear();
    
    wlenGenerators_=wlenGenerators;
}

void I3CLSimStepToPhotonConverterOpenCL::SetWlenBias(I3CLSimFunctionConstPtr wlenBias)
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
    
    if ((result.photons) && (!saveAllPhotons_)) {
        ReplaceStringDOMIndexWithStringDOMIDs(*result.photons,
                                              stringIndexToStringIDBuffer_,
                                              domIndexToDomIDBuffer_perStringIndex_);
        
    }
    
    return result;
}
