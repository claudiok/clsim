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
 * @file I3CLSimTesterBase.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "test/I3CLSimTesterBase.h"

#include <string>
#include <sstream>
#include <algorithm>
#include <limits>

#include <stdlib.h>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "opencl/I3CLSimHelperLoadProgramSource.h"

#define __CL_ENABLE_EXCEPTIONS
#include "clsim/cl.hpp"

I3CLSimTesterBase::I3CLSimTesterBase()
{
    ;
}


void I3CLSimTesterBase::DoSetup(const I3CLSimOpenCLDevice &device,
                                uint64_t workgroupSize_,
                                uint64_t workItemsPerIteration_,
                                const std::vector<std::string> &source,
                                const std::string compilerOptions)
{
    log_debug("I3CLSimTesterBase::DoSetup()");
    
    workgroupSize=workgroupSize_;
    workItemsPerIteration=workItemsPerIteration_;

    if (workgroupSize==0) {
        log_error("workgroupSize must not be 0!");
        throw std::runtime_error("workgroupSize must not be 0!");
    }
    
    if (workItemsPerIteration==0) {
        log_error("workItemsPerIteration must not be 0!");
        throw std::runtime_error("workItemsPerIteration must not be 0!");
    }

    if (workItemsPerIteration%workgroupSize != 0) {
        log_error("workItemsPerIteration is not a multiple of workgroupSize.");
        throw std::runtime_error("workItemsPerIteration is not a multiple of workgroupSize.");
    }
    
    
    // prepend the rng to the sources list
    sourceStrings_.clear();
    cl::Program::Sources source_;

    std::string combined_source;
    
    // copy the input source strings
    BOOST_FOREACH(const std::string &src, source)
    {
        combined_source += src + "\n";
    }

    sourceStrings_.push_back(combined_source);
    source_.push_back(std::make_pair(sourceStrings_.back().c_str(),sourceStrings_.back().size()));
    
    // get the device object
    boost::shared_ptr<cl::Platform> platformHandle = device.GetPlatformHandle();
    boost::shared_ptr<cl::Device> deviceHandle = device.GetDeviceHandle();
    
    // initialize things
    boost::shared_ptr<std::vector<cl::Device> > devices(new std::vector<cl::Device>(1, *deviceHandle));
    
    try {
        // prepare a device vector (containing a single device)
        cl_context_properties properties[] = 
        { CL_CONTEXT_PLATFORM, (cl_context_properties)(*platformHandle)(), 0};
        
        // create a context
        context = boost::shared_ptr<cl::Context>(new cl::Context(*devices, properties));
    } catch (cl::Error &err) {
        log_error("OpenCL error: could not set up context!");
        throw std::runtime_error("OpenCL error: could not set up context!");
    }

    log_debug("Compiling..");
    // accumulate the build options
    std::string BuildOptions;
    
    //BuildOptions += "-cl-mad-enable ";
    //BuildOptions += "-cl-fast-relaxed-math ";
    
    // only valid if extension "cl_nv_compiler_options" is present
    //BuildOptions += "-cl-nv-verbose ";          // Passed on to ptxas as --verbose
    //BuildOptions += "-cl-nv-maxrregcount=60 ";  // Passed on to ptxas as --maxrregcount <N>
    //BuildOptions += "-cl-nv-opt-level=3 ";     // Passed on to ptxas as --opt-level <N>
    if (device.GetUseNativeMath()) {BuildOptions += "-DUSE_NATIVE_MATH ";}
    
    BuildOptions += compilerOptions;
    
    try {
        program = boost::shared_ptr<cl::Program>(new cl::Program(*context, source_));
        log_debug("building...");
        program->build(*devices, BuildOptions.c_str());
        log_debug("...building finished.");
    } catch (cl::Error &err) {
        log_error("OpenCL ERROR (compile): %s (%i)", err.what(), err.err());
        
        std::string deviceName = deviceHandle->getInfo<CL_DEVICE_NAME>();
        log_error("  * build status on %s\"", deviceName.c_str());
        log_error("==============================");
        log_error("Build Status: %s", boost::lexical_cast<std::string>(program->getBuildInfo<CL_PROGRAM_BUILD_STATUS>(*deviceHandle)).c_str());
        log_error("Build Options: %s", boost::lexical_cast<std::string>(program->getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(*deviceHandle)).c_str());
        log_error("Build Log: %s", boost::lexical_cast<std::string>(program->getBuildInfo<CL_PROGRAM_BUILD_LOG>(*deviceHandle)).c_str());
        log_error("==============================");
        
        throw std::runtime_error("OpenCL error: could build the OpenCL program!");;
    }
    log_debug("code compiled.");
    
    // instantiate the command queue
    log_debug("Initializing..");
    try {
        //queue_ = boost::shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context_, device.GetDeviceHandle(), CL_QUEUE_PROFILING_ENABLE));
        queue = boost::shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context, *deviceHandle, 0));
    } catch (cl::Error &err) {
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw std::runtime_error("OpenCL error: could not set up command queue!");
    }
    log_debug("initialized.");
    
    // create the kernel
    log_debug("Creating kernel..");
    
    try {
        // instantiate the kernel object
        kernel = boost::shared_ptr<cl::Kernel>(new cl::Kernel(*program, "testKernel"));
        maxWorkgroupSize = kernel->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(*deviceHandle);
        
        log_debug("Maximum workgroup sizes for the kernel is %" PRIu64, maxWorkgroupSize);
    } catch (cl::Error &err) {
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw std::runtime_error("OpenCL error: could not create kernel!");
    }
    log_debug("created.");

    if (workgroupSize > maxWorkgroupSize)
    {
        log_error("Requested workgroup size is too large, maximum is %" PRIu64,
                  maxWorkgroupSize);
        std::string message("Requested workgroup size is too large, maximum is " + boost::lexical_cast<std::string>(maxWorkgroupSize));
        throw std::runtime_error(message.c_str());
    }
    
    
}

