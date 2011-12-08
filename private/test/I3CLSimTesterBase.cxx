//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module g4sim-interface.
//
//   this file is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

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

std::vector<std::pair<std::string, std::string> > I3CLSimTesterBase::GetDeviceNameList()
{
    std::vector<std::pair<std::string, std::string> > deviceNameList;
    
    // enumerate platforms and devices
    std::vector<cl::Platform> platforms;
    
    try {
        // get a list of available platforms
        cl::Platform::get(&platforms);
    } catch (cl::Error &err) {
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw std::runtime_error("OpenCL error.");
    }
    
    BOOST_FOREACH(cl::Platform &platform, platforms)
    {
        const std::string platformName = platform.getInfo<CL_PLATFORM_NAME>();
        
        std::vector<cl::Device> devices;
        
        try {
            platform.getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);
        } catch (cl::Error &err) {
            // an OpenCL error here most probably means that there are no devices of the
            // requested type. So just continue quietly.
            log_debug("OpenCL ERROR: %s (%i)", err.what(), err.err());
            continue;
        }
        
        BOOST_FOREACH(cl::Device &device, devices)
        {
            const std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
            
            deviceNameList.push_back(std::make_pair(platformName, deviceName));
            
            log_trace("PLATFORM: \"%s\" -> DEVICE: \"%s\"",
                      platformName.c_str(),
                      deviceName.c_str());
        }
    }
    
    return deviceNameList;
}


std::pair<cl::Platform, cl::Device> I3CLSimTesterBase::GetPlatformDeviceFromNames(const std::string &platformName, const std::string &deviceName) const
{
    // enumerate platforms and devices
    std::vector<cl::Platform> platforms;
    
    try {
        // get a list of available platforms
        cl::Platform::get(&platforms);
    } catch (cl::Error &err) {
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw std::runtime_error("OpenCL error.");
    }
    
    BOOST_FOREACH(cl::Platform &platform, platforms)
    {
        const std::string thisPlatformName = platform.getInfo<CL_PLATFORM_NAME>();
        
        if (thisPlatformName!=platformName) continue;
        
        std::vector<cl::Device> devices;
        
        try {
            platform.getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);
        } catch (cl::Error &err) {
            // an OpenCL error here most probably means that there are no devices of the
            // requested type. So just continue quietly.
            log_debug("OpenCL ERROR: %s (%i)", err.what(), err.err());
            continue;
        }
        
        BOOST_FOREACH(cl::Device &device, devices)
        {
            const std::string thisDeviceName = device.getInfo<CL_DEVICE_NAME>();
            
            if (thisDeviceName != deviceName) continue;
            
            // found it!
            return std::make_pair(platform, device);
        }
    }
    
    log_error("Platform \"%s\"/Device \"%s\" not found!",
              platformName.c_str(), deviceName.c_str());
    
    return std::make_pair(cl::Platform(), cl::Device());
}


void I3CLSimTesterBase::DoSetup(const std::pair<std::string, std::string> &platformAndDeviceName,
                                bool useNativeMath,
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
    
    // copy the input source strings
    BOOST_FOREACH(const std::string &src, source)
    {
        sourceStrings_.push_back(src);
        source_.push_back(std::make_pair(sourceStrings_.back().c_str(),sourceStrings_.back().size()));
    }
    
    // get the device object
    std::pair<cl::Platform, cl::Device> platformAndDevice = GetPlatformDeviceFromNames(platformAndDeviceName.first, platformAndDeviceName.second);
    cl::Platform &platform = platformAndDevice.first;
    cl::Device &device = platformAndDevice.second;
    
    // initialize things
    devices = shared_ptr<std::vector<cl::Device> >(new std::vector<cl::Device>(1, device));
    
    try {
        // prepare a device vector (containing a single device)
        cl_context_properties properties[] = 
        { CL_CONTEXT_PLATFORM, (cl_context_properties)(platform)(), 0};
        
        // create a context
        context = shared_ptr<cl::Context>(new cl::Context(*devices, properties));
    } catch (cl::Error &err) {
        log_error("OpenCL error: could not set up context!");
        throw std::runtime_error("OpenCL error: could not set up context!");
    }
    
    {
        log_debug("Running on: ");
        std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
        log_debug("  * \"%s\"", deviceName.c_str());
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
    
    //BuildOptions += "-cl-mad-enable ";
    //BuildOptions += "-cl-fast-relaxed-math ";
    
    // only valid if extension "cl_nv_compiler_options" is present
    //BuildOptions += "-cl-nv-verbose ";          // Passed on to ptxas as --verbose
    //BuildOptions += "-cl-nv-maxrregcount=60 ";  // Passed on to ptxas as --maxrregcount <N>
    //BuildOptions += "-cl-nv-opt-level=3 ";     // Passed on to ptxas as --opt-level <N>
    if (useNativeMath) {BuildOptions += "-DUSE_NATIVE_MATH ";}
    
    BuildOptions += compilerOptions;
    
    try {
        program = shared_ptr<cl::Program>(new cl::Program(*context, source_));
        log_debug("building...");
        program->build(*devices, BuildOptions.c_str());
        log_debug("...building finished.");
    } catch (cl::Error &err) {
        log_error("OpenCL ERROR (compile): %s (%i)", err.what(), err.err());
        
        std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
        log_error("  * build status on %s\"", deviceName.c_str());
        log_error("==============================");
        log_error("Build Status: %s", boost::lexical_cast<std::string>(program->getBuildInfo<CL_PROGRAM_BUILD_STATUS>(device)).c_str());
        log_error("Build Options: %s", boost::lexical_cast<std::string>(program->getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device)).c_str());
        log_error("Build Log: %s", boost::lexical_cast<std::string>(program->getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)).c_str());
        log_error("==============================");
        
        throw std::runtime_error("OpenCL error: could build the OpenCL program!");;
    }
    log_debug("code compiled.");
    
    // instantiate the command queue
    log_debug("Initializing..");
    try {
        //queue_ = shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context_, device, CL_QUEUE_PROFILING_ENABLE));
        queue = shared_ptr<cl::CommandQueue>(new cl::CommandQueue(*context, device, 0));
    } catch (cl::Error err) {
        log_error("OpenCL ERROR: %s (%i)", err.what(), err.err());
        throw std::runtime_error("OpenCL error: could not set up command queue!");
    }
    log_debug("initialized.");
    
    // create the kernel
    log_debug("Creating kernel..");
    
    try {
        // instantiate the kernel object
        kernel = shared_ptr<cl::Kernel>(new cl::Kernel(*program, "testKernel"));
        maxWorkgroupSize = kernel->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
        
        log_debug("Maximum workgroup sizes for the kernel is %" PRIu64, maxWorkgroupSize);
    } catch (cl::Error err) {
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

