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
 * @file I3CLSimOpenCLDevice.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <clsim/I3CLSimOpenCLDevice.h>

#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#define __CL_ENABLE_EXCEPTIONS
#include "clsim/cl.hpp"


// static data
bool I3CLSimOpenCLDevice::lists_initialized_=false;
std::vector<I3CLSimOpenCLDevice> I3CLSimOpenCLDevice::allDevices_;


// construction & destruction
I3CLSimOpenCLDevice::I3CLSimOpenCLDevice(const std::string &platformName,
                                         const std::string &deviceName)
{
    InitializeStaticStuff();
    DoInit(platformName, deviceName);
    
}

I3CLSimOpenCLDevice::I3CLSimOpenCLDevice(const std::string &platformName,
                                         const std::string &deviceName,
                                         bool useNativeMath,
                                         uint32_t approximateNumberOfWorkItems)
{
    InitializeStaticStuff();
    DoInit(platformName, deviceName);

    useNativeMath_=useNativeMath;
    approximateNumberOfWorkItems_=approximateNumberOfWorkItems;
}

void I3CLSimOpenCLDevice::DoInit(const std::string &platformName,
                                 const std::string &deviceName)
{
    std::size_t deviceIndex=0;
    bool deviceFound=false;
    
    for (std::size_t i=0; i<allDevices_.size(); ++i)
    {
        if ((allDevices_[i].platformName_ == platformName) &&
            (allDevices_[i].deviceName_ == deviceName))
        {
            deviceIndex=i;
            deviceFound=true;
            break;
        }
    }
    
    if (!deviceFound)
    {
        log_error("The selected OpenCL device was not found! [platform=\"%s\", device=\"%s\"]",
                  platformName.c_str(), deviceName.c_str());
        log_error("Here is a list of available devices:");
        for (std::size_t i=0; i<allDevices_.size(); ++i)
        {
            log_error("platform: \"%s\", device: \"%s\"",
                      allDevices_[i].platformName_.c_str(), allDevices_[i].deviceName_.c_str());
            
        }
        log_fatal("Could not find selected OpenCL device.");
    }
    else
    {
        // device found!
        
        // initialize from template
        (*this) = allDevices_[deviceIndex];
        
        if ((!platform_) || (!device_))
            log_fatal("Internal error: platform or device is NULL");
    }
}

I3CLSimOpenCLDevice::I3CLSimOpenCLDevice()
:
useNativeMath_(false),
approximateNumberOfWorkItems_(1024)
{
    
}

I3CLSimOpenCLDevice::~I3CLSimOpenCLDevice() 
{ 
    
}

I3CLSimOpenCLDeviceSeriesPtr I3CLSimOpenCLDevice::GetAllDevices()
{
    InitializeStaticStuff();

    I3CLSimOpenCLDeviceSeriesPtr retval(new I3CLSimOpenCLDeviceSeries(allDevices_));

    return retval;
}

// comparison
bool operator==(const I3CLSimOpenCLDevice &a, const I3CLSimOpenCLDevice &b)
{
    if (((!a.platform_) && (b.platform_)) || ((a.platform_) && (!b.platform_))) return false;
    if (((!a.device_) && (b.device_)) || ((a.device_) && (!b.device_))) return false;
    
    // currently compares the pointers only (not what they point at)
    if (a.platform_) {
        if (a.platform_ != b.platform_) return false;
    }

    if (a.device_) {
        if (a.device_ != b.device_) return false;
    }

    if (a.useNativeMath_ != b.useNativeMath_) return false;
    if (a.approximateNumberOfWorkItems_ != b.approximateNumberOfWorkItems_) return false;
    
    return true;
}

// device information
bool I3CLSimOpenCLDevice::IsCPU() const {return (device_->getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);}
bool I3CLSimOpenCLDevice::IsGPU() const {return (device_->getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_GPU);}
std::size_t I3CLSimOpenCLDevice::GetMaxComputeUnits() const {return device_->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();}
std::size_t I3CLSimOpenCLDevice::GetMaxWorkItemSize() const {return device_->getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0];}
std::size_t I3CLSimOpenCLDevice::GetMaxWorkGroupSize() const {return device_->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();}
std::size_t I3CLSimOpenCLDevice::GetMaxClockFrequencyMhz() const {return device_->getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>();}
std::size_t I3CLSimOpenCLDevice::GetGlobalMemSize() const {return device_->getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();}
std::size_t I3CLSimOpenCLDevice::GetMaxConstantBufferSize() const {return device_->getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>();}
std::size_t I3CLSimOpenCLDevice::GetLocalMemSize() const {return device_->getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();}
bool I3CLSimOpenCLDevice::HasDedicatedLocalMem() const {return (device_->getInfo<CL_DEVICE_LOCAL_MEM_TYPE>() == CL_LOCAL);}
bool I3CLSimOpenCLDevice::HasErrorCorrectionSupport() const {return device_->getInfo<CL_DEVICE_ERROR_CORRECTION_SUPPORT>();}
bool I3CLSimOpenCLDevice::IsAvailable() const {return device_->getInfo<CL_DEVICE_AVAILABLE>();}
std::string I3CLSimOpenCLDevice::GetVendor() const {return device_->getInfo<CL_DEVICE_VENDOR>();}
std::string I3CLSimOpenCLDevice::GetDriverVersion() const {return device_->getInfo<CL_DRIVER_VERSION>();}
std::string I3CLSimOpenCLDevice::GetDeviceVersion() const {return device_->getInfo<CL_DEVICE_VERSION>();}
std::string I3CLSimOpenCLDevice::GetExtensions() const {return device_->getInfo<CL_DEVICE_EXTENSIONS>();}



// static initialization
void I3CLSimOpenCLDevice::InitializeStaticStuff()
{
    if (lists_initialized_) return;
    
    // enumerate platforms and devices
    std::vector<std::pair<std::string, std::string> > deviceNameList;
    std::vector<std::pair<shared_ptr<cl::Platform>, shared_ptr<cl::Device> > > clPlatformDeviceList;
    
    std::vector<cl::Platform> platforms;
    
    try {
        // get a list of available platforms
        cl::Platform::get(&platforms);
    } catch (cl::Error &err) {
        log_fatal("OpenCL ERROR: %s (%i)", err.what(), err.err());
    }
    
    BOOST_FOREACH(cl::Platform &platform, platforms)
    {
        const std::string platformName = platform.getInfo<CL_PLATFORM_NAME>();
        
        std::vector<cl::Device> devices;
        
        try {
            platform.getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR: %s (%i)", err.what(), err.err());
        }
        
        BOOST_FOREACH(cl::Device &device, devices)
        {
            const std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
            
            deviceNameList.push_back(std::make_pair(platformName, deviceName));
            clPlatformDeviceList.push_back(std::make_pair(shared_ptr<cl::Platform>(new cl::Platform(platform)), shared_ptr<cl::Device>(new cl::Device(device))));
            
            log_trace("raw: PLATFORM: \"%s\" -> DEVICE: \"%s\"",
                      platformName.c_str(),
                      deviceName.c_str());
        }
    }
    
    // make sure devices with exactly the same platform&deviceName pair get distinct names
    if (deviceNameList.size() > 1)
    {
        std::vector<bool> checked(deviceNameList.size(), false);
        for (std::size_t i=0; i<deviceNameList.size()-1; ++i)
        {
            std::size_t occurence=0;
            for (std::size_t j=i+1; j<deviceNameList.size(); ++j)
            {
                if (checked[j]) continue;
                
                if (deviceNameList[i].first != deviceNameList[j].first) continue; // different platforms
                if (deviceNameList[i].second != deviceNameList[j].second) continue; // different device names
                ++occurence; // found one
                
                deviceNameList[j].second += " (" + boost::lexical_cast<std::string>(occurence+1) + ")";
                checked[j] = true;
            }
            
            if (occurence>0) {
                deviceNameList[i].second += " (1)";
            }
            
            checked[i]=true;
        }
    }
    
    for (std::size_t i=0;i<deviceNameList.size();++i)
    {
        log_debug("PLATFORM: \"%s\" -> DEVICE: \"%s\"",
                  deviceNameList[i].first.c_str(),
                  deviceNameList[i].second.c_str());
    }
    
    // generate template I3CLSimOpenCLDevice objects
    allDevices_.clear();
    for (std::size_t i=0;i<deviceNameList.size();++i)
    {
        allDevices_.push_back(I3CLSimOpenCLDevice());
        I3CLSimOpenCLDevice &newDevice = allDevices_.back();
        
        newDevice.platformName_ = deviceNameList[i].first;
        newDevice.deviceName_ = deviceNameList[i].second;
        
        newDevice.platform_ = clPlatformDeviceList[i].first;
        newDevice.device_ = clPlatformDeviceList[i].second;
    }
    
    
    lists_initialized_=true;
}
