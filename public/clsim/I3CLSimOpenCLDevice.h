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
 * @file I3CLSimOpenCLDevice.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMOPENCLDEVICE_H_INCLUDED
#define I3CLSIMOPENCLDEVICE_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

// forward declarations
namespace cl {
    class Platform;
    class Device;
};

/**
 * @brief Describes an OpenCL platform, device and
 * a few device-specific parameters
 */
class I3CLSimOpenCLDevice
{
public:
    ~I3CLSimOpenCLDevice();

    I3CLSimOpenCLDevice(const std::string &platformName,
                        const std::string &deviceName);

    I3CLSimOpenCLDevice(const std::string &platformName,
                        const std::string &deviceName,
                        bool useNativeMath,
                        uint32_t approximateNumberOfWorkItems);

    static boost::shared_ptr<std::vector<I3CLSimOpenCLDevice> > GetAllDevices();
    
    // split a device using "OpenCL device fission"
    // (might be useful to get a single core of a CPU)
    boost::shared_ptr<std::vector<I3CLSimOpenCLDevice> > SplitDevice() const;
    
    inline const std::string &GetPlatformName() const {return platformName_;}
    inline const std::string &GetDeviceName() const {return deviceName_;}

    // device settings
    inline void SetUseNativeMath(bool useit) {useNativeMath_=useit;}
    inline bool GetUseNativeMath() const {return useNativeMath_;}
    inline void SetApproximateNumberOfWorkItems(uint32_t newnumber) {approximateNumberOfWorkItems_=newnumber;}
    inline uint32_t GetApproximateNumberOfWorkItems() const {return approximateNumberOfWorkItems_;}

    // device properties (from OpenCL)
    bool IsCPU() const;
    bool IsGPU() const;
    std::size_t GetMaxComputeUnits() const;
    std::size_t GetMaxWorkItemSize() const;
    std::size_t GetMaxWorkGroupSize() const;
    std::size_t GetMaxClockFrequencyMhz() const;
    std::size_t GetGlobalMemSize() const;
    std::size_t GetMaxConstantBufferSize() const;
    std::size_t GetLocalMemSize() const;
    bool HasDedicatedLocalMem() const;
    bool HasErrorCorrectionSupport() const;
    bool IsAvailable() const;
    std::string GetVendor() const;
    std::string GetDriverVersion() const;
    std::string GetDeviceVersion() const;
    std::string GetExtensions() const;
    
    // get platform and device handles
    inline boost::shared_ptr<cl::Platform> GetPlatformHandle() const {return platform_;}
    inline boost::shared_ptr<cl::Device> GetDeviceHandle() const {return device_;}
    
private:
    I3CLSimOpenCLDevice(); // no default construction
    void DoInit(const std::string &platformName,
                const std::string &deviceName);
    
    boost::shared_ptr<cl::Platform> platform_;
    boost::shared_ptr<cl::Device> device_;

    std::string platformName_;
    std::string deviceName_;
    
    bool useNativeMath_;
    uint32_t approximateNumberOfWorkItems_;
    
private: // static stuff
    
    static bool lists_initialized_;
    static void InitializeStaticStuff();
    
    static std::vector<I3CLSimOpenCLDevice> allDevices_;
    
    friend bool operator==(const I3CLSimOpenCLDevice &, const I3CLSimOpenCLDevice &);
};
bool operator==(const I3CLSimOpenCLDevice &a, const I3CLSimOpenCLDevice &b);

typedef std::vector<I3CLSimOpenCLDevice> I3CLSimOpenCLDeviceSeries;

I3_POINTER_TYPEDEFS(I3CLSimOpenCLDevice);
I3_POINTER_TYPEDEFS(I3CLSimOpenCLDeviceSeries);

#endif //I3CLSIMOPENCLDEVICE_H_INCLUDED
