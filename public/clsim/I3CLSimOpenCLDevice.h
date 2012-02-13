//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module clsim
//
//   clsim is free software; you can redistribute it and/or modify
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

    static shared_ptr<std::vector<I3CLSimOpenCLDevice> > GetAllDevices();
    
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
    inline shared_ptr<cl::Platform> GetPlatformHandle() const {return platform_;}
    inline shared_ptr<cl::Device> GetDeviceHandle() const {return device_;}
    
private:
    I3CLSimOpenCLDevice(); // no default construction
    void DoInit(const std::string &platformName,
                const std::string &deviceName);
    
    shared_ptr<cl::Platform> platform_;
    shared_ptr<cl::Device> device_;

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
