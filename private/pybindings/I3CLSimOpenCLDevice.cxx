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

#include <sstream>

#include "clsim/I3CLSimOpenCLDevice.h"

#include <icetray/python/list_indexing_suite.hpp>
#include <icetray/python/copy_suite.hpp>

using namespace boost::python;
namespace bp = boost::python;

bool I3CLSimOpenCLDevice_equalWrap(const I3CLSimOpenCLDevice &a, const I3CLSimOpenCLDevice &b)
{
    return a==b;
}

static std::string 
I3CLSimOpenCLDevice_repr(const I3CLSimOpenCLDevice& s)
{
    std::ostringstream oss;

    oss << "I3CLSimOpenCLDevice('" << s.GetPlatformName() << "', '" << s.GetDeviceName() << "', useNativeMath=" << (s.GetUseNativeMath()?"True":"False") << ", approximateNumberOfWorkItems=" << s.GetApproximateNumberOfWorkItems() << ")";
    
    return oss.str();
}

static std::string 
I3CLSimOpenCLDevice_prettyprint(const I3CLSimOpenCLDevice& s)
{
    std::ostringstream oss;

    oss << "[ I3CLSimOpenCLDevice platform : "  << s.GetPlatformName() << std::endl
        << "                        device : "  << s.GetDeviceName() << std::endl
        << std::endl
        << "                           cpu : "  << (s.IsCPU()?"YES":"NO") << std::endl
        << "                           gpu : "  << (s.IsGPU()?"YES":"NO") << std::endl
        << "               maxComputeUnits : "  << s.GetMaxComputeUnits() << std::endl
        << "               maxWorkItemSize : "  << s.GetMaxWorkItemSize() << std::endl
        << "              maxWorkGroupSize : "  << s.GetMaxWorkGroupSize() << std::endl
        << "             maxClockFrequency : "  << s.GetMaxClockFrequencyMhz() << "MHz" << std::endl
        << "                 globalMemSize : "  << static_cast<double>(s.GetGlobalMemSize())/1024./1024. << "MiB" << std::endl
        << "         maxConstantBufferSize : "  << static_cast<double>(s.GetMaxConstantBufferSize())/1024. << "kiB" << std::endl
        << "                  localMemSize : "  << static_cast<double>(s.GetLocalMemSize())/1024. << "kiB" << std::endl
        << "             dedicatedLocalMem : "  << (s.HasDedicatedLocalMem()?"YES":"NO") << std::endl
        << "        errorCorrectionSupport : "  << (s.HasErrorCorrectionSupport()?"YES":"NO") << std::endl
        << "                     available : "  << (s.IsAvailable()?"YES":"NO") << std::endl
        << "                        vendor : "  << s.GetVendor() << std::endl
        << "                 driverVersion : "  << s.GetDriverVersion() << std::endl
        << "                 deviceVersion : "  << s.GetDeviceVersion() << std::endl
        << "                    extensions : "  << s.GetExtensions() << std::endl
        << std::endl
        << "                 useNativeMath : "  << (s.GetUseNativeMath()?"YES":"NO") << std::endl
        << "  approximateNumberOfWorkItems : "  << s.GetApproximateNumberOfWorkItems() << std::endl
        << "]" ;

    return oss.str();
}

void register_I3CLSimOpenCLDevice()
{
    {
        bp::scope I3CLSimOpenCLDevice_scope = 
        bp::class_<I3CLSimOpenCLDevice, boost::shared_ptr<I3CLSimOpenCLDevice> >
        ("I3CLSimOpenCLDevice", 
         bp::init<const std::string &, const std::string &>
         (
          (
           bp::arg("platformName"),
           bp::arg("deviceName")
          )
         )
        )
        .def(init<const std::string &, const std::string &, bool, uint32_t>
             (
              (
               bp::arg("platformName"),
               bp::arg("deviceName"),
               bp::arg("useNativeMath"),
               bp::arg("approximateNumberOfWorkItems")
               )
              )
             ) 

        .def("__eq__", &I3CLSimOpenCLDevice_equalWrap)
        .def("__str__", &I3CLSimOpenCLDevice_prettyprint)
        .def("__repr__", &I3CLSimOpenCLDevice_repr)

        .def(copy_suite<I3CLSimOpenCLDevice>())

        .def("GetAllDevices", &I3CLSimOpenCLDevice::GetAllDevices)
        .staticmethod("GetAllDevices")

        .def("SplitDevice", &I3CLSimOpenCLDevice::SplitDevice)

        .def("GetPlatformName", &I3CLSimOpenCLDevice::GetPlatformName, bp::return_value_policy<bp::copy_const_reference>())
        .def("GetDeviceName", &I3CLSimOpenCLDevice::GetDeviceName, bp::return_value_policy<bp::copy_const_reference>())
        .add_property("platform", 
                      make_function(&I3CLSimOpenCLDevice::GetPlatformName, return_value_policy<copy_const_reference>()))
        .add_property("device", 
                      make_function(&I3CLSimOpenCLDevice::GetDeviceName, return_value_policy<copy_const_reference>()))

        
        .def("IsCPU", &I3CLSimOpenCLDevice::IsCPU)
        .add_property("cpu", &I3CLSimOpenCLDevice::IsCPU)
        .def("IsGPU", &I3CLSimOpenCLDevice::IsGPU)
        .add_property("gpu", &I3CLSimOpenCLDevice::IsGPU)
        .def("GetMaxComputeUnits", &I3CLSimOpenCLDevice::GetMaxComputeUnits)
        .add_property("maxComputeUnits", &I3CLSimOpenCLDevice::GetMaxComputeUnits)
        .def("GetMaxWorkItemSize", &I3CLSimOpenCLDevice::GetMaxWorkItemSize)
        .add_property("maxWorkItemSize", &I3CLSimOpenCLDevice::GetMaxWorkItemSize)
        .def("GetMaxWorkGroupSize", &I3CLSimOpenCLDevice::GetMaxWorkGroupSize)
        .add_property("maxWorkGroupSize", &I3CLSimOpenCLDevice::GetMaxWorkGroupSize)
        .def("GetMaxClockFrequencyMhz", &I3CLSimOpenCLDevice::GetMaxClockFrequencyMhz)
        .add_property("maxClockFrequencyMhz", &I3CLSimOpenCLDevice::GetMaxClockFrequencyMhz)
        .def("GetGlobalMemSize", &I3CLSimOpenCLDevice::GetGlobalMemSize)
        .add_property("globalMemSize", &I3CLSimOpenCLDevice::GetGlobalMemSize)
        .def("GetMaxConstantBufferSize", &I3CLSimOpenCLDevice::GetMaxConstantBufferSize)
        .add_property("maxConstantBufferSize", &I3CLSimOpenCLDevice::GetMaxConstantBufferSize)
        .def("GetLocalMemSize", &I3CLSimOpenCLDevice::GetLocalMemSize)
        .add_property("localMemSize", &I3CLSimOpenCLDevice::GetLocalMemSize)
        .def("HasDedicatedLocalMem", &I3CLSimOpenCLDevice::HasDedicatedLocalMem)
        .add_property("dedicatedLocalMem", &I3CLSimOpenCLDevice::HasDedicatedLocalMem)
        .def("HasErrorCorrectionSupport", &I3CLSimOpenCLDevice::HasErrorCorrectionSupport)
        .add_property("errorCorrectionSupport", &I3CLSimOpenCLDevice::HasErrorCorrectionSupport)
        .def("IsAvailable", &I3CLSimOpenCLDevice::IsAvailable)
        .add_property("available", &I3CLSimOpenCLDevice::IsAvailable)
        .def("GetVendor", &I3CLSimOpenCLDevice::GetVendor)
        .add_property("vendor", &I3CLSimOpenCLDevice::GetVendor)
        .def("GetDriverVersion", &I3CLSimOpenCLDevice::GetDriverVersion)
        .add_property("driverVersion", &I3CLSimOpenCLDevice::GetDriverVersion)
        .def("GetDeviceVersion", &I3CLSimOpenCLDevice::GetDeviceVersion)
        .add_property("deviceVersion", &I3CLSimOpenCLDevice::GetDeviceVersion)
        .def("GetExtensions", &I3CLSimOpenCLDevice::GetExtensions)
        .add_property("extensions", &I3CLSimOpenCLDevice::GetExtensions)

        
        .def("GetUseNativeMath", &I3CLSimOpenCLDevice::GetUseNativeMath)
        .def("SetUseNativeMath", &I3CLSimOpenCLDevice::SetUseNativeMath)
        .add_property("useNativeMath", &I3CLSimOpenCLDevice::GetExtensions, &I3CLSimOpenCLDevice::SetUseNativeMath)

        .def("GetApproximateNumberOfWorkItems", &I3CLSimOpenCLDevice::GetApproximateNumberOfWorkItems)
        .def("SetApproximateNumberOfWorkItems", &I3CLSimOpenCLDevice::SetApproximateNumberOfWorkItems)
        .add_property("approximateNumberOfWorkItems", &I3CLSimOpenCLDevice::GetApproximateNumberOfWorkItems, &I3CLSimOpenCLDevice::SetApproximateNumberOfWorkItems)
        ;
    }
    

    class_<I3CLSimOpenCLDeviceSeries, I3CLSimOpenCLDeviceSeriesPtr>("I3CLSimOpenCLDeviceSeries")
    .def(list_indexing_suite<I3CLSimOpenCLDeviceSeries>())
    ;
    
    bp::implicitly_convertible<shared_ptr<I3CLSimOpenCLDevice>, shared_ptr<const I3CLSimOpenCLDevice> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimOpenCLDeviceSeries>, shared_ptr<const I3CLSimOpenCLDeviceSeries> >();

    from_python_sequence<I3CLSimOpenCLDeviceSeries, variable_capacity_policy>();
    
}
