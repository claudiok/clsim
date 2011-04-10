//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module g4sim-interface.
//
//   g4sim-intrface is free software; you can redistribute it and/or modify
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

#include <sstream>

#include <clsim/I3CLSimStepToPhotonConverter.h>
#include <clsim/I3CLSimStepToPhotonConverterOpenCL.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_const.hpp>

#include <boost/foreach.hpp>


using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimStepToPhotonConverterWrapper : I3CLSimStepToPhotonConverter, bp::wrapper<I3CLSimStepToPhotonConverter>
{
    // pure virtual
    virtual void SetWlenGenerator(I3CLSimRandomValueConstPtr wlenGenerator) {this->get_override("SetWlenGenerator")(wlenGenerator);}
    virtual void SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias) {this->get_override("SetWlenBias")(wlenBias);}

    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) {this->get_override("SetMediumProperties")(mediumProperties);}
    virtual void SetGeometry(I3CLSimSimpleGeometryConstPtr geometry) {this->get_override("SetGeometry")(geometry);}

    virtual void Initialize() {this->get_override("Initialize")();}
    virtual bool IsInitialized() const {return this->get_override("IsInitialized")();}

    virtual void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier) {this->get_override("EnqueueSteps")(steps, identifier);}
    virtual bool MorePhotonsAvailable() const {return this->get_override("MorePhotonsAvailable")();}
    virtual I3CLSimStepToPhotonConverter::ConversionResult_t GetConversionResult() {return this->get_override("GetConversionResult")();}
};

namespace I3CLSimStepToPhotonConverter_utils
{
    template <typename T1, typename T2> 
    struct std_pair_to_tuple 
    { 
        static PyObject* convert(std::pair<T1, T2> const& p) 
        { 
            return boost::python::incref(boost::python::make_tuple(p.first, p.second).ptr()); 
        } 
    }; 

    template <typename T1, typename T2> 
    struct std_pair_to_python_converter 
    { 
        std_pair_to_python_converter() 
        { 
            boost::python::to_python_converter< 
            std::pair<T1, T2>, 
            std_pair_to_tuple<T1, T2> >(); 
        } 
    };
};

using namespace I3CLSimStepToPhotonConverter_utils;

bp::object I3CLSimStepToPhotonConverterOpenCL_GetDeviceList_python(I3CLSimStepToPhotonConverterOpenCL &ref)
{
    shared_ptr<const std::vector<std::pair<std::string, std::string> > > ret = ref.GetDeviceList();
    if (!ret) return bp::object();

    const std::vector<std::pair<std::string, std::string> > &vect = *ret;
    
    bp::list retList = bp::list();

    for (std::size_t i = 0; i < vect.size(); ++i)
    {
        retList.append(bp::make_tuple(vect[i].first, vect[i].second));
    }
    
    return retList;
}

void register_I3CLSimStepToPhotonConverter()
{
    {
        bp::scope I3CLSimStepToPhotonConverter_scope = 
        bp::class_<I3CLSimStepToPhotonConverterWrapper, boost::shared_ptr<I3CLSimStepToPhotonConverterWrapper>, boost::noncopyable>
        ("I3CLSimStepToPhotonConverter", bp::no_init)
        .def("SetWlenGenerator", bp::pure_virtual(&I3CLSimStepToPhotonConverter::SetWlenGenerator))
        .def("SetWlenBias", bp::pure_virtual(&I3CLSimStepToPhotonConverter::SetWlenBias))
        .def("SetMediumProperties", bp::pure_virtual(&I3CLSimStepToPhotonConverter::SetMediumProperties))
        .def("SetGeometry", bp::pure_virtual(&I3CLSimStepToPhotonConverter::SetGeometry))
        .def("Initialize", bp::pure_virtual(&I3CLSimStepToPhotonConverter::Initialize))
        .def("IsInitialized", bp::pure_virtual(&I3CLSimStepToPhotonConverter::IsInitialized))
        .def("EnqueueSteps", bp::pure_virtual(&I3CLSimStepToPhotonConverter::EnqueueSteps))
        .def("MorePhotonsAvailable", bp::pure_virtual(&I3CLSimStepToPhotonConverter::MorePhotonsAvailable))
        .def("GetConversionResult", bp::pure_virtual(&I3CLSimStepToPhotonConverter::GetConversionResult))
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimStepToPhotonConverterWrapper>, shared_ptr<const I3CLSimStepToPhotonConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimStepToPhotonConverterWrapper>, shared_ptr<I3CLSimStepToPhotonConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimStepToPhotonConverterWrapper>, shared_ptr<const I3CLSimStepToPhotonConverterWrapper> >();
    std_pair_to_python_converter<uint32_t, I3CLSimPhotonSeriesPtr>();
    
    // I3CLSimStepToPhotonConverterOpenCL
    {
        bp::class_<
        I3CLSimStepToPhotonConverterOpenCL, 
        boost::shared_ptr<I3CLSimStepToPhotonConverterOpenCL>, 
        bases<I3CLSimStepToPhotonConverter>,
        boost::noncopyable
        >
        (
         "I3CLSimStepToPhotonConverterOpenCL",
         bp::init<
         I3RandomServicePtr,bool,bool,bool
         >(
           (
            bp::arg("RandomService"),
            bp::arg("UseNativeMath")=I3CLSimStepToPhotonConverterOpenCL::default_useNativeMath,
            bp::arg("CPUOnly")=I3CLSimStepToPhotonConverterOpenCL::default_cpuOnly,
            bp::arg("GPUOnly")=I3CLSimStepToPhotonConverterOpenCL::default_gpuOnly
           )
          )
        )
        .def("Compile", &I3CLSimStepToPhotonConverterOpenCL::Compile)
        .def("GetFullSource", &I3CLSimStepToPhotonConverterOpenCL::GetFullSource)
        
        .def("GetDeviceList", &I3CLSimStepToPhotonConverterOpenCL_GetDeviceList_python)
        .def("SetDeviceIndex", &I3CLSimStepToPhotonConverterOpenCL::SetDeviceIndex)
        .def("SetDeviceName", &I3CLSimStepToPhotonConverterOpenCL::SetDeviceName)
        .def("GetMaxWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCL::GetMaxWorkgroupSize)
        .add_property("maxWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCL::GetMaxWorkgroupSize)

        .def("GetWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCL::GetWorkgroupSize)
        .def("SetWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCL::SetWorkgroupSize)
        .def("GetMaxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCL::GetMaxNumWorkitems)
        .def("SetMaxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCL::SetMaxNumWorkitems)

        .add_property("workgroupSize", &I3CLSimStepToPhotonConverterOpenCL::GetWorkgroupSize, &I3CLSimStepToPhotonConverterOpenCL::SetWorkgroupSize)
        .add_property("maxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCL::GetMaxNumWorkitems, &I3CLSimStepToPhotonConverterOpenCL::SetMaxNumWorkitems)
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimStepToPhotonConverterOpenCL>, shared_ptr<const I3CLSimStepToPhotonConverterOpenCL> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimStepToPhotonConverterOpenCL>, shared_ptr<I3CLSimStepToPhotonConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimStepToPhotonConverterOpenCL>, shared_ptr<const I3CLSimStepToPhotonConverter> >();
    
}
