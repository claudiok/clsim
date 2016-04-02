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
 * @file I3CLSimStepToPhotonConverter.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/I3CLSimStepToPhotonConverter.h>
#include <clsim/I3CLSimStepToPhotonConverterOpenCL.h>

#include <boost/preprocessor/seq.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_const.hpp>

#include <boost/foreach.hpp>

#include "python_gil_holder.h"

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimStepToPhotonConverterWrapper : I3CLSimStepToPhotonConverter, bp::wrapper<I3CLSimStepToPhotonConverter>
{
    // pure virtual
    virtual void SetWlenGenerators(const std::vector<I3CLSimRandomValueConstPtr> &wlenGenerators) {utils::python_gil_holder gil; this->get_override("SetWlenGenerators")(wlenGenerators);}
    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias) {utils::python_gil_holder gil; this->get_override("SetWlenBias")(wlenBias);}

    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) {utils::python_gil_holder gil; this->get_override("SetMediumProperties")(mediumProperties);}
    virtual void SetGeometry(I3CLSimSimpleGeometryConstPtr geometry) {utils::python_gil_holder gil; this->get_override("SetGeometry")(geometry);}

    virtual void Initialize() {utils::python_gil_holder gil; this->get_override("Initialize")();}
    virtual bool IsInitialized() const {utils::python_gil_holder gil; return this->get_override("IsInitialized")();}

    virtual void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier) {utils::python_gil_holder gil; this->get_override("EnqueueSteps")(steps, identifier);}
    virtual std::size_t QueueSize() const {utils::python_gil_holder gil; return this->get_override("QueueSize")();}
    virtual bool MorePhotonsAvailable() const {utils::python_gil_holder gil; return this->get_override("MorePhotonsAvailable")();}
    virtual I3CLSimStepToPhotonConverter::ConversionResult_t GetConversionResult() {utils::python_gil_holder gil; return this->get_override("GetConversionResult")();}
};

struct I3CLSimStepToPhotonConverterOpenCLWrapper : I3CLSimStepToPhotonConverterOpenCL, bp::wrapper<I3CLSimStepToPhotonConverterOpenCL> {
    I3CLSimStepToPhotonConverterOpenCLWrapper(I3RandomServicePtr rng, bool nm)
        : I3CLSimStepToPhotonConverterOpenCL(rng, nm) {}
    virtual std::string GetGeometrySource() {
        utils::python_gil_holder gil;
        if (this->get_override("GetGeometrySource"))
            return this->get_override("GetGeometrySource")();
        else
            return I3CLSimStepToPhotonConverterOpenCL::GetGeometrySource();
    }
    virtual std::string GetCollisionDetectionSource(bool header) {
        utils::python_gil_holder gil;
        if (this->get_override("GetCollisionDetectionSource"))
            return this->get_override("GetCollisionDetectionSource")(header);
        else
            return I3CLSimStepToPhotonConverterOpenCL::GetCollisionDetectionSource(header);
    }
};


void register_I3CLSimStepToPhotonConverter()
{
    {
        bp::scope I3CLSimStepToPhotonConverter_scope = 
        bp::class_<I3CLSimStepToPhotonConverterWrapper, boost::shared_ptr<I3CLSimStepToPhotonConverterWrapper>, boost::noncopyable>
        ("I3CLSimStepToPhotonConverter", bp::no_init)
        .def("SetWlenGenerators", bp::pure_virtual(&I3CLSimStepToPhotonConverter::SetWlenGenerators))
        .def("SetWlenBias", bp::pure_virtual(&I3CLSimStepToPhotonConverter::SetWlenBias))
        .def("SetMediumProperties", bp::pure_virtual(&I3CLSimStepToPhotonConverter::SetMediumProperties))
        .def("SetGeometry", bp::pure_virtual(&I3CLSimStepToPhotonConverter::SetGeometry))
        .def("Initialize", bp::pure_virtual(&I3CLSimStepToPhotonConverter::Initialize))
        .def("IsInitialized", bp::pure_virtual(&I3CLSimStepToPhotonConverter::IsInitialized))
        .def("EnqueueSteps", bp::pure_virtual(&I3CLSimStepToPhotonConverter::EnqueueSteps))
        .def("QueueSize", bp::pure_virtual(&I3CLSimStepToPhotonConverter::QueueSize))
        .def("MorePhotonsAvailable", bp::pure_virtual(&I3CLSimStepToPhotonConverter::MorePhotonsAvailable))
        .def("GetConversionResult", bp::pure_virtual(&I3CLSimStepToPhotonConverter::GetConversionResult))
        ;
        
        
        bp::class_<I3CLSimStepToPhotonConverter::ConversionResult_t>
        ("ConversionResult_t", 
         bp::init<uint32_t, I3CLSimPhotonSeriesPtr, I3CLSimPhotonHistorySeriesPtr>
         (
          (
           bp::arg("identifier"),
           bp::arg("photons")=I3CLSimPhotonSeriesPtr(),
           bp::arg("photonHistories")=I3CLSimPhotonHistorySeriesPtr()
          )
         )
        )
        .def(bp::init<>()) // add a default constructor, too
        .def_readwrite("identifier", &I3CLSimStepToPhotonConverter::ConversionResult_t::identifier)
        .def_readwrite("photons", &I3CLSimStepToPhotonConverter::ConversionResult_t::photons)
        .def_readwrite("photonHistories", &I3CLSimStepToPhotonConverter::ConversionResult_t::photonHistories)
        ;
        
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterWrapper>, boost::shared_ptr<const I3CLSimStepToPhotonConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterWrapper>, boost::shared_ptr<I3CLSimStepToPhotonConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterWrapper>, boost::shared_ptr<const I3CLSimStepToPhotonConverterWrapper> >();
    

    
    // I3CLSimStepToPhotonConverterOpenCL
    {

        
        bp::class_<
        I3CLSimStepToPhotonConverterOpenCLWrapper, 
        boost::shared_ptr<I3CLSimStepToPhotonConverterOpenCLWrapper>, 
        bases<I3CLSimStepToPhotonConverter>,
        boost::noncopyable
        >
        (
         "I3CLSimStepToPhotonConverterOpenCL",
         bp::init<
         I3RandomServicePtr,bool
         >(
           (
            bp::arg("RandomService"),
            bp::arg("UseNativeMath")=I3CLSimStepToPhotonConverterOpenCLWrapper::default_useNativeMath
           )
          )
        )
        .def("Compile", &I3CLSimStepToPhotonConverterOpenCLWrapper::Compile)
        .def("GetFullSource", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetFullSource)
                
        .def("GetGeometrySource", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetGeometrySource)
        .def("GetCollisionDetectionSource", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetCollisionDetectionSource)
        
        .def("SetDevice", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetDevice)
        .def("GetMaxWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetMaxWorkgroupSize)
        .add_property("maxWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetMaxWorkgroupSize)

        .def("GetWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetWorkgroupSize)
        .def("SetWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetWorkgroupSize)
        .def("GetMaxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetMaxNumWorkitems)
        .def("SetMaxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetMaxNumWorkitems)

        .def("SetEnableDoubleBuffering", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetEnableDoubleBuffering)
        .def("GetEnableDoubleBuffering", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetEnableDoubleBuffering)

        .def("SetDoublePrecision", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetDoublePrecision)
        .def("GetDoublePrecision", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetDoublePrecision)

        .def("SetStopDetectedPhotons", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetStopDetectedPhotons)
        .def("GetStopDetectedPhotons", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetStopDetectedPhotons)

        .def("SetSaveAllPhotons", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetSaveAllPhotons)
        .def("GetSaveAllPhotons", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetSaveAllPhotons)

        .def("SetSaveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetSaveAllPhotonsPrescale)
        .def("GetSaveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetSaveAllPhotonsPrescale)

        .def("SetPhotonHistoryEntries", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetPhotonHistoryEntries)
        .def("GetPhotonHistoryEntries", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetPhotonHistoryEntries)

        .def("SetFixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetFixedNumberOfAbsorptionLengths)
        .def("GetFixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetFixedNumberOfAbsorptionLengths)

        .def("SetDOMPancakeFactor", &I3CLSimStepToPhotonConverterOpenCLWrapper::SetDOMPancakeFactor)
        .def("GetDOMPancakeFactor", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetDOMPancakeFactor)

        
        .add_property("workgroupSize", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetWorkgroupSize, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetWorkgroupSize)
        .add_property("maxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetMaxNumWorkitems, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetMaxNumWorkitems)
        .add_property("enableDoubleBuffering", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetEnableDoubleBuffering, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetEnableDoubleBuffering)
        .add_property("doublePrecision", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetDoublePrecision, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetDoublePrecision)
        .add_property("stopDetectedPhotons", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetStopDetectedPhotons, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetStopDetectedPhotons)
        .add_property("saveAllPhotons", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetSaveAllPhotons, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetSaveAllPhotons)
        .add_property("saveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetSaveAllPhotonsPrescale, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetSaveAllPhotonsPrescale)
        .add_property("photonHistoryEntries", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetPhotonHistoryEntries, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetPhotonHistoryEntries)
        .add_property("fixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetFixedNumberOfAbsorptionLengths, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetFixedNumberOfAbsorptionLengths)
        .add_property("DOMPancakeFactor", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetDOMPancakeFactor, &I3CLSimStepToPhotonConverterOpenCLWrapper::SetDOMPancakeFactor)
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterOpenCLWrapper>, boost::shared_ptr<const I3CLSimStepToPhotonConverterOpenCLWrapper> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterOpenCLWrapper>, boost::shared_ptr<I3CLSimStepToPhotonConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterOpenCLWrapper>, boost::shared_ptr<const I3CLSimStepToPhotonConverter> >();
    
}
