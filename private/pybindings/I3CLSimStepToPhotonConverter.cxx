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

#include <clsim/I3CLSimServer.h>
#include <clsim/I3CLSimStepToPhotonConverter.h>
#include <clsim/I3CLSimStepToPhotonConverterOpenCL.h>

#include <boost/preprocessor/seq.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_const.hpp>

#include <boost/foreach.hpp>

#include "icetray/python/list_indexing_suite.hpp"

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
    virtual std::size_t GetWorkgroupSize() const { utils::python_gil_holder gil; return this->get_override("GetWorkgroupSize")(); };
    virtual std::size_t GetMaxNumWorkitems() const { utils::python_gil_holder gil; return this->get_override("GetMaxNumWorkitems")(); };
    virtual std::size_t QueueSize() const {utils::python_gil_holder gil; return this->get_override("QueueSize")();}
    virtual bool MorePhotonsAvailable() const {utils::python_gil_holder gil; return this->get_override("MorePhotonsAvailable")();}
    virtual I3CLSimStepToPhotonConverter::ConversionResult_t GetConversionResult() {
        utils::python_gil_holder gil;
        I3CLSimStepToPhotonConverter::ConversionResult_t result = this->get_override("GetConversionResult")();
        // Make a deep copy so that the return value can be safely destroyed without holding the GIL
        if (result.photons)
            result.photons = boost::make_shared<I3CLSimPhotonSeries>(*result.photons);
        if (result.photonHistories)
            result.photonHistories = boost::make_shared<I3CLSimPhotonHistorySeries>(*result.photonHistories);
        
        return result;
    }
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
        ("I3CLSimStepToPhotonConverter")
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
        .add_property("workgroupSize", &I3CLSimStepToPhotonConverter::GetWorkgroupSize)
        .add_property("maxNumWorkitems", &I3CLSimStepToPhotonConverter::GetMaxNumWorkitems)
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
         "I3CLSimStepToPhotonConverterOpenCLBase",
         bp::init<
         I3RandomServicePtr,bool
         >(
           (
            bp::arg("RandomService"),
            bp::arg("UseNativeMath")=I3CLSimStepToPhotonConverterOpenCLWrapper::default_useNativeMath
           )
          )
        )
        
        .def("GetGeometrySource", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetGeometrySource)
        .def("GetCollisionDetectionSource", &I3CLSimStepToPhotonConverterOpenCLWrapper::GetCollisionDetectionSource)
        
        ;
    }
    
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
         I3RandomServicePtr,bool
         >(
           (
            bp::arg("RandomService"),
            bp::arg("UseNativeMath")=I3CLSimStepToPhotonConverterOpenCL::default_useNativeMath
           )
          )
        )
        .def("Compile", &I3CLSimStepToPhotonConverterOpenCL::Compile)
        .def("GetFullSource", &I3CLSimStepToPhotonConverterOpenCL::GetFullSource)
                
        .def("GetGeometrySource", &I3CLSimStepToPhotonConverterOpenCL::GetGeometrySource)
        .def("GetCollisionDetectionSource", &I3CLSimStepToPhotonConverterOpenCL::GetCollisionDetectionSource)
        
        .def("SetDevice", &I3CLSimStepToPhotonConverterOpenCL::SetDevice)
        .def("GetMaxWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCL::GetMaxWorkgroupSize)
        .add_property("maxWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCL::GetMaxWorkgroupSize)

        .def("GetWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCL::GetWorkgroupSize)
        .def("SetWorkgroupSize", &I3CLSimStepToPhotonConverterOpenCL::SetWorkgroupSize)
        .def("GetMaxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCL::GetMaxNumWorkitems)
        .def("SetMaxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCL::SetMaxNumWorkitems)

        .def("SetEnableDoubleBuffering", &I3CLSimStepToPhotonConverterOpenCL::SetEnableDoubleBuffering)
        .def("GetEnableDoubleBuffering", &I3CLSimStepToPhotonConverterOpenCL::GetEnableDoubleBuffering)

        .def("SetDoublePrecision", &I3CLSimStepToPhotonConverterOpenCL::SetDoublePrecision)
        .def("GetDoublePrecision", &I3CLSimStepToPhotonConverterOpenCL::GetDoublePrecision)

        .def("SetStopDetectedPhotons", &I3CLSimStepToPhotonConverterOpenCL::SetStopDetectedPhotons)
        .def("GetStopDetectedPhotons", &I3CLSimStepToPhotonConverterOpenCL::GetStopDetectedPhotons)

        .def("SetSaveAllPhotons", &I3CLSimStepToPhotonConverterOpenCL::SetSaveAllPhotons)
        .def("GetSaveAllPhotons", &I3CLSimStepToPhotonConverterOpenCL::GetSaveAllPhotons)

        .def("SetSaveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterOpenCL::SetSaveAllPhotonsPrescale)
        .def("GetSaveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterOpenCL::GetSaveAllPhotonsPrescale)

        .def("SetPhotonHistoryEntries", &I3CLSimStepToPhotonConverterOpenCL::SetPhotonHistoryEntries)
        .def("GetPhotonHistoryEntries", &I3CLSimStepToPhotonConverterOpenCL::GetPhotonHistoryEntries)

        .def("SetFixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterOpenCL::SetFixedNumberOfAbsorptionLengths)
        .def("GetFixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterOpenCL::GetFixedNumberOfAbsorptionLengths)

        .def("SetDOMPancakeFactor", &I3CLSimStepToPhotonConverterOpenCL::SetDOMPancakeFactor)
        .def("GetDOMPancakeFactor", &I3CLSimStepToPhotonConverterOpenCL::GetDOMPancakeFactor)

        
        .add_property("workgroupSize", &I3CLSimStepToPhotonConverterOpenCL::GetWorkgroupSize, &I3CLSimStepToPhotonConverterOpenCL::SetWorkgroupSize)
        .add_property("maxNumWorkitems", &I3CLSimStepToPhotonConverterOpenCL::GetMaxNumWorkitems, &I3CLSimStepToPhotonConverterOpenCL::SetMaxNumWorkitems)
        .add_property("enableDoubleBuffering", &I3CLSimStepToPhotonConverterOpenCL::GetEnableDoubleBuffering, &I3CLSimStepToPhotonConverterOpenCL::SetEnableDoubleBuffering)
        .add_property("doublePrecision", &I3CLSimStepToPhotonConverterOpenCL::GetDoublePrecision, &I3CLSimStepToPhotonConverterOpenCL::SetDoublePrecision)
        .add_property("stopDetectedPhotons", &I3CLSimStepToPhotonConverterOpenCL::GetStopDetectedPhotons, &I3CLSimStepToPhotonConverterOpenCL::SetStopDetectedPhotons)
        .add_property("saveAllPhotons", &I3CLSimStepToPhotonConverterOpenCL::GetSaveAllPhotons, &I3CLSimStepToPhotonConverterOpenCL::SetSaveAllPhotons)
        .add_property("saveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterOpenCL::GetSaveAllPhotonsPrescale, &I3CLSimStepToPhotonConverterOpenCL::SetSaveAllPhotonsPrescale)
        .add_property("photonHistoryEntries", &I3CLSimStepToPhotonConverterOpenCL::GetPhotonHistoryEntries, &I3CLSimStepToPhotonConverterOpenCL::SetPhotonHistoryEntries)
        .add_property("fixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterOpenCL::GetFixedNumberOfAbsorptionLengths, &I3CLSimStepToPhotonConverterOpenCL::SetFixedNumberOfAbsorptionLengths)
        .add_property("DOMPancakeFactor", &I3CLSimStepToPhotonConverterOpenCL::GetDOMPancakeFactor, &I3CLSimStepToPhotonConverterOpenCL::SetDOMPancakeFactor)
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterOpenCLWrapper>, boost::shared_ptr<const I3CLSimStepToPhotonConverterOpenCLWrapper> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterOpenCLWrapper>, boost::shared_ptr<I3CLSimStepToPhotonConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterOpenCLWrapper>, boost::shared_ptr<const I3CLSimStepToPhotonConverter> >();
    
}

namespace {

bp::dict GetStatistics(const I3CLSimServer &server)
{
    bp::dict summary;
    
    for (auto &pair : server.GetStatistics())
        summary[pair.first] = pair.second;
    
    return summary;
}

}

void register_I3CLSimServer()
{
    typedef std::vector<I3CLSimStepToPhotonConverterPtr> ConverterSeries;
    bp::class_<ConverterSeries, boost::shared_ptr<ConverterSeries> >("I3CLSimStepToPhotonConverterSeries")
        .def(bp::list_indexing_suite<ConverterSeries>())
    ;
    
    bp::class_<I3CLSimServer, boost::shared_ptr<I3CLSimServer>, boost::noncopyable>("I3CLSimServer", bp::init<const std::string&, const ConverterSeries&>())
        .def("GetStatistics", &GetStatistics)
    ;
    
    bp::class_<I3CLSimClient, boost::shared_ptr<I3CLSimClient>, boost::noncopyable>("I3CLSimClient", bp::init<const std::string&>())
        .def("EnqueueSteps", &I3CLSimClient::EnqueueSteps)
        .def("GetConversionResult", &I3CLSimClient::GetConversionResult)
        .add_property("workgroupSize", &I3CLSimClient::GetWorkgroupSize)
        .add_property("maxNumWorkitems", &I3CLSimClient::GetMaxNumWorkitems)
    ;
}
