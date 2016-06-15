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
 * @file I3CLSimLightSourceToStepConverter.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/I3CLSimLightSourceToStepConverter.h>
#include <clsim/I3CLSimLightSourceToStepConverterGeant4.h>
#include <clsim/I3CLSimLightSourceToStepConverterPPC.h>
#include <clsim/I3CLSimLightSourceToStepConverterFlasher.h>

#include <boost/preprocessor/seq.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_const.hpp>

#include "python_gil_holder.h"

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimLightSourceToStepConverterWrapper : I3CLSimLightSourceToStepConverter, bp::wrapper<I3CLSimLightSourceToStepConverter>
{
    // pure virtual
    virtual void SetBunchSizeGranularity(uint64_t num) {utils::python_gil_holder gil; this->get_override("SetBunchSizeGranularity")(num);}
    virtual void SetMaxBunchSize(uint64_t num) {utils::python_gil_holder gil; this->get_override("SetMaxBunchSize")(num);}
    virtual void SetRandomService(I3RandomServicePtr random) {utils::python_gil_holder gil; this->get_override("SetRandomService")(random);}
    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias) {utils::python_gil_holder gil; this->get_override("SetWlenBias")(wlenBias);}
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) {utils::python_gil_holder gil; this->get_override("SetMediumProperties")(mediumProperties);}
    virtual void Initialize() {utils::python_gil_holder gil; this->get_override("Initialize")();}
    virtual bool IsInitialized() const {utils::python_gil_holder gil; return this->get_override("IsInitialized")();}
    virtual void EnqueueLightSource(const I3CLSimLightSource &lightSource, uint32_t identifier) {utils::python_gil_holder gil; this->get_override("EnqueueLightSource")(lightSource, identifier);}
    virtual void EnqueueBarrier() {utils::python_gil_holder gil; this->get_override("EnqueueLightSource")();}
    virtual bool BarrierActive() const {utils::python_gil_holder gil; return this->get_override("BarrierActive")();}
    virtual bool MoreStepsAvailable() const {utils::python_gil_holder gil; return this->get_override("MoreStepsAvailable")();}
    
    virtual I3CLSimStepSeriesConstPtr GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
    {utils::python_gil_holder gil; return this->get_override("GetConversionResultWithBarrierInfo")(barrierWasReset, timeout);}

    virtual I3CLSimStepSeriesConstPtr GetConversionResult(double timeout)
    {
        utils::python_gil_holder gil;
        if (override f = this->get_override("GetConversionResult")) {
            return f(timeout);
        } else {
            return I3CLSimLightSourceToStepConverter::GetConversionResult(timeout);
        }
    }
    I3CLSimStepSeriesConstPtr default_GetConversionResult(double timeout) {utils::python_gil_holder gil; return this->get_override("GetConversionResult")(timeout);}

    virtual void SetLightSourceParameterizationSeries(const I3CLSimLightSourceParameterizationSeries &parameterizationSeries_)
    {
        utils::python_gil_holder gil;
        if (override f = this->get_override("SetLightSourceParameterizationSeries")) {
            f(parameterizationSeries_);
        } else {
            I3CLSimLightSourceToStepConverter::SetLightSourceParameterizationSeries(parameterizationSeries_);
        }
    }
    void default_SetLightSourceParameterizationSeries(const I3CLSimLightSourceParameterizationSeries &parameterizationSeries_) 
    {this->I3CLSimLightSourceToStepConverter::SetLightSourceParameterizationSeries(parameterizationSeries_);}


    virtual const I3CLSimLightSourceParameterizationSeries &GetLightSourceParameterizationSeries() const
    {
        utils::python_gil_holder gil;
        if (override f = this->get_override("GetLightSourceParameterizationSeries")) {
            return f();
        } else {
            return I3CLSimLightSourceToStepConverter::GetLightSourceParameterizationSeries();
        }
    }
    const I3CLSimLightSourceParameterizationSeries &default_GetLightSourceParameterizationSeries() const
    {return this->I3CLSimLightSourceToStepConverter::GetLightSourceParameterizationSeries();}

};

/*
//
// visitor for converting contents of variant to object
//
struct vc : boost::static_visitor<bp::object>
{
    template <typename T>
    bp::object operator()(const boost::shared_ptr<const T>& v) const
    {
        return bp::object(boost::const_pointer_cast<T>(v));
    }

    // this is done for non-const shared pointers
    template <typename T>
    bp::object operator()(const boost::shared_ptr<T>& v) const
    {
        return bp::object(v);
    }
    
    // do this for everything but boost::shared_ptrs
    template <typename T>
    bp::object operator()(const T& v,
                          typename boost::disable_if<is_shared_ptr<T> >::type * = 0) const
    {
        return bp::object(v);
    }

    bp::object operator()(const std::pair<uint32_t, I3ParticleConstPtr>& v) const
    {
        return bp::make_tuple(v.first, boost::const_pointer_cast<I3Particle>(v.second));
    }
    
};

/// to-python convert to I3CLSimLightSourceToStepConverter::ConversionResult_t
struct I3CLSimLightSourceToStepConverter_ConversionResult_t_to_python
{
    static PyObject *convert(const I3CLSimLightSourceToStepConverter::ConversionResult_t& val)
    {
        return boost::python::incref(boost::apply_visitor(vc(), val).ptr()); 
    }
    
};
*/

void register_I3CLSimLightSourceToStepConverter()
{
    {
        bp::scope I3CLSimLightSourceToStepConverter_scope = 
        bp::class_<I3CLSimLightSourceToStepConverterWrapper, boost::shared_ptr<I3CLSimLightSourceToStepConverterWrapper>, boost::noncopyable>("I3CLSimLightSourceToStepConverter", bp::no_init)
        .def("SetBunchSizeGranularity", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::SetBunchSizeGranularity))
        .def("SetMaxBunchSize", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::SetMaxBunchSize))
        .def("SetRandomService", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::SetRandomService))
        .def("SetWlenBias", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::SetWlenBias))
        .def("SetMediumProperties", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::SetMediumProperties))
        .def("Initialize", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::Initialize))
        .def("IsInitialized", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::IsInitialized))
        .def("EnqueueLightSource", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::EnqueueLightSource))
        .def("EnqueueBarrier", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::EnqueueBarrier))
        .def("BarrierActive", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::BarrierActive))
        .def("MoreStepsAvailable", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::MoreStepsAvailable))
        .def("GetConversionResult", bp::pure_virtual(&I3CLSimLightSourceToStepConverter::GetConversionResult), bp::arg("timeout")=NAN)
        .def("GetConversionResultWithBarrierInfo", 
             &I3CLSimLightSourceToStepConverter::GetConversionResultWithBarrierInfo,
             &I3CLSimLightSourceToStepConverterWrapper::GetConversionResultWithBarrierInfo,
             bp::arg("barrierWasReset"), bp::arg("timeout")=NAN)

        .def("SetLightSourceParameterizationSeries", 
             &I3CLSimLightSourceToStepConverter::SetLightSourceParameterizationSeries,
             &I3CLSimLightSourceToStepConverterWrapper::default_SetLightSourceParameterizationSeries)

        .def("GetLightSourceParameterizationSeries", 
             &I3CLSimLightSourceToStepConverter::GetLightSourceParameterizationSeries,
             &I3CLSimLightSourceToStepConverterWrapper::default_GetLightSourceParameterizationSeries,
             bp::return_value_policy<bp::copy_const_reference>())
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterWrapper>, boost::shared_ptr<const I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterWrapper>, boost::shared_ptr<I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterWrapper>, boost::shared_ptr<const I3CLSimLightSourceToStepConverterWrapper> >();
    
    //bp::to_python_converter<I3CLSimLightSourceToStepConverter::ConversionResult_t, I3CLSimLightSourceToStepConverter_ConversionResult_t_to_python>();
    
    // I3CLSimLightSourceToStepConverterGeant4
    {
        bp::class_<
        I3CLSimLightSourceToStepConverterGeant4, 
        boost::shared_ptr<I3CLSimLightSourceToStepConverterGeant4>, 
        bases<I3CLSimLightSourceToStepConverter>,
        boost::noncopyable
        >
        (
         "I3CLSimLightSourceToStepConverterGeant4",
         bp::init<
         std::string,
         double,
         uint32_t,
         uint32_t
         >(
           (
            bp::arg("physicsListName") = I3CLSimLightSourceToStepConverterGeant4::default_physicsListName,
            bp::arg("maxBetaChangePerStep") = I3CLSimLightSourceToStepConverterGeant4::default_maxBetaChangePerStep,
            bp::arg("maxNumPhotonsPerStep") = I3CLSimLightSourceToStepConverterGeant4::default_maxNumPhotonsPerStep,
            bp::arg("maxQueueItems") = I3CLSimLightSourceToStepConverterGeant4::default_maxQueueItems
           )
          )
        )
        .add_static_property("can_use_geant4",bp::make_getter(I3CLSimLightSourceToStepConverterGeant4::canUseGeant4))
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterGeant4>, boost::shared_ptr<const I3CLSimLightSourceToStepConverterGeant4> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterGeant4>, boost::shared_ptr<I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterGeant4>, boost::shared_ptr<const I3CLSimLightSourceToStepConverter> >();
    
    
    // I3CLSimLightSourceToStepConverterPPC
    {
        bp::class_<
        I3CLSimLightSourceToStepConverterPPC, 
        boost::shared_ptr<I3CLSimLightSourceToStepConverterPPC>, 
        bases<I3CLSimLightSourceToStepConverter>,
        boost::noncopyable
        >
        (
         "I3CLSimLightSourceToStepConverterPPC",
         bp::init<
         uint32_t, uint32_t, double
         >(
           (
            bp::arg("photonsPerStep") = I3CLSimLightSourceToStepConverterPPC::default_photonsPerStep,
            bp::arg("highPhotonsPerStep") = I3CLSimLightSourceToStepConverterPPC::default_highPhotonsPerStep,
            bp::arg("useHighPhotonsPerStepStartingFromNumPhotons") = I3CLSimLightSourceToStepConverterPPC::default_useHighPhotonsPerStepStartingFromNumPhotons
            )
           )
         )
        .def("SetUseCascadeExtension", &I3CLSimLightSourceToStepConverterPPC::SetUseCascadeExtension)
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterPPC>, boost::shared_ptr<const I3CLSimLightSourceToStepConverterPPC> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterPPC>, boost::shared_ptr<I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterPPC>, boost::shared_ptr<const I3CLSimLightSourceToStepConverter> >();

    
    // I3CLSimLightSourceToStepConverterFlasher
    {
        bp::class_<
        I3CLSimLightSourceToStepConverterFlasher, 
        boost::shared_ptr<I3CLSimLightSourceToStepConverterFlasher>,
        bases<I3CLSimLightSourceToStepConverter>,
        boost::noncopyable
        >
        (
         "I3CLSimLightSourceToStepConverterFlasher",
         bp::init<
         I3CLSimFunctionConstPtr,
         I3CLSimSpectrumTablePtr,
         I3CLSimRandomValueConstPtr,
         I3CLSimRandomValueConstPtr,
         I3CLSimRandomValueConstPtr,
         bool,
         uint32_t
         >(
           (
            bp::arg("flasherSpectrumNoBias"),
            bp::arg("spectrumTable"),
            bp::arg("angularProfileDistributionPolar"),
            bp::arg("angularProfileDistributionAzimuthal"),
            bp::arg("timeDelayDistribution"),
            bp::arg("interpretAngularDistributionsInPolarCoordinates") = I3CLSimLightSourceToStepConverterFlasher::default_interpretAngularDistributionsInPolarCoordinates,
            bp::arg("photonsPerStep") = I3CLSimLightSourceToStepConverterFlasher::default_photonsPerStep
            )
           )
         )
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterFlasher>, boost::shared_ptr<const I3CLSimLightSourceToStepConverterFlasher> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterFlasher>, boost::shared_ptr<I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourceToStepConverterFlasher>, boost::shared_ptr<const I3CLSimLightSourceToStepConverter> >();

}
