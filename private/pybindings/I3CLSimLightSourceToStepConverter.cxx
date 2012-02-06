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

#include <clsim/I3CLSimLightSourceToStepConverter.h>
#include <clsim/I3CLSimLightSourceToStepConverterGeant4.h>
#include <clsim/I3CLSimLightSourceToStepConverterPPC.h>
#include <clsim/I3CLSimLightSourceToStepConverterFlasher.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_const.hpp>

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimLightSourceToStepConverterWrapper : I3CLSimLightSourceToStepConverter, bp::wrapper<I3CLSimLightSourceToStepConverter>
{
    // pure virtual
    virtual void SetBunchSizeGranularity(uint64_t num) {this->get_override("SetBunchSizeGranularity")(num);}
    virtual void SetMaxBunchSize(uint64_t num) {this->get_override("SetMaxBunchSize")(num);}
    virtual void SetRandomService(I3RandomServicePtr random) {this->get_override("SetRandomService")(random);}
    virtual void SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias) {this->get_override("SetWlenBias")(wlenBias);}
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) {this->get_override("SetMediumProperties")(mediumProperties);}
    virtual void Initialize() {this->get_override("Initialize")();}
    virtual bool IsInitialized() const {return this->get_override("IsInitialized")();}
    virtual void EnqueueLightSource(const I3Particle &particle, uint32_t identifier) {this->get_override("EnqueueLightSource")(particle, identifier);}
    virtual void EnqueueBarrier() {this->get_override("EnqueueLightSource")();}
    virtual bool BarrierActive() const {return this->get_override("BarrierActive")();}
    virtual bool MoreStepsAvailable() const {return this->get_override("MoreStepsAvailable")();}
    
    virtual I3CLSimStepSeriesConstPtr GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
    {return this->get_override("GetConversionResultWithBarrierInfo")(barrierWasReset, timeout);}

    virtual I3CLSimStepSeriesConstPtr GetConversionResult(double timeout)
    {
        if (override f = this->get_override("GetConversionResult")) {
            return f(timeout);
        } else {
            return I3CLSimLightSourceToStepConverter::GetConversionResult(timeout);
        }
    }
    I3CLSimStepSeriesConstPtr default_GetConversionResult(double timeout) {return this->get_override("GetConversionResult")(timeout);}

    virtual void SetLightSourceParameterizationSeries(const I3CLSimLightSourceParameterizationSeries &parameterizationSeries_)
    {
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
    bp::object operator()(const shared_ptr<const T>& v) const
    {
        return bp::object(boost::const_pointer_cast<T>(v));
    }

    // this is done for non-const shared pointers
    template <typename T>
    bp::object operator()(const shared_ptr<T>& v) const
    {
        return bp::object(v);
    }
    
    // do this for everything but shared_ptrs
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
             bp::return_internal_reference<>())
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterWrapper>, shared_ptr<const I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterWrapper>, shared_ptr<I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterWrapper>, shared_ptr<const I3CLSimLightSourceToStepConverterWrapper> >();
    
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
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterGeant4>, shared_ptr<const I3CLSimLightSourceToStepConverterGeant4> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterGeant4>, shared_ptr<I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterGeant4>, shared_ptr<const I3CLSimLightSourceToStepConverter> >();
    
    
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
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterPPC>, shared_ptr<const I3CLSimLightSourceToStepConverterPPC> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterPPC>, shared_ptr<I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterPPC>, shared_ptr<const I3CLSimLightSourceToStepConverter> >();

    
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
         I3CLSimWlenDependentValueConstPtr, I3CLSimSpectrumTablePtr, uint32_t
         >(
           (
            bp::arg("flasherSpectrumNoBias"),
            bp::arg("spectrumTable"),
            bp::arg("photonsPerStep") = I3CLSimLightSourceToStepConverterPPC::default_photonsPerStep
            )
           )
         )
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterFlasher>, shared_ptr<const I3CLSimLightSourceToStepConverterFlasher> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterFlasher>, shared_ptr<I3CLSimLightSourceToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceToStepConverterFlasher>, shared_ptr<const I3CLSimLightSourceToStepConverter> >();

}
