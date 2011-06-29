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

#include <clsim/I3CLSimParticleToStepConverter.h>
#include <clsim/I3CLSimParticleToStepConverterGeant4.h>
#include <clsim/I3CLSimParticleToStepConverterPPC.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_const.hpp>

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimParticleToStepConverterWrapper : I3CLSimParticleToStepConverter, bp::wrapper<I3CLSimParticleToStepConverter>
{
    // pure virtual
    virtual void SetBunchSizeGranularity(uint64_t num) {this->get_override("SetBunchSizeGranularity")(num);}
    virtual void SetMaxBunchSize(uint64_t num) {this->get_override("SetMaxBunchSize")(num);}
    virtual void SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias) {this->get_override("SetWlenBias")(wlenBias);}
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) {this->get_override("SetMediumProperties")(mediumProperties);}
    virtual void Initialize() {this->get_override("Initialize")();}
    virtual bool IsInitialized() const {return this->get_override("IsInitialized")();}
    virtual void EnqueueParticle(const I3Particle &particle, uint32_t identifier) {this->get_override("EnqueueParticle")(particle, identifier);}
    virtual void EnqueueBarrier() {this->get_override("EnqueueParticle")();}
    virtual bool BarrierActive() const {return this->get_override("BarrierActive")();}
    virtual bool MoreStepsAvailable() const {return this->get_override("MoreStepsAvailable")();}
    
    virtual I3CLSimStepSeriesConstPtr GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
    {return this->get_override("GetConversionResultWithBarrierInfo")(barrierWasReset, timeout);}

    virtual I3CLSimStepSeriesConstPtr GetConversionResult(double timeout)
    {
        if (override f = this->get_override("GetConversionResult")) {
            return f(timeout);
        } else {
            return I3CLSimParticleToStepConverter::GetConversionResult(timeout);
        }
    }
    I3CLSimStepSeriesConstPtr default_GetConversionResult(double timeout) {return this->get_override("GetConversionResult")(timeout);}

    virtual void SetParticleParameterizationSeries(const I3CLSimParticleParameterizationSeries &parameterizationSeries_)
    {
        if (override f = this->get_override("SetParticleParameterizationSeries")) {
            f(parameterizationSeries_);
        } else {
            I3CLSimParticleToStepConverter::SetParticleParameterizationSeries(parameterizationSeries_);
        }
    }
    void default_SetParticleParameterizationSeries(const I3CLSimParticleParameterizationSeries &parameterizationSeries_) 
    {this->I3CLSimParticleToStepConverter::SetParticleParameterizationSeries(parameterizationSeries_);}


    virtual const I3CLSimParticleParameterizationSeries &GetParticleParameterizationSeries() const
    {
        if (override f = this->get_override("GetParticleParameterizationSeries")) {
            return f();
        } else {
            return I3CLSimParticleToStepConverter::GetParticleParameterizationSeries();
        }
    }
    const I3CLSimParticleParameterizationSeries &default_GetParticleParameterizationSeries() const
    {return this->I3CLSimParticleToStepConverter::GetParticleParameterizationSeries();}

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

/// to-python convert to I3CLSimParticleToStepConverter::ConversionResult_t
struct I3CLSimParticleToStepConverter_ConversionResult_t_to_python
{
    static PyObject *convert(const I3CLSimParticleToStepConverter::ConversionResult_t& val)
    {
        return boost::python::incref(boost::apply_visitor(vc(), val).ptr()); 
    }
    
};
*/

void register_I3CLSimParticleToStepConverter()
{
    {
        bp::scope I3CLSimParticleToStepConverter_scope = 
        bp::class_<I3CLSimParticleToStepConverterWrapper, boost::shared_ptr<I3CLSimParticleToStepConverterWrapper>, boost::noncopyable>("I3CLSimParticleToStepConverter", bp::no_init)
        .def("SetBunchSizeGranularity", bp::pure_virtual(&I3CLSimParticleToStepConverter::SetBunchSizeGranularity))
        .def("SetMaxBunchSize", bp::pure_virtual(&I3CLSimParticleToStepConverter::SetMaxBunchSize))
        .def("SetWlenBias", bp::pure_virtual(&I3CLSimParticleToStepConverter::SetWlenBias))
        .def("SetMediumProperties", bp::pure_virtual(&I3CLSimParticleToStepConverter::SetMediumProperties))
        .def("Initialize", bp::pure_virtual(&I3CLSimParticleToStepConverter::Initialize))
        .def("IsInitialized", bp::pure_virtual(&I3CLSimParticleToStepConverter::IsInitialized))
        .def("EnqueueParticle", bp::pure_virtual(&I3CLSimParticleToStepConverter::EnqueueParticle))
        .def("EnqueueBarrier", bp::pure_virtual(&I3CLSimParticleToStepConverter::EnqueueBarrier))
        .def("BarrierActive", bp::pure_virtual(&I3CLSimParticleToStepConverter::BarrierActive))
        .def("MoreStepsAvailable", bp::pure_virtual(&I3CLSimParticleToStepConverter::MoreStepsAvailable))
        .def("GetConversionResult", bp::pure_virtual(&I3CLSimParticleToStepConverter::GetConversionResult), bp::arg("timeout")=NAN)
        .def("GetConversionResultWithBarrierInfo", 
             &I3CLSimParticleToStepConverter::GetConversionResultWithBarrierInfo,
             &I3CLSimParticleToStepConverterWrapper::GetConversionResultWithBarrierInfo,
             bp::arg("barrierWasReset"), bp::arg("timeout")=NAN)

        .def("SetParticleParameterizationSeries", 
             &I3CLSimParticleToStepConverter::SetParticleParameterizationSeries,
             &I3CLSimParticleToStepConverterWrapper::default_SetParticleParameterizationSeries)

        .def("GetParticleParameterizationSeries", 
             &I3CLSimParticleToStepConverter::GetParticleParameterizationSeries,
             &I3CLSimParticleToStepConverterWrapper::default_GetParticleParameterizationSeries,
             bp::return_internal_reference<>())
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterWrapper>, shared_ptr<const I3CLSimParticleToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterWrapper>, shared_ptr<I3CLSimParticleToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterWrapper>, shared_ptr<const I3CLSimParticleToStepConverterWrapper> >();
    
    //bp::to_python_converter<I3CLSimParticleToStepConverter::ConversionResult_t, I3CLSimParticleToStepConverter_ConversionResult_t_to_python>();
    
    // I3CLSimParticleToStepConverterGeant4
    {
        bp::class_<
        I3CLSimParticleToStepConverterGeant4, 
        boost::shared_ptr<I3CLSimParticleToStepConverterGeant4>, 
        bases<I3CLSimParticleToStepConverter>,
        boost::noncopyable
        >
        (
         "I3CLSimParticleToStepConverterGeant4",
         bp::init<
         uint32_t,
         std::string,
         double,
         uint32_t,
         uint32_t
         >(
           (
            bp::arg("randomSeed"),
            bp::arg("physicsListName") = I3CLSimParticleToStepConverterGeant4::default_physicsListName,
            bp::arg("maxBetaChangePerStep") = I3CLSimParticleToStepConverterGeant4::default_maxBetaChangePerStep,
            bp::arg("maxNumPhotonsPerStep") = I3CLSimParticleToStepConverterGeant4::default_maxNumPhotonsPerStep,
            bp::arg("maxQueueItems") = I3CLSimParticleToStepConverterGeant4::default_maxQueueItems
           )
          )
        )
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterGeant4>, shared_ptr<const I3CLSimParticleToStepConverterGeant4> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterGeant4>, shared_ptr<I3CLSimParticleToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterGeant4>, shared_ptr<const I3CLSimParticleToStepConverter> >();
    
    
    // I3CLSimParticleToStepConverterPPC
    {
        bp::class_<
        I3CLSimParticleToStepConverterPPC, 
        boost::shared_ptr<I3CLSimParticleToStepConverterPPC>, 
        bases<I3CLSimParticleToStepConverter>,
        boost::noncopyable
        >
        (
         "I3CLSimParticleToStepConverterPPC",
         bp::init<
         I3RandomServicePtr,
         uint32_t, uint32_t, double
         >(
           (
            bp::arg("randomService"),
            bp::arg("photonsPerStep") = I3CLSimParticleToStepConverterPPC::default_photonsPerStep,
            bp::arg("highPhotonsPerStep") = I3CLSimParticleToStepConverterPPC::default_highPhotonsPerStep,
            bp::arg("useHighPhotonsPerStepStartingFromNumPhotons") = I3CLSimParticleToStepConverterPPC::default_useHighPhotonsPerStepStartingFromNumPhotons
            )
           )
         )
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterPPC>, shared_ptr<const I3CLSimParticleToStepConverterPPC> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterPPC>, shared_ptr<I3CLSimParticleToStepConverter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleToStepConverterPPC>, shared_ptr<const I3CLSimParticleToStepConverter> >();

}
