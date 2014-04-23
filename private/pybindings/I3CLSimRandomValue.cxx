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
 * @file I3CLSimRandomValue.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <dataclasses/I3Vector.h>

#include <clsim/random_value/I3CLSimRandomValue.h>
#include <clsim/random_value/I3CLSimRandomValueHenyeyGreenstein.h>
#include <clsim/random_value/I3CLSimRandomValueRayleighScatteringCosAngle.h>
#include <clsim/random_value/I3CLSimRandomValueInterpolatedDistribution.h>
#include <clsim/random_value/I3CLSimRandomValueSimplifiedLiu.h>
#include <clsim/random_value/I3CLSimRandomValueMixed.h>
#include <clsim/random_value/I3CLSimRandomValueApplyFunction.h>
#include <clsim/random_value/I3CLSimRandomValueWlenCherenkovNoDispersion.h>
#include <clsim/random_value/I3CLSimRandomValueNormalDistribution.h>
#include <clsim/random_value/I3CLSimRandomValueFixParameter.h>
#include <clsim/random_value/I3CLSimRandomValueConstant.h>
#include <clsim/random_value/I3CLSimRandomValueUniform.h>

#include <icetray/python/list_indexing_suite.hpp>

#include <boost/preprocessor/seq.hpp>
#include "const_ptr_helpers.h"
#include "python_gil_holder.h"

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimRandomValueWrapper : I3CLSimRandomValue, bp::wrapper<I3CLSimRandomValue>
{
    // pure virtual
    virtual std::size_t NumberOfParameters() const
    {
        utils::python_gil_holder gil;
        return this->get_override("NumberOfParameters")();
    }

    // pure virtual
    virtual double SampleFromDistribution(const I3RandomServicePtr &random,
                                          const std::vector<double> &parameters) const
    {
        utils::python_gil_holder gil;
        return this->get_override("SampleFromDistribution")(random, parameters);
    }

    // pure virtual
    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const 
    {
        utils::python_gil_holder gil;
        return this->get_override("OpenCLFunctionWillOnlyUseASingleRandomNumber")();
    }

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc) const 
    {
        utils::python_gil_holder gil;
        return this->get_override("GetOpenCLFunction")(functionName,
                                                       functionArgs,
                                                       functionArgsToCall,
                                                       uniformRandomCall_co,
                                                       uniformRandomCall_oc);
    }

    virtual bool CompareTo(const I3CLSimRandomValue &other) const 
    {
        utils::python_gil_holder gil;
        return this->get_override("CompareTo")(other);
    }

};

bool I3CLSimRandomValue_equalWrap(const I3CLSimRandomValue &a, const I3CLSimRandomValue &b)
{
    return a==b;
}

void register_I3CLSimRandomValue()
{
    {
        bp::scope I3CLSimRandomValue_scope =
        bp::class_<I3CLSimRandomValueWrapper, shared_ptr<I3CLSimRandomValueWrapper>, boost::noncopyable>("I3CLSimRandomValue")
        .def("NumberOfParameters", bp::pure_virtual(&I3CLSimRandomValue::NumberOfParameters))
        .def("SampleFromDistribution", bp::pure_virtual(&I3CLSimRandomValue::SampleFromDistribution))
        .def("OpenCLFunctionWillOnlyUseASingleRandomNumber", bp::pure_virtual(&I3CLSimRandomValue::OpenCLFunctionWillOnlyUseASingleRandomNumber))
        .def("GetOpenCLFunction", bp::pure_virtual(&I3CLSimRandomValue::GetOpenCLFunction))
        .def("CompareTo", bp::pure_virtual(&I3CLSimRandomValue::CompareTo))
        .def("__eq__", &I3CLSimRandomValue_equalWrap)
        ;
    }

    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWrapper>, shared_ptr<const I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWrapper>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValue>, shared_ptr<const I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWrapper>, shared_ptr<const I3CLSimRandomValueWrapper> >();
    utils::register_const_ptr<I3CLSimRandomValue>();

    // Henyey-Greenstein
    {
        bp::class_<
        I3CLSimRandomValueHenyeyGreenstein, 
        boost::shared_ptr<I3CLSimRandomValueHenyeyGreenstein>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueHenyeyGreenstein",
         bp::init<
         double
         >(
           (
            bp::arg("meanCosine")
           )
          )
        )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueHenyeyGreenstein>, shared_ptr<const I3CLSimRandomValueHenyeyGreenstein> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueHenyeyGreenstein>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueHenyeyGreenstein>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueHenyeyGreenstein>();

    
    // Rayleigh scattering
    {
        bp::class_<
        I3CLSimRandomValueRayleighScatteringCosAngle, 
        boost::shared_ptr<I3CLSimRandomValueRayleighScatteringCosAngle>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueRayleighScatteringCosAngle",
         bp::init<>()
        )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueRayleighScatteringCosAngle>, shared_ptr<const I3CLSimRandomValueRayleighScatteringCosAngle> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueRayleighScatteringCosAngle>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueRayleighScatteringCosAngle>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueRayleighScatteringCosAngle>();


    // simplified Liu
    {
        bp::class_<
        I3CLSimRandomValueSimplifiedLiu, 
        boost::shared_ptr<I3CLSimRandomValueSimplifiedLiu>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueSimplifiedLiu",
         bp::init<
         double
         >(
           (
            bp::arg("meanCosine")
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueSimplifiedLiu>, shared_ptr<const I3CLSimRandomValueSimplifiedLiu> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueSimplifiedLiu>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueSimplifiedLiu>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueSimplifiedLiu>();


    // table of (x,value) pairs
    {
        bp::class_<
        I3CLSimRandomValueInterpolatedDistribution, 
        boost::shared_ptr<I3CLSimRandomValueInterpolatedDistribution>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueInterpolatedDistribution",
         bp::init<
         const std::vector<double> &,
         const std::vector<double> &
         >(
           (
            bp::arg("x"),
            bp::arg("y")
           )
          )
         )
        .def(init<double, double, const std::vector<double> & >
             (
              (
               bp::arg("xFirst"),
               bp::arg("xSpacing"),
               bp::arg("y")
               )
              )
             )

        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueInterpolatedDistribution>, shared_ptr<const I3CLSimRandomValueInterpolatedDistribution> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueInterpolatedDistribution>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueInterpolatedDistribution>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueInterpolatedDistribution>();


    // mix of two distributions
    {
        bp::class_<
        I3CLSimRandomValueMixed, 
        boost::shared_ptr<I3CLSimRandomValueMixed>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueMixed",
         bp::init<
         double,
         I3CLSimRandomValueConstPtr,
         I3CLSimRandomValueConstPtr
         >(
           (
            bp::arg("fractionOfFirstDistribution"),
            bp::arg("firstDistribution"),
            bp::arg("secondDistribution")
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueMixed>, shared_ptr<const I3CLSimRandomValueMixed> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueMixed>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueMixed>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueMixed>();


    // apply function to generated random value
    {
        bp::class_<
        I3CLSimRandomValueApplyFunction, 
        boost::shared_ptr<I3CLSimRandomValueApplyFunction>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueApplyFunction",
         bp::init<
         const I3CLSimRandomValuePtr &,
         const std::string &
         >(
           (
            bp::arg("randomDistUsed"),
            bp::arg("functionName")
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueApplyFunction>, shared_ptr<const I3CLSimRandomValueApplyFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueApplyFunction>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueApplyFunction>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueApplyFunction>();


    // wavelength distributed according to a Cherenkov spectrum (no dispersion)
    {
        bp::class_<
        I3CLSimRandomValueWlenCherenkovNoDispersion, 
        boost::shared_ptr<I3CLSimRandomValueWlenCherenkovNoDispersion>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueWlenCherenkovNoDispersion",
         bp::init<
         double,
         double
         >(
           (
            bp::arg("fromWlen"),
            bp::arg("toWlen")
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWlenCherenkovNoDispersion>, shared_ptr<const I3CLSimRandomValueWlenCherenkovNoDispersion> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWlenCherenkovNoDispersion>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWlenCherenkovNoDispersion>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueWlenCherenkovNoDispersion>();


    // sampling from a gaussian
    {
        bp::class_<
        I3CLSimRandomValueNormalDistribution, 
        boost::shared_ptr<I3CLSimRandomValueNormalDistribution>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueNormalDistribution",
         bp::init<
         >(
           //(
           //)
          )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueNormalDistribution>, shared_ptr<const I3CLSimRandomValueNormalDistribution> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueNormalDistribution>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueNormalDistribution>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueNormalDistribution>();

    
    // set a run-time parameter in a distribution to a fixed value
    {
        bp::class_<
        I3CLSimRandomValueFixParameter, 
        boost::shared_ptr<I3CLSimRandomValueFixParameter>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueFixParameter",
         bp::init<
         const I3CLSimRandomValuePtr &,
         std::size_t,
         double
         >(
           (
            bp::arg("randomDistUsed"),
            bp::arg("parameterIndex"),
            bp::arg("parameterValue")
           )
          )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueFixParameter>, shared_ptr<const I3CLSimRandomValueFixParameter> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueFixParameter>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueFixParameter>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueFixParameter>();

    
    // a single, constant value
    {
        bp::class_<
        I3CLSimRandomValueConstant, 
        boost::shared_ptr<I3CLSimRandomValueConstant>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueConstant",
         bp::init<
         double
         >(
           (
            bp::arg("value")
            )
           )
         )
        .def(init<>()) // this one also has a default constructor
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueConstant>, shared_ptr<const I3CLSimRandomValueConstant> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueConstant>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueConstant>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueConstant>();

    
    // a uniformly distributed value
    {
        bp::class_<
        I3CLSimRandomValueUniform, 
        boost::shared_ptr<I3CLSimRandomValueUniform>, 
        bases<I3CLSimRandomValue>,
        boost::noncopyable
        >
        (
         "I3CLSimRandomValueUniform",
         bp::init<
         double, double
         >(
           (
            bp::arg("from"),
            bp::arg("to")
            )
           )
         )
        .def(init<>()) // this one also has a default constructor
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueUniform>, shared_ptr<const I3CLSimRandomValueUniform> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueUniform>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueUniform>, shared_ptr<const I3CLSimRandomValue> >();
    utils::register_const_ptr<I3CLSimRandomValueUniform>();


    // a vector of distributions
    {
        typedef std::vector<I3CLSimRandomValueConstPtr> Series;
        bp::class_<Series, shared_ptr<Series> >("I3CLSimRandomValuePtrSeries")
            .def(bp::list_indexing_suite<Series>())
        ;
    }

}
