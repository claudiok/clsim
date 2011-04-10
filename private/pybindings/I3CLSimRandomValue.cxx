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

#include <clsim/I3CLSimRandomValue.h>
#include <clsim/I3CLSimRandomValueHenyeyGreenstein.h>
#include <clsim/I3CLSimRandomValueRayleighScatteringCosAngle.h>
#include <clsim/I3CLSimRandomValueInterpolatedDistribution.h>
#include <clsim/I3CLSimRandomValueSimplifiedLiu.h>
#include <clsim/I3CLSimRandomValueMixed.h>
#include <clsim/I3CLSimRandomValueApplyFunction.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimRandomValueWrapper : I3CLSimRandomValue, bp::wrapper<I3CLSimRandomValue>
{
    // pure virtual
    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const 
    {
        return this->get_override("OpenCLFunctionWillOnlyUseASingleRandomNumber")();
    }
    
    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc) const 
    {
        return this->get_override("GetOpenCLFunction")(functionName,
                                                       functionArgs,
                                                       functionArgsToCall,
                                                       uniformRandomCall_co,
                                                       uniformRandomCall_oc);
    }
    
    virtual bool CompareTo(const I3CLSimRandomValue &other) const 
    {
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
        bp::class_<I3CLSimRandomValueWrapper, boost::shared_ptr<I3CLSimRandomValueWrapper>, boost::noncopyable>("I3CLSimRandomValue")
        .def("OpenCLFunctionWillOnlyUseASingleRandomNumber", bp::pure_virtual(&I3CLSimRandomValue::OpenCLFunctionWillOnlyUseASingleRandomNumber))
        .def("GetOpenCLFunction", bp::pure_virtual(&I3CLSimRandomValue::GetOpenCLFunction))
        .def("CompareTo", bp::pure_virtual(&I3CLSimRandomValue::CompareTo))
        .def("__eq__", &I3CLSimRandomValue_equalWrap)
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWrapper>, shared_ptr<const I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWrapper>, shared_ptr<I3CLSimRandomValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimRandomValueWrapper>, shared_ptr<const I3CLSimRandomValueWrapper> >();
    
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
         I3CLSimRandomValuePtr,
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

}
