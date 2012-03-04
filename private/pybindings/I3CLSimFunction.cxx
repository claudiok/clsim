//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module clsim.
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

#include <sstream>

#include <clsim/function/I3CLSimFunction.h>

#include <clsim/function/I3CLSimFunctionConstant.h>

#include <clsim/function/I3CLSimFunctionFromTable.h>
#include <clsim/function/I3CLSimFunctionRefIndexQuanFry.h>
#include <clsim/function/I3CLSimFunctionScatLenPartic.h>

#include <clsim/function/I3CLSimFunctionRefIndexIceCube.h>
#include <clsim/function/I3CLSimFunctionAbsLenIceCube.h>
#include <clsim/function/I3CLSimFunctionScatLenIceCube.h>

#include <clsim/function/I3CLSimFunctionPolynomial.h>

#include <boost/preprocessor/seq.hpp>
#include "const_ptr_helpers.h"


using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimFunctionWrapper : I3CLSimFunction, bp::wrapper<I3CLSimFunction>
{
    // pure virtual
    virtual bool HasNativeImplementation() const {return this->get_override("HasNativeImplementation")();}
    virtual bool HasDerivative() const {return this->get_override("HasDerivative")();}
    virtual double GetValue(double wlen) const {return this->get_override("GetValue")(wlen);}
    virtual double GetMinWlen() const {return this->get_override("GetMinWlen")();}
    virtual double GetMaxWlen() const {return this->get_override("GetMaxWlen")();}
    virtual std::string GetOpenCLFunction(const std::string &functionName) const {return this->get_override("GetOpenCLFunction")(functionName);}
    virtual bool CompareTo(const I3CLSimFunction &other) const {return this->get_override("CompareTo")(other);}
    
    // default implementation
    virtual double GetDerivative(double wlen) const {if (override f = this->get_override("GetDerivative")) {return f(wlen);} else {return I3CLSimFunction::GetDerivative(wlen);}}
    virtual std::string GetOpenCLFunctionDerivative(const std::string &functionName) const {if (override f = this->get_override("GetOpenCLFunctionDerivative")) {return f(functionName);} else {return I3CLSimFunction::GetOpenCLFunctionDerivative(functionName);}}
    
    double default_GetDerivative(double wlen) const {return this->I3CLSimFunction::GetDerivative(wlen);}
    std::string default_GetOpenCLFunctionDerivative(const std::string &functionName) const {return this->I3CLSimFunction::GetOpenCLFunctionDerivative(functionName);}
    
};

bool I3CLSimFunction_equalWrap(const I3CLSimFunction &a, const I3CLSimFunction &b)
{
    return a==b;
}

void register_I3CLSimFunction()
{
    {
        bp::scope I3CLSimFunction_scope = 
        bp::class_<I3CLSimFunctionWrapper, boost::shared_ptr<I3CLSimFunctionWrapper>, boost::noncopyable>("I3CLSimFunction")
        .def("HasNativeImplementation", bp::pure_virtual(&I3CLSimFunction::HasNativeImplementation))
        .def("HasDerivative", bp::pure_virtual(&I3CLSimFunction::HasDerivative))
        .def("GetValue", bp::pure_virtual(&I3CLSimFunction::GetValue))
        .def("GetMinWlen", bp::pure_virtual(&I3CLSimFunction::GetMinWlen))
        .def("GetMaxWlen", bp::pure_virtual(&I3CLSimFunction::GetMaxWlen))
        .def("GetOpenCLFunction", bp::pure_virtual(&I3CLSimFunction::GetOpenCLFunction))
        .def("CompareTo", bp::pure_virtual(&I3CLSimFunction::CompareTo))
        .def("__eq__", &I3CLSimFunction_equalWrap)
        
        .def("GetDerivative", &I3CLSimFunction::GetDerivative, &I3CLSimFunctionWrapper::default_GetDerivative)
        .def("GetOpenCLFunctionDerivative", &I3CLSimFunction::GetOpenCLFunctionDerivative, &I3CLSimFunctionWrapper::default_GetOpenCLFunctionDerivative)
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionWrapper>, shared_ptr<const I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionWrapper>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionWrapper>, shared_ptr<const I3CLSimFunctionWrapper> >();
    utils::register_const_ptr<I3CLSimFunction>();

    // constant value
    {
        bp::class_<
        I3CLSimFunctionConstant, 
        boost::shared_ptr<I3CLSimFunctionConstant>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionConstant",
         bp::init<
         double
         >(
           (
            bp::arg("value")
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionConstant>, shared_ptr<const I3CLSimFunctionConstant> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionConstant>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionConstant>, shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionConstant>();

    
    // from table
    {
        bp::class_<
        I3CLSimFunctionFromTable, 
        boost::shared_ptr<I3CLSimFunctionFromTable>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionFromTable",
         bp::init<
         double,
         double,
         const std::vector<double> &
         >(
           (
            bp::arg("startWlen"), 
            bp::arg("wlenStep"), 
            bp::arg("values")
           )
          )
        )
        .def(bp::init<
             const std::vector<double> &,
             const std::vector<double> &
             >(
               (
                bp::arg("wlens"),
                bp::arg("values")
               )
              )
            )
        .def("GetFirstWavelength", &I3CLSimFunctionFromTable::GetFirstWavelength)
        .def("GetWavelengthStepping", &I3CLSimFunctionFromTable::GetWavelengthStepping)
        .def("GetNumEntries", &I3CLSimFunctionFromTable::GetNumEntries)
        .def("GetEntryValue", &I3CLSimFunctionFromTable::GetEntryValue)
        .def("GetEntryWavelength", &I3CLSimFunctionFromTable::GetEntryWavelength)
        .def("GetInEqualSpacingMode", &I3CLSimFunctionFromTable::GetInEqualSpacingMode)
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionFromTable>, shared_ptr<const I3CLSimFunctionFromTable> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionFromTable>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionFromTable>, shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionFromTable>();

    
    // scattering length (partic)
    {
        bp::class_<
        I3CLSimFunctionScatLenPartic, 
        boost::shared_ptr<I3CLSimFunctionScatLenPartic>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionScatLenPartic",
         bp::init<
         double,
         double
         >(
           (
            bp::arg("volumeConcentrationSmallParticles") = I3CLSimFunctionScatLenPartic::default_volumeConcentrationSmallParticles, 
            bp::arg("volumeConcentrationLargeParticles") = I3CLSimFunctionScatLenPartic::default_volumeConcentrationLargeParticles
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionScatLenPartic>, shared_ptr<const I3CLSimFunctionScatLenPartic> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionScatLenPartic>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionScatLenPartic>, shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionScatLenPartic>();


    // refractive index (Quan&Fry)
    {
        bp::class_<
        I3CLSimFunctionRefIndexQuanFry, 
        boost::shared_ptr<I3CLSimFunctionRefIndexQuanFry>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionRefIndexQuanFry",
         bp::init<
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double
         >(
           (
            bp::arg("salinity") = I3CLSimFunctionRefIndexQuanFry::default_salinity, 
            bp::arg("temperature") = I3CLSimFunctionRefIndexQuanFry::default_temperature,
            bp::arg("pressure") = I3CLSimFunctionRefIndexQuanFry::default_pressure,
            bp::arg("n0") = I3CLSimFunctionRefIndexQuanFry::default_n0,
            bp::arg("n1") = I3CLSimFunctionRefIndexQuanFry::default_n1,
            bp::arg("n2") = I3CLSimFunctionRefIndexQuanFry::default_n2,
            bp::arg("n3") = I3CLSimFunctionRefIndexQuanFry::default_n3,
            bp::arg("n4") = I3CLSimFunctionRefIndexQuanFry::default_n4,
            bp::arg("n5") = I3CLSimFunctionRefIndexQuanFry::default_n5,
            bp::arg("n6") = I3CLSimFunctionRefIndexQuanFry::default_n6,
            bp::arg("n7") = I3CLSimFunctionRefIndexQuanFry::default_n7,
            bp::arg("n8") = I3CLSimFunctionRefIndexQuanFry::default_n8,
            bp::arg("n9") = I3CLSimFunctionRefIndexQuanFry::default_n9,
            bp::arg("n10") = I3CLSimFunctionRefIndexQuanFry::default_n10
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionRefIndexQuanFry>, shared_ptr<const I3CLSimFunctionRefIndexQuanFry> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionRefIndexQuanFry>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionRefIndexQuanFry>, shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionRefIndexQuanFry>();

    
    
    // refractive index (IceCube)
    {
        bp::class_<
        I3CLSimFunctionRefIndexIceCube, 
        boost::shared_ptr<I3CLSimFunctionRefIndexIceCube>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionRefIndexIceCube",
         bp::init<
         std::string,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double,
         double
         >(
           (
            bp::arg("mode") = I3CLSimFunctionRefIndexIceCube::default_mode,
            bp::arg("n0") = I3CLSimFunctionRefIndexIceCube::default_n0,
            bp::arg("n1") = I3CLSimFunctionRefIndexIceCube::default_n1,
            bp::arg("n2") = I3CLSimFunctionRefIndexIceCube::default_n2,
            bp::arg("n3") = I3CLSimFunctionRefIndexIceCube::default_n3,
            bp::arg("n4") = I3CLSimFunctionRefIndexIceCube::default_n4,
            bp::arg("g0") = I3CLSimFunctionRefIndexIceCube::default_g0,
            bp::arg("g1") = I3CLSimFunctionRefIndexIceCube::default_g1,
            bp::arg("g2") = I3CLSimFunctionRefIndexIceCube::default_g2,
            bp::arg("g3") = I3CLSimFunctionRefIndexIceCube::default_g3,
            bp::arg("g4") = I3CLSimFunctionRefIndexIceCube::default_g4
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionRefIndexIceCube>, shared_ptr<const I3CLSimFunctionRefIndexIceCube> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionRefIndexIceCube>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionRefIndexIceCube>, shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionRefIndexIceCube>();

    
    // absorption length (IceCube)
    {
        bp::class_<
        I3CLSimFunctionAbsLenIceCube, 
        boost::shared_ptr<I3CLSimFunctionAbsLenIceCube>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionAbsLenIceCube",
         bp::init<
         double,
         double,
         double,
         double,
         double,
         double,
         double
         >(
           (
            bp::arg("kappa"),
            bp::arg("A"),
            bp::arg("B"),
            bp::arg("D"),
            bp::arg("E"),
            bp::arg("aDust400"),
            bp::arg("deltaTau")
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionAbsLenIceCube>, shared_ptr<const I3CLSimFunctionAbsLenIceCube> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionAbsLenIceCube>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionAbsLenIceCube>, shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionAbsLenIceCube>();


    // scattering length (IceCube)
    {
        bp::class_<
        I3CLSimFunctionScatLenIceCube, 
        boost::shared_ptr<I3CLSimFunctionScatLenIceCube>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionScatLenIceCube",
         bp::init<
         double,
         double
         >(
           (
            bp::arg("alpha"),
            bp::arg("b400")
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionScatLenIceCube>, shared_ptr<const I3CLSimFunctionScatLenIceCube> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionScatLenIceCube>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionScatLenIceCube>, shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionScatLenIceCube>();

    
    // polynomial
    {
        bp::class_<
        I3CLSimFunctionPolynomial, 
        boost::shared_ptr<I3CLSimFunctionPolynomial>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionPolynomial",
         bp::init<
         const std::vector<double> &
         >(
           (
            bp::arg("coeffs")
            )
           )
         )
        .def(init<const std::vector<double> &, double, double>
             (
              (
               bp::arg("coeffs"),
               bp::arg("rangemin"),
               bp::arg("rangemax")
              )
             )
            )
        .def(init<const std::vector<double> &, double, double, double, double>
             (
              (
               bp::arg("coeffs"),
               bp::arg("rangemin"),
               bp::arg("rangemax"),
               bp::arg("underflow"),
               bp::arg("overflow")
               )
              )
             )
        .def("GetCoefficients", &I3CLSimFunctionPolynomial::GetCoefficients, bp::return_value_policy<bp::copy_const_reference>())
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionPolynomial>, shared_ptr<const I3CLSimFunctionPolynomial> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionPolynomial>, shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimFunctionPolynomial>, shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionPolynomial>();

}
