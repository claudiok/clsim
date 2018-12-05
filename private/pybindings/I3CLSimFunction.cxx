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
 * @file I3CLSimFunction.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <icetray/OMKey.h>
#include <icetray/python/std_map_indexing_suite.hpp>

#include <clsim/function/I3CLSimFunction.h>

#include <clsim/function/I3CLSimFunctionConstant.h>
#include <clsim/function/I3CLSimFunctionDeltaPeak.h>

#include <clsim/function/I3CLSimFunctionFromTable.h>
#include <clsim/function/I3CLSimFunctionRefIndexQuanFry.h>
#include <clsim/function/I3CLSimFunctionScatLenPartic.h>

#include <clsim/function/I3CLSimFunctionRefIndexIceCube.h>
#include <clsim/function/I3CLSimFunctionAbsLenIceCube.h>
#include <clsim/function/I3CLSimFunctionScatLenIceCube.h>

#include <clsim/function/I3CLSimFunctionPolynomial.h>

#include <boost/preprocessor/seq.hpp>
#include "const_ptr_helpers.h"
#include "python_gil_holder.h"


using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimFunctionWrapper : I3CLSimFunction, bp::wrapper<I3CLSimFunction>
{
    // pure virtual
    virtual bool HasNativeImplementation() const {utils::python_gil_holder gil; return this->get_override("HasNativeImplementation")();}
    virtual bool HasDerivative() const {utils::python_gil_holder gil; return this->get_override("HasDerivative")();}
    virtual double GetValue(double wlen) const {utils::python_gil_holder gil; return this->get_override("GetValue")(wlen);}
    virtual double GetMinWlen() const {utils::python_gil_holder gil; return this->get_override("GetMinWlen")();}
    virtual double GetMaxWlen() const {utils::python_gil_holder gil; return this->get_override("GetMaxWlen")();}
    virtual std::string GetOpenCLFunction(const std::string &functionName) const {utils::python_gil_holder gil; return this->get_override("GetOpenCLFunction")(functionName);}
    virtual bool CompareTo(const I3CLSimFunction &other) const {utils::python_gil_holder gil; return this->get_override("CompareTo")(other);}
    
    // default implementation
    virtual double GetDerivative(double wlen) const {utils::python_gil_holder gil; if (override f = this->get_override("GetDerivative")) {return f(wlen);} else {return I3CLSimFunction::GetDerivative(wlen);}}
    virtual std::string GetOpenCLFunctionDerivative(const std::string &functionName) const {utils::python_gil_holder gil; if (override f = this->get_override("GetOpenCLFunctionDerivative")) {return f(functionName);} else {return I3CLSimFunction::GetOpenCLFunctionDerivative(functionName);}}
    
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
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionWrapper>, boost::shared_ptr<const I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionWrapper>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionWrapper>, boost::shared_ptr<const I3CLSimFunctionWrapper> >();
    utils::register_const_ptr<I3CLSimFunction>();

    // values that are different for each DOM

    { 
      typedef std::map<OMKey, I3CLSimFunctionConstPtr> I3CLSimFunctionPtrMap; 
      I3_POINTER_TYPEDEFS(I3CLSimFunctionPtrMap); 
      bp::class_<I3CLSimFunctionPtrMap, I3CLSimFunctionPtrMapPtr>("I3CLSimFunctionMap") 
	.def(bp::std_map_indexing_suite<I3CLSimFunctionPtrMap>()) 
	; 
      bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionPtrMap>, boost::shared_ptr<const I3CLSimFunctionPtrMap> >(); 
      utils::register_const_ptr<I3CLSimFunctionPtrMap>(); 
    } 

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
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionConstant>, boost::shared_ptr<const I3CLSimFunctionConstant> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionConstant>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionConstant>, boost::shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionConstant>();


    // a delta peak
    {
        bp::class_<
        I3CLSimFunctionDeltaPeak, 
        boost::shared_ptr<I3CLSimFunctionDeltaPeak>, 
        bases<I3CLSimFunction>,
        boost::noncopyable
        >
        (
         "I3CLSimFunctionDeltaPeak",
         bp::init<
         double
         >(
           (
            bp::arg("peakPosition")
            )
           )
         )
        .def("GetPeakPosition", &I3CLSimFunctionDeltaPeak::GetPeakPosition)
        ;
    }
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionDeltaPeak>, boost::shared_ptr<const I3CLSimFunctionDeltaPeak> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionDeltaPeak>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionConstant>, boost::shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionDeltaPeak>();

    
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
         const std::vector<double> &,
         bool
         >(
           (
            bp::arg("startWlen"), 
            bp::arg("wlenStep"), 
            bp::arg("values"),
            bp::arg("storeDataAsHalfPrecision") = I3CLSimFunctionFromTable::default_storeDataAsHalfPrecision
           )
          )
        )
        .def(bp::init<
             const std::vector<double> &,
             const std::vector<double> &,
             bool
             >(
               (
                bp::arg("wlens"),
                bp::arg("values"),
                bp::arg("storeDataAsHalfPrecision") = I3CLSimFunctionFromTable::default_storeDataAsHalfPrecision
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
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionFromTable>, boost::shared_ptr<const I3CLSimFunctionFromTable> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionFromTable>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionFromTable>, boost::shared_ptr<const I3CLSimFunction> >();
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
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionScatLenPartic>, boost::shared_ptr<const I3CLSimFunctionScatLenPartic> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionScatLenPartic>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionScatLenPartic>, boost::shared_ptr<const I3CLSimFunction> >();
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
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionRefIndexQuanFry>, boost::shared_ptr<const I3CLSimFunctionRefIndexQuanFry> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionRefIndexQuanFry>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionRefIndexQuanFry>, boost::shared_ptr<const I3CLSimFunction> >();
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
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionRefIndexIceCube>, boost::shared_ptr<const I3CLSimFunctionRefIndexIceCube> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionRefIndexIceCube>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionRefIndexIceCube>, boost::shared_ptr<const I3CLSimFunction> >();
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
        .def("GetKappa", &I3CLSimFunctionAbsLenIceCube::GetKappa)
        .def("GetA", &I3CLSimFunctionAbsLenIceCube::GetA)
        .def("GetB", &I3CLSimFunctionAbsLenIceCube::GetB)
        .def("GetD", &I3CLSimFunctionAbsLenIceCube::GetD)
        .def("GetE", &I3CLSimFunctionAbsLenIceCube::GetE)
        .def("GetADust400", &I3CLSimFunctionAbsLenIceCube::GetADust400)
        .def("GetDeltaTau", &I3CLSimFunctionAbsLenIceCube::GetDeltaTau)

        .add_property("kappa", &I3CLSimFunctionAbsLenIceCube::GetKappa)
        .add_property("A", &I3CLSimFunctionAbsLenIceCube::GetA)
        .add_property("B", &I3CLSimFunctionAbsLenIceCube::GetB)
        .add_property("D", &I3CLSimFunctionAbsLenIceCube::GetD)
        .add_property("E", &I3CLSimFunctionAbsLenIceCube::GetE)
        .add_property("aDust400", &I3CLSimFunctionAbsLenIceCube::GetADust400)
        .add_property("deltaTau", &I3CLSimFunctionAbsLenIceCube::GetDeltaTau)
        ;
    }
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionAbsLenIceCube>, boost::shared_ptr<const I3CLSimFunctionAbsLenIceCube> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionAbsLenIceCube>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionAbsLenIceCube>, boost::shared_ptr<const I3CLSimFunction> >();
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
        .def("GetAlpha", &I3CLSimFunctionScatLenIceCube::GetAlpha)
        .def("GetB400", &I3CLSimFunctionScatLenIceCube::GetB400)
        .add_property("alpha", &I3CLSimFunctionScatLenIceCube::GetAlpha)
        .add_property("b400", &I3CLSimFunctionScatLenIceCube::GetB400)
        ;
    }
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionScatLenIceCube>, boost::shared_ptr<const I3CLSimFunctionScatLenIceCube> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionScatLenIceCube>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionScatLenIceCube>, boost::shared_ptr<const I3CLSimFunction> >();
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
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionPolynomial>, boost::shared_ptr<const I3CLSimFunctionPolynomial> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionPolynomial>, boost::shared_ptr<I3CLSimFunction> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimFunctionPolynomial>, boost::shared_ptr<const I3CLSimFunction> >();
    utils::register_const_ptr<I3CLSimFunctionPolynomial>();

}
