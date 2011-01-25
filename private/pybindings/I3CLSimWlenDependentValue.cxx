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

#include <clsim/I3CLSimWlenDependentValue.h>
#include <clsim/I3CLSimWlenDependentValueFromTable.h>
#include <clsim/I3CLSimWlenDependentValueRefIndexQuanFry.h>
#include <clsim/I3CLSimWlenDependentValueScatLenPartic.h>

#include <clsim/I3CLSimWlenDependentValueRefIndexIceCube.h>
#include <clsim/I3CLSimWlenDependentValueAbsLenIceCube.h>
#include <clsim/I3CLSimWlenDependentValueScatLenIceCube.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimWlenDependentValueWrapper : I3CLSimWlenDependentValue, bp::wrapper<I3CLSimWlenDependentValue>
{
    // pure virtual
    virtual bool HasNativeImplementation() const {return this->get_override("HasNativeImplementation")();}
    virtual bool HasDerivative() const {return this->get_override("HasDerivative")();}
    virtual double GetValue(double wlen) const {return this->get_override("GetValue")(wlen);}
    virtual double GetMinWlen() const {return this->get_override("GetMinWlen")();}
    virtual double GetMaxWlen() const {return this->get_override("GetMaxWlen")();}
    virtual std::string GetOpenCLFunction(const std::string &functionName) const {return this->get_override("GetOpenCLFunction")(functionName);}
    virtual bool CompareTo(const I3CLSimWlenDependentValue &other) const {return this->get_override("CompareTo")(other);}
    
    // default implementation
    virtual double GetDerivative(double wlen) const {if (override f = this->get_override("GetDerivative")) {return f(wlen);} else {return I3CLSimWlenDependentValue::GetDerivative(wlen);}}
    virtual std::string GetOpenCLFunctionDerivative(const std::string &functionName) const {if (override f = this->get_override("GetOpenCLFunctionDerivative")) {return f(functionName);} else {return I3CLSimWlenDependentValue::GetOpenCLFunctionDerivative(functionName);}}
    
    double default_GetDerivative(double wlen) const {return this->I3CLSimWlenDependentValue::GetDerivative(wlen);}
    std::string default_GetOpenCLFunctionDerivative(const std::string &functionName) const {return this->I3CLSimWlenDependentValue::GetOpenCLFunctionDerivative(functionName);}
    
};

bool I3CLSimWlenDependentValue_equalWrap(const I3CLSimWlenDependentValue &a, const I3CLSimWlenDependentValue &b)
{
    return a==b;
}

void register_I3CLSimWlenDependentValue()
{
    {
        bp::scope I3CLSimWlenDependentValue_scope = 
        bp::class_<I3CLSimWlenDependentValueWrapper, boost::shared_ptr<I3CLSimWlenDependentValueWrapper>, boost::noncopyable>("I3CLSimWlenDependentValue")
        .def("HasNativeImplementation", bp::pure_virtual(&I3CLSimWlenDependentValue::HasNativeImplementation))
        .def("HasDerivative", bp::pure_virtual(&I3CLSimWlenDependentValue::HasDerivative))
        .def("GetValue", bp::pure_virtual(&I3CLSimWlenDependentValue::GetValue))
        .def("GetMinWlen", bp::pure_virtual(&I3CLSimWlenDependentValue::GetMinWlen))
        .def("GetMaxWlen", bp::pure_virtual(&I3CLSimWlenDependentValue::GetMaxWlen))
        .def("GetOpenCLFunction", bp::pure_virtual(&I3CLSimWlenDependentValue::GetOpenCLFunction))
        .def("CompareTo", bp::pure_virtual(&I3CLSimWlenDependentValue::CompareTo))
        .def("__eq__", &I3CLSimWlenDependentValue_equalWrap)
        
        .def("GetDerivative", &I3CLSimWlenDependentValue::GetDerivative, &I3CLSimWlenDependentValueWrapper::default_GetDerivative)
        .def("GetOpenCLFunctionDerivative", &I3CLSimWlenDependentValue::GetOpenCLFunctionDerivative, &I3CLSimWlenDependentValueWrapper::default_GetOpenCLFunctionDerivative)
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueWrapper>, shared_ptr<const I3CLSimWlenDependentValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueWrapper>, shared_ptr<I3CLSimWlenDependentValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueWrapper>, shared_ptr<const I3CLSimWlenDependentValueWrapper> >();
    
    // from table
    {
        bp::class_<
        I3CLSimWlenDependentValueFromTable, 
        boost::shared_ptr<I3CLSimWlenDependentValueFromTable>, 
        bases<I3CLSimWlenDependentValue>,
        boost::noncopyable
        >
        (
         "I3CLSimWlenDependentValueFromTable",
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
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueFromTable>, shared_ptr<const I3CLSimWlenDependentValueFromTable> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueFromTable>, shared_ptr<I3CLSimWlenDependentValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueFromTable>, shared_ptr<const I3CLSimWlenDependentValue> >();

    
    // scattering length (partic)
    {
        bp::class_<
        I3CLSimWlenDependentValueScatLenPartic, 
        boost::shared_ptr<I3CLSimWlenDependentValueScatLenPartic>, 
        bases<I3CLSimWlenDependentValue>,
        boost::noncopyable
        >
        (
         "I3CLSimWlenDependentValueScatLenPartic",
         bp::init<
         double,
         double
         >(
           (
            bp::arg("volumeConcentrationSmallParticles") = I3CLSimWlenDependentValueScatLenPartic::default_volumeConcentrationSmallParticles, 
            bp::arg("volumeConcentrationLargeParticles") = I3CLSimWlenDependentValueScatLenPartic::default_volumeConcentrationLargeParticles
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueScatLenPartic>, shared_ptr<const I3CLSimWlenDependentValueScatLenPartic> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueScatLenPartic>, shared_ptr<I3CLSimWlenDependentValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueScatLenPartic>, shared_ptr<const I3CLSimWlenDependentValue> >();


    // refractive index (Quan&Fry)
    {
        bp::class_<
        I3CLSimWlenDependentValueRefIndexQuanFry, 
        boost::shared_ptr<I3CLSimWlenDependentValueRefIndexQuanFry>, 
        bases<I3CLSimWlenDependentValue>,
        boost::noncopyable
        >
        (
         "I3CLSimWlenDependentValueRefIndexQuanFry",
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
            bp::arg("salinity") = I3CLSimWlenDependentValueRefIndexQuanFry::default_salinity, 
            bp::arg("temperature") = I3CLSimWlenDependentValueRefIndexQuanFry::default_temperature,
            bp::arg("pressure") = I3CLSimWlenDependentValueRefIndexQuanFry::default_pressure,
            bp::arg("n0") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n0,
            bp::arg("n1") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n1,
            bp::arg("n2") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n2,
            bp::arg("n3") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n3,
            bp::arg("n4") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n4,
            bp::arg("n5") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n5,
            bp::arg("n6") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n6,
            bp::arg("n7") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n7,
            bp::arg("n8") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n8,
            bp::arg("n9") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n9,
            bp::arg("n10") = I3CLSimWlenDependentValueRefIndexQuanFry::default_n10
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueRefIndexQuanFry>, shared_ptr<const I3CLSimWlenDependentValueRefIndexQuanFry> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueRefIndexQuanFry>, shared_ptr<I3CLSimWlenDependentValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueRefIndexQuanFry>, shared_ptr<const I3CLSimWlenDependentValue> >();

    // refractive index (IceCube)
    {
        bp::class_<
        I3CLSimWlenDependentValueRefIndexIceCube, 
        boost::shared_ptr<I3CLSimWlenDependentValueRefIndexIceCube>, 
        bases<I3CLSimWlenDependentValue>,
        boost::noncopyable
        >
        (
         "I3CLSimWlenDependentValueRefIndexIceCube",
         bp::init<
         double,
         double,
         double,
         double,
         double
         >(
           (
            bp::arg("n0") = I3CLSimWlenDependentValueRefIndexIceCube::default_n0,
            bp::arg("n1") = I3CLSimWlenDependentValueRefIndexIceCube::default_n1,
            bp::arg("n2") = I3CLSimWlenDependentValueRefIndexIceCube::default_n2,
            bp::arg("n3") = I3CLSimWlenDependentValueRefIndexIceCube::default_n3,
            bp::arg("n4") = I3CLSimWlenDependentValueRefIndexIceCube::default_n4
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueRefIndexIceCube>, shared_ptr<const I3CLSimWlenDependentValueRefIndexIceCube> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueRefIndexIceCube>, shared_ptr<I3CLSimWlenDependentValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueRefIndexIceCube>, shared_ptr<const I3CLSimWlenDependentValue> >();

    
    // absorption length (IceCube)
    {
        bp::class_<
        I3CLSimWlenDependentValueAbsLenIceCube, 
        boost::shared_ptr<I3CLSimWlenDependentValueAbsLenIceCube>, 
        bases<I3CLSimWlenDependentValue>,
        boost::noncopyable
        >
        (
         "I3CLSimWlenDependentValueAbsLenIceCube",
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
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueAbsLenIceCube>, shared_ptr<const I3CLSimWlenDependentValueAbsLenIceCube> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueAbsLenIceCube>, shared_ptr<I3CLSimWlenDependentValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueAbsLenIceCube>, shared_ptr<const I3CLSimWlenDependentValue> >();
    

    // scattering length (IceCube)
    {
        bp::class_<
        I3CLSimWlenDependentValueScatLenIceCube, 
        boost::shared_ptr<I3CLSimWlenDependentValueScatLenIceCube>, 
        bases<I3CLSimWlenDependentValue>,
        boost::noncopyable
        >
        (
         "I3CLSimWlenDependentValueScatLenIceCube",
         bp::init<
         double,
         double
         >(
           (
            bp::arg("alpha"),
            bp::arg("be400")
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueScatLenIceCube>, shared_ptr<const I3CLSimWlenDependentValueScatLenIceCube> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueScatLenIceCube>, shared_ptr<I3CLSimWlenDependentValue> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimWlenDependentValueScatLenIceCube>, shared_ptr<const I3CLSimWlenDependentValue> >();
    
    
    
}
