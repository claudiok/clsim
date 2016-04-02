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
 * @file I3CLSimScalarField.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/function/I3CLSimScalarField.h>

#include <clsim/function/I3CLSimScalarFieldConstant.h>
#include <clsim/function/I3CLSimScalarFieldAnisotropyAbsLenScaling.h>
#include <clsim/function/I3CLSimScalarFieldIceTiltZShift.h>

#include <boost/preprocessor/seq.hpp>
#include "const_ptr_helpers.h"
#include "python_gil_holder.h"


using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimScalarFieldWrapper : I3CLSimScalarField, bp::wrapper<I3CLSimScalarField>
{
    // pure virtual
    virtual bool HasNativeImplementation() const {utils::python_gil_holder gil; return this->get_override("HasNativeImplementation")();}
    virtual double GetValue(double x, double y, double z) const {utils::python_gil_holder gil; return this->get_override("GetValue")(x,y,z);}
    virtual std::string GetOpenCLFunction(const std::string &functionName) const {utils::python_gil_holder gil; return this->get_override("GetOpenCLFunction")(functionName);}
    virtual bool CompareTo(const I3CLSimScalarField &other) const {utils::python_gil_holder gil; return this->get_override("CompareTo")(other);}
};

bool I3CLSimScalarField_equalWrap(const I3CLSimScalarField &a, const I3CLSimScalarField &b)
{
    return a==b;
}

void register_I3CLSimScalarField()
{
    {
        double (I3CLSimScalarField::*GetValue_three)(double, double, double) const = &I3CLSimScalarField::GetValue;
        double (I3CLSimScalarField::*GetValue_vec)(const std::vector<double> &) const = &I3CLSimScalarField::GetValue;

        bp::scope I3CLSimScalarField_scope = 
        bp::class_<I3CLSimScalarFieldWrapper, boost::shared_ptr<I3CLSimScalarFieldWrapper>, boost::noncopyable>("I3CLSimScalarField")
        .def("HasNativeImplementation", bp::pure_virtual(&I3CLSimScalarField::HasNativeImplementation))
        .def("GetValue", bp::pure_virtual(GetValue_three))
        .def("GetValue", GetValue_vec)
        .def("GetOpenCLFunction", bp::pure_virtual(&I3CLSimScalarField::GetOpenCLFunction))
        .def("CompareTo", bp::pure_virtual(&I3CLSimScalarField::CompareTo))
        .def("__eq__", &I3CLSimScalarField_equalWrap)
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldWrapper>, boost::shared_ptr<const I3CLSimScalarField> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldWrapper>, boost::shared_ptr<I3CLSimScalarField> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldWrapper>, boost::shared_ptr<const I3CLSimScalarFieldWrapper> >();
    utils::register_const_ptr<I3CLSimScalarField>();

    // constant value
    {
        bp::class_<
        I3CLSimScalarFieldConstant, 
        boost::shared_ptr<I3CLSimScalarFieldConstant>, 
        bases<I3CLSimScalarField>,
        boost::noncopyable
        >
        (
         "I3CLSimScalarFieldConstant",
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
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldConstant>, boost::shared_ptr<const I3CLSimScalarFieldConstant> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldConstant>, boost::shared_ptr<I3CLSimScalarField> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldConstant>, boost::shared_ptr<const I3CLSimScalarField> >();
    utils::register_const_ptr<I3CLSimScalarFieldConstant>();


    // Spice-Lea anisotropy
    {
        bp::class_<
        I3CLSimScalarFieldAnisotropyAbsLenScaling, 
        boost::shared_ptr<I3CLSimScalarFieldAnisotropyAbsLenScaling>, 
        bases<I3CLSimScalarField>,
        boost::noncopyable
        >
        (
         "I3CLSimScalarFieldAnisotropyAbsLenScaling",
         bp::init<
         double, double, double
         >(
           (
            bp::arg("anisotropyDirAzimuth") = I3CLSimScalarFieldAnisotropyAbsLenScaling::default_anisotropyDirAzimuth,
            bp::arg("magnitudeAlongDir") = I3CLSimScalarFieldAnisotropyAbsLenScaling::default_magnitudeAlongDir,
            bp::arg("magnitudePerpToDir") = I3CLSimScalarFieldAnisotropyAbsLenScaling::default_magnitudePerpToDir
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldAnisotropyAbsLenScaling>, boost::shared_ptr<const I3CLSimScalarFieldAnisotropyAbsLenScaling> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldAnisotropyAbsLenScaling>, boost::shared_ptr<I3CLSimScalarField> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldConstant>, boost::shared_ptr<const I3CLSimScalarField> >();
    utils::register_const_ptr<I3CLSimScalarFieldAnisotropyAbsLenScaling>();


    // ice tilt z-shifts
    {
        bp::class_<
        I3CLSimScalarFieldIceTiltZShift, 
        boost::shared_ptr<I3CLSimScalarFieldIceTiltZShift>, 
        bases<I3CLSimScalarField>,
        boost::noncopyable
        >
        (
         "I3CLSimScalarFieldIceTiltZShift",
         bp::init<
         const std::vector<double> &, const std::vector<double> &, const I3Matrix &, double
         >(
           (
            bp::arg("distancesFromOriginAlongTilt"),
            bp::arg("zCoordinates"),
            bp::arg("zCorrections"),
            bp::arg("directionOfTiltAzimuth") = I3CLSimScalarFieldIceTiltZShift::default_directionOfTiltAzimuth
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldIceTiltZShift>, boost::shared_ptr<const I3CLSimScalarFieldIceTiltZShift> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldIceTiltZShift>, boost::shared_ptr<I3CLSimScalarField> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimScalarFieldConstant>, boost::shared_ptr<const I3CLSimScalarField> >();
    utils::register_const_ptr<I3CLSimScalarFieldIceTiltZShift>();

}
