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
 * @file I3CLSimVectorTransform.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/function/I3CLSimVectorTransform.h>

#include <clsim/function/I3CLSimVectorTransformConstant.h>
#include <clsim/function/I3CLSimVectorTransformMatrix.h>

#include <boost/preprocessor/seq.hpp>
#include "const_ptr_helpers.h"
#include "python_gil_holder.h"


using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimVectorTransformWrapper : I3CLSimVectorTransform, bp::wrapper<I3CLSimVectorTransform>
{
    // pure virtual
    virtual bool HasNativeImplementation() const {utils::python_gil_holder gil; return this->get_override("HasNativeImplementation")();}
    virtual std::vector<double> ApplyTransform(const std::vector<double> &vec) const {utils::python_gil_holder gil; return this->get_override("ApplyTransform")(vec);}
    virtual std::string GetOpenCLFunction(const std::string &functionName) const {utils::python_gil_holder gil; return this->get_override("GetOpenCLFunction")(functionName);}
    virtual bool CompareTo(const I3CLSimVectorTransform &other) const {utils::python_gil_holder gil; return this->get_override("CompareTo")(other);}
};

bool I3CLSimVectorTransform_equalWrap(const I3CLSimVectorTransform &a, const I3CLSimVectorTransform &b)
{
    return a==b;
}

void register_I3CLSimVectorTransform()
{
    {
        bp::scope I3CLSimVectorTransform_scope = 
        bp::class_<I3CLSimVectorTransformWrapper, boost::shared_ptr<I3CLSimVectorTransformWrapper>, boost::noncopyable>("I3CLSimVectorTransform")
        .def("HasNativeImplementation", bp::pure_virtual(&I3CLSimVectorTransform::HasNativeImplementation))
        .def("ApplyTransform", bp::pure_virtual(&I3CLSimVectorTransform::ApplyTransform))
        .def("GetOpenCLFunction", bp::pure_virtual(&I3CLSimVectorTransform::GetOpenCLFunction))
        .def("CompareTo", bp::pure_virtual(&I3CLSimVectorTransform::CompareTo))
        .def("__eq__", &I3CLSimVectorTransform_equalWrap)
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformWrapper>, boost::shared_ptr<const I3CLSimVectorTransform> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformWrapper>, boost::shared_ptr<I3CLSimVectorTransform> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformWrapper>, boost::shared_ptr<const I3CLSimVectorTransformWrapper> >();
    utils::register_const_ptr<I3CLSimVectorTransform>();

    // constant value
    {
        bp::class_<
        I3CLSimVectorTransformConstant, 
        boost::shared_ptr<I3CLSimVectorTransformConstant>, 
        bases<I3CLSimVectorTransform>,
        boost::noncopyable
        >
        (
         "I3CLSimVectorTransformConstant",
         bp::init<>()
         )
        ;
    }
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformConstant>, boost::shared_ptr<const I3CLSimVectorTransformConstant> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformConstant>, boost::shared_ptr<I3CLSimVectorTransform> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformConstant>, boost::shared_ptr<const I3CLSimVectorTransform> >();
    utils::register_const_ptr<I3CLSimVectorTransformConstant>();


    // Matrix transform
    {
        bp::class_<
        I3CLSimVectorTransformMatrix, 
        boost::shared_ptr<I3CLSimVectorTransformMatrix>, 
        bases<I3CLSimVectorTransform>,
        boost::noncopyable
        >
        (
         "I3CLSimVectorTransformMatrix",
         bp::init<
         const I3Matrix &, bool
         >(
           (
            bp::arg("matrix"),
            bp::arg("renormalize") = false
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformMatrix>, boost::shared_ptr<const I3CLSimVectorTransformMatrix> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformMatrix>, boost::shared_ptr<I3CLSimVectorTransform> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimVectorTransformConstant>, boost::shared_ptr<const I3CLSimVectorTransform> >();
    utils::register_const_ptr<I3CLSimVectorTransformMatrix>();

}
