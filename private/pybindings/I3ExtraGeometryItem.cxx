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
 * @file I3ExtraGeometryItem.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/shadow/I3ExtraGeometryItem.h>
#include <clsim/shadow/I3ExtraGeometryItemUnion.h>
#include <clsim/shadow/I3ExtraGeometryItemMove.h>
#include <clsim/shadow/I3ExtraGeometryItemCylinder.h>

#include <icetray/python/list_indexing_suite.hpp>
#include <icetray/python/dataclass_suite.hpp>

#include <boost/preprocessor/seq.hpp>
#include "const_ptr_helpers.h"

#include "python_gil_holder.h"


using namespace boost::python;
namespace bp = boost::python;

struct I3ExtraGeometryItemWrapper : I3ExtraGeometryItem, bp::wrapper<I3ExtraGeometryItem>
{
    // pure virtual
    virtual bool DoesLineIntersect(const I3Position &lineStart,
                                   const I3Position &lineEnd) const
    {
        utils::python_gil_holder gil;
        return this->get_override("DoesLineIntersect")(lineStart, lineEnd);
    }

    virtual std::pair<I3Position, I3Position> GetBoundingBox() const
    {
        utils::python_gil_holder gil;
        return this->get_override("GetBoundingBox")();
    }
};

namespace I3ExtraGeometryItem_utils
{
    template <typename T1, typename T2> 
    struct std_pair_to_tuple 
    { 
        static PyObject* convert(std::pair<T1, T2> const& p) 
        { 
            return boost::python::incref(boost::python::make_tuple(p.first, p.second).ptr()); 
        } 
    }; 

    template <typename T1, typename T2> 
    struct std_pair_to_python_converter 
    { 
        std_pair_to_python_converter() 
        { 
            boost::python::to_python_converter< 
            std::pair<T1, T2>, 
            std_pair_to_tuple<T1, T2> >(); 
        } 
    };
    
    
    static boost::shared_ptr<I3ExtraGeometryItemUnion> makeI3ExtraGeometryItemUnion(const bp::list& data)
    {
        std::vector<I3ExtraGeometryItemConstPtr> arg;
        
        for (unsigned int i = 0; i < bp::len(data); ++i)
        {
            arg.push_back(boost::python::extract<I3ExtraGeometryItemConstPtr>(data[i]));
        }
        
        return boost::shared_ptr<I3ExtraGeometryItemUnion>(new I3ExtraGeometryItemUnion(arg));
    }

};

using namespace I3ExtraGeometryItem_utils;


void register_I3ExtraGeometryItem()
{
    {
        bp::scope I3ExtraGeometryItem_scope = 
        bp::class_<I3ExtraGeometryItemWrapper, boost::shared_ptr<I3ExtraGeometryItemWrapper>, boost::noncopyable>
        ("I3ExtraGeometryItem",
         bp::init<>()
        )
        .def("DoesLineIntersect", bp::pure_virtual(&I3ExtraGeometryItem::DoesLineIntersect))
        .def("GetBoundingBox", bp::pure_virtual(&I3ExtraGeometryItem::GetBoundingBox))
        .def_pickle(bp::boost_serializable_pickle_suite<I3ExtraGeometryItem>())
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemWrapper>, shared_ptr<const I3ExtraGeometryItem> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemWrapper>, shared_ptr<I3ExtraGeometryItem> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemWrapper>, shared_ptr<const I3ExtraGeometryItemWrapper> >();
    utils::register_const_ptr<I3ExtraGeometryItem>();

    std_pair_to_python_converter<I3Position, I3Position>();
    
    class_<std::vector<I3ExtraGeometryItemConstPtr> >("I3ExtraGeometryItemConstPtrVect")
    .def(bp::list_indexing_suite<std::vector<I3ExtraGeometryItemConstPtr> >())
    ;

    
    // union
    {
        bp::class_<
        I3ExtraGeometryItemUnion, 
        boost::shared_ptr<I3ExtraGeometryItemUnion>, 
        bases<I3ExtraGeometryItem>
        >
        (
         "I3ExtraGeometryItemUnion",
         bp::init<
         const std::vector<I3ExtraGeometryItemConstPtr> &
         >(
           (
            bp::arg("elements")
           )
          )
        )
        .def(bp::init<>())
        .def("__init__", bp::make_constructor(makeI3ExtraGeometryItemUnion))
        .def(bp::dataclass_suite<I3ExtraGeometryItemUnion>())
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemUnion>, shared_ptr<const I3ExtraGeometryItemUnion> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemUnion>, shared_ptr<I3ExtraGeometryItem> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemUnion>, shared_ptr<const I3ExtraGeometryItem> >();
    utils::register_const_ptr<I3ExtraGeometryItemUnion>();

    
    // move
    {
        bp::class_<
        I3ExtraGeometryItemMove,
        boost::shared_ptr<I3ExtraGeometryItemMove>,
        bases<I3ExtraGeometryItem>
        >
        (
         "I3ExtraGeometryItemMove",
         bp::init<
         I3ExtraGeometryItemConstPtr, const I3Position &
         >(
           (
            bp::arg("element"),
            bp::arg("offset")
            )
           )
         )
        .def(bp::init<>())
        .def(bp::dataclass_suite<I3ExtraGeometryItemMove>())
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemMove>, shared_ptr<const I3ExtraGeometryItemMove> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemMove>, shared_ptr<I3ExtraGeometryItem> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemMove>, shared_ptr<const I3ExtraGeometryItem> >();
    utils::register_const_ptr<I3ExtraGeometryItemMove>();

    
    // cylinder
    {
        bp::class_<
        I3ExtraGeometryItemCylinder,
        boost::shared_ptr<I3ExtraGeometryItemCylinder>,
        bases<I3ExtraGeometryItem>
        >
        (
         "I3ExtraGeometryItemCylinder",
         bp::init<
         const I3Position &, const I3Position &, double
         >(
           (
            bp::arg("from"),
            bp::arg("to"),
            bp::arg("radius")
            )
           )
         )
        .def(bp::init<>())
        .def(bp::dataclass_suite<I3ExtraGeometryItemCylinder>())
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemCylinder>, shared_ptr<const I3ExtraGeometryItemCylinder> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemCylinder>, shared_ptr<I3ExtraGeometryItem> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemCylinder>, shared_ptr<const I3ExtraGeometryItem> >();
    utils::register_const_ptr<I3ExtraGeometryItemCylinder>();


}
