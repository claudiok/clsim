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

#include <clsim/shadow/I3ExtraGeometryItem.h>
#include <clsim/shadow/I3ExtraGeometryItemUnion.h>
#include <clsim/shadow/I3ExtraGeometryItemMove.h>

#include <icetray/python/std_vector_indexing_suite.hpp>
#include <icetray/python/dataclass_suite.hpp>

#include <boost/preprocessor/seq.hpp>
#include "const_ptr_helpers.h"


using namespace boost::python;
namespace bp = boost::python;

struct I3ExtraGeometryItemWrapper : I3ExtraGeometryItem, bp::wrapper<I3ExtraGeometryItem>
{
    // pure virtual
    virtual bool DoesLineIntersect(const I3Position &lineStart,
                                   const I3Position &lineEnd) const
    {
        return this->get_override("DoesLineIntersect")(lineStart, lineEnd);
    }

    virtual std::pair<I3Position, I3Position> GetBoundingBox() const
    {
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
        .def(freeze())
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemWrapper>, shared_ptr<const I3ExtraGeometryItem> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemWrapper>, shared_ptr<I3ExtraGeometryItem> >();
    bp::implicitly_convertible<shared_ptr<I3ExtraGeometryItemWrapper>, shared_ptr<const I3ExtraGeometryItemWrapper> >();
    utils::register_const_ptr<I3ExtraGeometryItem>();

    std_pair_to_python_converter<I3Position, I3Position>();
    
    class_<std::vector<I3ExtraGeometryItemConstPtr> >("I3ExtraGeometryItemConstPtrVect")
    .def(bp::std_vector_indexing_suite<std::vector<I3ExtraGeometryItemConstPtr> >())
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


}
