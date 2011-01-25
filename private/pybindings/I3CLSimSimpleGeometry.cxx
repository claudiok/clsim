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

#include <clsim/I3CLSimSimpleGeometry.h>
#include <clsim/I3CLSimSimpleGeometryUserConfigurable.h>
#include <clsim/I3CLSimSimpleGeometryTextFile.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimSimpleGeometryWrapper : I3CLSimSimpleGeometry, bp::wrapper<I3CLSimSimpleGeometry>
{
    // pure virtual
    virtual std::size_t size() const {return this->get_override("size")();}
    virtual double GetOMRadius() const {return this->get_override("GetOMRadius")();}
    
    virtual const std::vector<int32_t> &GetStringIDVector() const {return this->get_override("GetStringIDVector")();}
    virtual const std::vector<uint32_t> &GetDomIDVector() const {return this->get_override("GetDomIDVector")();}
    virtual const std::vector<double> &GetPosXVector() const {return this->get_override("GetPosXVector")();}
    virtual const std::vector<double> &GetPosYVector() const {return this->get_override("GetPosYVector")();}
    virtual const std::vector<double> &GetPosZVector() const {return this->get_override("GetPosZVector")();}
    
    virtual int32_t GetStringID(std::size_t pos) const {return this->get_override("GetStringID")();}
    virtual uint32_t GetDomID(std::size_t pos) const {return this->get_override("GetDomID")();}
    virtual double GetPosX(std::size_t pos) const {return this->get_override("GetPosX")();}
    virtual double GetPosY(std::size_t pos) const {return this->get_override("GetPosY")();}
    virtual double GetPosZ(std::size_t pos) const {return this->get_override("GetPosZ")();}
};

void register_I3CLSimSimpleGeometry()
{
    {
        bp::scope I3CLSimSimpleGeometry_scope = 
        bp::class_<I3CLSimSimpleGeometryWrapper, boost::shared_ptr<I3CLSimSimpleGeometryWrapper>, boost::noncopyable>("I3CLSimSimpleGeometry", bp::no_init)
        .def("size", bp::pure_virtual(&I3CLSimSimpleGeometry::size))
        .def("__len__", bp::pure_virtual(&I3CLSimSimpleGeometry::size))
        
        .def("GetOMRadius", bp::pure_virtual(&I3CLSimSimpleGeometry::GetOMRadius))

        .def("GetStringIDVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetStringIDVector), bp::return_value_policy<bp::copy_const_reference>())
        .def("GetDomIDVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetDomIDVector), bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosXVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosXVector), bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosYVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosYVector), bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosZVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosZVector), bp::return_value_policy<bp::copy_const_reference>())

        .add_property("stringIDs", bp::make_function(&I3CLSimSimpleGeometry::GetStringIDVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("domIDs", bp::make_function(&I3CLSimSimpleGeometry::GetDomIDVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("posX", bp::make_function(&I3CLSimSimpleGeometry::GetPosXVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("posY", bp::make_function(&I3CLSimSimpleGeometry::GetPosYVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("posZ", bp::make_function(&I3CLSimSimpleGeometry::GetPosZVector, bp::return_value_policy<bp::copy_const_reference>()))

        .def("GetDomIDVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetDomIDVector), bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosXVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosXVector), bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosYVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosYVector), bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosZVector", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosZVector), bp::return_value_policy<bp::copy_const_reference>())
        
        .def("GetStringID", bp::pure_virtual(&I3CLSimSimpleGeometry::GetStringID))
        .def("GetDomID", bp::pure_virtual(&I3CLSimSimpleGeometry::GetDomID))
        .def("GetPosX", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosX))
        .def("GetPosY", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosY))
        .def("GetPosZ", bp::pure_virtual(&I3CLSimSimpleGeometry::GetPosZ))
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryWrapper>, shared_ptr<const I3CLSimSimpleGeometry> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryWrapper>, shared_ptr<I3CLSimSimpleGeometry> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryWrapper>, shared_ptr<const I3CLSimSimpleGeometryWrapper> >();
    
    // I3CLSimSimpleGeometryUserConfigurable
    {
        bp::class_<
        I3CLSimSimpleGeometryUserConfigurable, 
        boost::shared_ptr<I3CLSimSimpleGeometryUserConfigurable>, 
        bases<I3CLSimSimpleGeometry>,
        boost::noncopyable
        >
        (
         "I3CLSimSimpleGeometryUserConfigurable",
         bp::init<
         double, std::size_t
         >(
           (
            bp::arg("OMRadius"),
            bp::arg("numOMs")
           )
          )
        )
        .def("SetStringID", &I3CLSimSimpleGeometryUserConfigurable::SetStringID)
        .def("SetDomID", &I3CLSimSimpleGeometryUserConfigurable::SetDomID)
        .def("SetPosX", &I3CLSimSimpleGeometryUserConfigurable::SetPosX)
        .def("SetPosY", &I3CLSimSimpleGeometryUserConfigurable::SetPosY)
        .def("SetPosZ", &I3CLSimSimpleGeometryUserConfigurable::SetPosZ)
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryUserConfigurable>, shared_ptr<const I3CLSimSimpleGeometryUserConfigurable> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryUserConfigurable>, shared_ptr<I3CLSimSimpleGeometry> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryUserConfigurable>, shared_ptr<const I3CLSimSimpleGeometry> >();

    
    // I3CLSimSimpleGeometryTextFile
    {
        bp::class_<
        I3CLSimSimpleGeometryTextFile, 
        boost::shared_ptr<I3CLSimSimpleGeometryTextFile>, 
        bases<I3CLSimSimpleGeometry>,
        boost::noncopyable
        >
        (
         "I3CLSimSimpleGeometryTextFile",
         bp::init<
         double, const std::string &,
         int32_t, int32_t,
         uint32_t, uint32_t
         >(
           (
            bp::arg("OMRadius"),
            bp::arg("filename"),
            bp::arg("ignoreStringIDsSmallerThan")=I3CLSimSimpleGeometryTextFile::default_ignoreStringIDsSmallerThan,
            bp::arg("ignoreStringIDsLargerThan")=I3CLSimSimpleGeometryTextFile::default_ignoreStringIDsLargerThan,
            bp::arg("ignoreDomIDsSmallerThan")=I3CLSimSimpleGeometryTextFile::default_ignoreDomIDsSmallerThan,
            bp::arg("ignoreDomIDsLargerThan")=I3CLSimSimpleGeometryTextFile::default_ignoreDomIDsLargerThan
            )
           )
         )
        ;
    }
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryTextFile>, shared_ptr<const I3CLSimSimpleGeometryTextFile> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryTextFile>, shared_ptr<I3CLSimSimpleGeometry> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimSimpleGeometryTextFile>, shared_ptr<const I3CLSimSimpleGeometry> >();
    
}
