//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is free software; you can redistribute it and/or modify
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

#include <icetray/I3Units.h>

#include <clsim/I3CLSimLightSource.h>
#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>
#include <icetray/python/std_map_indexing_suite.hpp>

namespace bp = boost::python;

#define ENUM_DEF(r,data,T) .value(BOOST_PP_STRINGIZE(T), data::T)


void register_I3CLSimLightSource()
{
    {
        bp::scope i3clsimlightsource_scope = 
        bp::class_<I3CLSimLightSource, boost::shared_ptr<I3CLSimLightSource> >
        ("I3CLSimLightSource",
         bp::init<const I3Particle &>
         (
          (
           bp::arg("particle")
          )
         )
        )
        .def(bp::init<const I3FlasherInfo &>
             (
              (
               bp::arg("flasher")
              )
             )
            )
        
        .add_property("type", &I3CLSimLightSource::GetType)
        .def("GetType", &I3CLSimLightSource::GetType)
        
        .add_property("particle", make_function(&I3CLSimLightSource::GetParticle, bp::return_value_policy<bp::copy_const_reference>()))
        .def("GetParticle", &I3CLSimLightSource::GetParticle, bp::return_value_policy<bp::copy_const_reference>())

        .add_property("flasher_info", make_function(&I3CLSimLightSource::GetFlasherInfo, bp::return_value_policy<bp::copy_const_reference>()))
        .def("GetFlasherInfo", &I3CLSimLightSource::GetFlasherInfo, bp::return_value_policy<bp::copy_const_reference>())
        ;
        
        
        bp::enum_<I3CLSimLightSource::LightSourceType>("LightSourceType")
        BOOST_PP_SEQ_FOR_EACH(ENUM_DEF,I3CLSimLightSource,I3CLSIMLIGHTSOURCE_H_I3CLSimLightSource_LightSourceType)
        .export_values()
        ;

    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSource>, shared_ptr<const I3CLSimLightSource> >();
}
