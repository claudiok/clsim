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
 * @file I3CLSimLightSource.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <icetray/I3Units.h>

#include <clsim/I3CLSimLightSource.h>
#include <boost/preprocessor/seq.hpp>
     
#include <icetray/python/list_indexing_suite.hpp>
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
        .def(bp::init<const I3CLSimFlasherPulse &>
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

        .add_property("flasherPulse", make_function(&I3CLSimLightSource::GetFlasherPulse, bp::return_value_policy<bp::copy_const_reference>()))
        .def("GetFlasherPulse", &I3CLSimLightSource::GetFlasherPulse, bp::return_value_policy<bp::copy_const_reference>())
        ;
        
        
        bp::enum_<I3CLSimLightSource::LightSourceType>("LightSourceType")
        BOOST_PP_SEQ_FOR_EACH(ENUM_DEF,I3CLSimLightSource,I3CLSIMLIGHTSOURCE_H_I3CLSimLightSource_LightSourceType)
        .export_values()
        ;

    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSource>, shared_ptr<const I3CLSimLightSource> >();
}
