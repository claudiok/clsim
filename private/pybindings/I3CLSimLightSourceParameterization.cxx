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
 * @file I3CLSimLightSourceParameterization.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/I3CLSimLightSourceParameterization.h>
#include "clsim/I3CLSimLightSourceToStepConverter.h"

#include <icetray/python/list_indexing_suite.hpp>
#include <icetray/python/copy_suite.hpp>

using namespace boost::python;
namespace bp = boost::python;

void register_I3CLSimLightSourceParameterization()
{
    {
        bp::scope I3CLSimLightSourceParameterization_scope = 
        bp::class_<I3CLSimLightSourceParameterization, boost::shared_ptr<I3CLSimLightSourceParameterization> >
        ("I3CLSimLightSourceParameterization",
         bp::init<I3CLSimLightSourceToStepConverterPtr, I3Particle::ParticleType, double, double, bool>
         (
          (
           bp::arg("converter"),
           bp::arg("forParticleType"),
           bp::arg("fromEnergy"),
           bp::arg("toEnergy"),
           bp::arg("needsLength")=false
          )
         )
        )
        .def(init<I3CLSimLightSourceToStepConverterPtr, const I3Particle &, double, double, bool>
             (
              (
               bp::arg("converter"),
               bp::arg("forParticleType"),
               bp::arg("fromEnergy"),
               bp::arg("toEnergy"),
               bp::arg("needsLength")=false
              )
             )
            ) 
        .def(init<I3CLSimLightSourceToStepConverterPtr, const I3CLSimLightSourceParameterization::AllParticles_t &, double, double, bool>
             (
              (
               bp::arg("converter"),
               bp::arg("forParticleType"),
               bp::arg("fromEnergy"),
               bp::arg("toEnergy"),
               bp::arg("needsLength")=false
              )
             )
            ) 
        .def(init<I3CLSimLightSourceToStepConverterPtr, I3CLSimFlasherPulse::FlasherPulseType>
             (
              (
               bp::arg("converter"),
               bp::arg("forFlasherPulseType")
               )
              )
             ) 
        .def(init<>()) // the class also has a default constructor
        .def(copy_suite<I3CLSimLightSourceParameterization>())
        .def_readwrite("converter", &I3CLSimLightSourceParameterization::converter)
#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        .def_readwrite("forParticleType", &I3CLSimLightSourceParameterization::forParticleType)
#else
        .def_readwrite("forPdgEncoding", &I3CLSimLightSourceParameterization::forPdgEncoding)
#endif
        .def_readwrite("fromEnergy", &I3CLSimLightSourceParameterization::fromEnergy)
        .def_readwrite("toEnergy", &I3CLSimLightSourceParameterization::toEnergy)
        .def_readwrite("needsLength", &I3CLSimLightSourceParameterization::needsLength)
        .def_readwrite("catchAll", &I3CLSimLightSourceParameterization::catchAll)

        .def_readwrite("flasherMode", &I3CLSimLightSourceParameterization::flasherMode)
        .def_readwrite("forFlasherPulseType", &I3CLSimLightSourceParameterization::forFlasherPulseType)

        .def("IsValidForParticle", &I3CLSimLightSourceParameterization::IsValidForParticle, bp::arg("particle"))
        .def("IsValid", &I3CLSimLightSourceParameterization::IsValid, bp::arg("type"), bp::arg("energy"), bp::arg("length")=NAN)
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        .def("IsValidForPdgEncoding", &I3CLSimLightSourceParameterization::IsValidForPdgEncoding, bp::arg("encoding"), bp::arg("energy"), bp::arg("length")=NAN)
#endif
        .def("IsValidForLightSource", &I3CLSimLightSourceParameterization::IsValidForLightSource, bp::arg("lightSource"))
        ;

        bp::class_<I3CLSimLightSourceParameterization::AllParticles_t>("AllParticles_t");
    }

    class_<I3CLSimLightSourceParameterizationSeries, I3CLSimLightSourceParameterizationSeriesPtr>("I3CLSimLightSourceParameterizationSeries")
    .def(list_indexing_suite<I3CLSimLightSourceParameterizationSeries>())
    ;
    
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceParameterization>, shared_ptr<const I3CLSimLightSourceParameterization> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceParameterizationSeries>, shared_ptr<const I3CLSimLightSourceParameterizationSeries> >();

    from_python_sequence<I3CLSimLightSourceParameterizationSeries, variable_capacity_policy>();
    
}
