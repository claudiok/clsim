//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
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

#include <clsim/I3CLSimLightSourceParameterization.h>
#include "clsim/I3CLSimLightSourceToStepConverter.h"

#include <icetray/python/std_vector_indexing_suite.hpp>
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
    .def(std_vector_indexing_suite<I3CLSimLightSourceParameterizationSeries>())
    ;
    
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceParameterization>, shared_ptr<const I3CLSimLightSourceParameterization> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimLightSourceParameterizationSeries>, shared_ptr<const I3CLSimLightSourceParameterizationSeries> >();

    from_python_sequence<I3CLSimLightSourceParameterizationSeries, variable_capacity_policy>();
    
}
