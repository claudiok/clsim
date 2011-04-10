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

#include <clsim/I3CLSimParticleParameterization.h>
#include "clsim/I3CLSimParticleToStepConverter.h"

#include <icetray/python/std_vector_indexing_suite.hpp>
#include <icetray/python/copy_suite.hpp>

using namespace boost::python;
namespace bp = boost::python;

void register_I3CLSimParticleParameterization()
{
    {
        bp::scope I3CLSimParticleParameterization_scope = 
        bp::class_<I3CLSimParticleParameterization, boost::shared_ptr<I3CLSimParticleParameterization> >
        ("I3CLSimParticleParameterization",
         bp::init<I3CLSimParticleToStepConverterPtr, I3Particle::ParticleType, double, double, bool>
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
        .def(init<>()) // the class also has a defualt constructor
        .def(copy_suite<I3CLSimParticleParameterization>())
        .def_readwrite("converter", &I3CLSimParticleParameterization::converter)
        .def_readwrite("forParticleType", &I3CLSimParticleParameterization::forParticleType)
        .def_readwrite("fromEnergy", &I3CLSimParticleParameterization::fromEnergy)
        .def_readwrite("toEnergy", &I3CLSimParticleParameterization::toEnergy)
        .def_readwrite("needsLength", &I3CLSimParticleParameterization::needsLength)
        .def("IsValidForParticle", &I3CLSimParticleParameterization::IsValidForParticle, bp::arg("particle"))
        .def("IsValid", &I3CLSimParticleParameterization::IsValid, bp::arg("type"), bp::arg("energy"), bp::arg("length")=NAN)
        ;
    }
    

    class_<I3CLSimParticleParameterizationSeries, I3CLSimParticleParameterizationSeriesPtr>("I3CLSimParticleParameterizationSeries")
    .def(std_vector_indexing_suite<I3CLSimParticleParameterizationSeries>())
    ;
    
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleParameterization>, shared_ptr<const I3CLSimParticleParameterization> >();
    bp::implicitly_convertible<shared_ptr<I3CLSimParticleParameterizationSeries>, shared_ptr<const I3CLSimParticleParameterizationSeries> >();

    from_python_sequence<I3CLSimParticleParameterizationSeries, variable_capacity_policy>();
    
}
