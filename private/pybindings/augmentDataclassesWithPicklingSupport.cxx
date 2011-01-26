//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of IceTray.
//
//   IceTray is free software; you can redistribute it and/or modify
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

#include "boost_serializable_pickle_suite.hpp"

#include <boost/preprocessor.hpp>

#include <icetray/I3Frame.h>

#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/I3Direction.h>
#include <dataclasses/I3Position.h>

namespace bp=boost::python;

void register_augmentDataclassesWithPicklingSupport()
{
    // this is a bit of hack to add pickling support to some of the dataclasses
    // without having to modify the dataclasses project.
    // If we want pickling support for dataclasses, this support should
    // be integrated in there directy. -ck

    bp::def("I3Frame___getstate__", &boost_serializable_pickle_suite<I3Frame>::getstate);
    bp::def("I3Frame___setstate__", &boost_serializable_pickle_suite<I3Frame>::setstate);

    bp::def("I3Particle___getstate__", &boost_serializable_pickle_suite<I3Particle>::getstate);
    bp::def("I3Particle___setstate__", &boost_serializable_pickle_suite<I3Particle>::setstate);

    bp::def("I3Direction___getstate__", &boost_serializable_pickle_suite<I3Direction>::getstate);
    bp::def("I3Direction___setstate__", &boost_serializable_pickle_suite<I3Direction>::setstate);

    bp::def("I3Position___getstate__", &boost_serializable_pickle_suite<I3Position>::getstate);
    bp::def("I3Position___setstate__", &boost_serializable_pickle_suite<I3Position>::setstate);
}
