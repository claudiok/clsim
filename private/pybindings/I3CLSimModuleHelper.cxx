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

#include <clsim/I3CLSimModuleHelper.h>
#include <clsim/I3CLSimRandomValueInterpolatedDistribution.h>
#include <clsim/I3CLSimRandomValueWlenCherenkovNoDispersion.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

#include <boost/foreach.hpp>

using namespace boost::python;
namespace bp = boost::python;


void register_I3ModuleHelper()
{
    // this can be used for testing purposes
    bp::def("makeCherenkovWavelengthGenerator", &I3CLSimModuleHelper::makeCherenkovWavelengthGenerator);
    bp::def("makeWavelengthGenerator", &I3CLSimModuleHelper::makeWavelengthGenerator);
    
    
}
