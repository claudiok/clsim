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

#include <clsim/I3CLSimPhotonSeriesToPhotonSeriesMapConverter.h>

using namespace boost::python;
namespace bp = boost::python;

void register_I3CLSimPhotonSeriesToPhotonSeriesMapConverter()
{
    {
        bp::scope I3CLSimPhotonSeriesToPhotonSeriesMapConverter_scope = 
        bp::class_<I3CLSimPhotonSeriesToPhotonSeriesMapConverter, boost::shared_ptr<I3CLSimPhotonSeriesToPhotonSeriesMapConverter>, boost::noncopyable>("I3CLSimPhotonSeriesToPhotonSeriesMapConverter")
        .def("SetPhotonSeries", &I3CLSimPhotonSeriesToPhotonSeriesMapConverter::SetPhotonSeries)
        .def("SetGeometry", &I3CLSimPhotonSeriesToPhotonSeriesMapConverter::SetGeometry)
        .def("GetPhotonSeriesMapForIdentifier", &I3CLSimPhotonSeriesToPhotonSeriesMapConverter::GetPhotonSeriesMapForIdentifier)
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimPhotonSeriesToPhotonSeriesMapConverter>, shared_ptr<const I3CLSimPhotonSeriesToPhotonSeriesMapConverter> >();

}
