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

#include <clsim/I3CLSimMediumProperties.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

// disabled for now until this gets released into offline-software
//#include "boost_serializable_pickle_suite.hpp"

using namespace boost::python;
namespace bp = boost::python;

void register_I3CLSimMediumProperties()
{
    {
        bp::scope I3CLSimMediumProperties_scope = 
        bp::class_<I3CLSimMediumProperties, boost::shared_ptr<I3CLSimMediumProperties>, boost::noncopyable>
        ("I3CLSimMediumProperties",
         bp::init<double, uint32_t, double, double, double, double>
         (
          (
           bp::arg("mediumDensity")=I3CLSimMediumProperties::default_mediumDensity,
           bp::arg("layersNum")=I3CLSimMediumProperties::default_layersNum,
           bp::arg("layersZStart")=I3CLSimMediumProperties::default_layersZStart,
           bp::arg("layersHeight")=I3CLSimMediumProperties::default_layersHeight,
           bp::arg("rockZCoordinate")=I3CLSimMediumProperties::default_rockZCoordinate,
           bp::arg("airZCoordinate")=I3CLSimMediumProperties::default_airZCoordinate
           
          )
         )
        )
        .def("IsReady", &I3CLSimMediumProperties::IsReady)
        .def("GetAbsorptionLength", &I3CLSimMediumProperties::GetAbsorptionLength)
        .def("GetScatteringLength", &I3CLSimMediumProperties::GetScatteringLength)
        .def("GetPhaseRefractiveIndex", &I3CLSimMediumProperties::GetPhaseRefractiveIndex)
        .def("GetGroupRefractiveIndexOverride", &I3CLSimMediumProperties::GetGroupRefractiveIndexOverride)
        .def("GetScatteringCosAngleDistribution", &I3CLSimMediumProperties::GetScatteringCosAngleDistribution)
        .def("SetAbsorptionLength", &I3CLSimMediumProperties::SetAbsorptionLength)
        .def("SetScatteringLength", &I3CLSimMediumProperties::SetScatteringLength)
        .def("SetPhaseRefractiveIndex", &I3CLSimMediumProperties::SetPhaseRefractiveIndex)
        .def("SetGroupRefractiveIndexOverride", &I3CLSimMediumProperties::SetGroupRefractiveIndexOverride)
        .def("SetScatteringCosAngleDistribution", &I3CLSimMediumProperties::SetScatteringCosAngleDistribution)
        .add_property("AbsorptionLength", &I3CLSimMediumProperties::GetAbsorptionLength, &I3CLSimMediumProperties::SetAbsorptionLength)
        .add_property("ScatteringLength", &I3CLSimMediumProperties::GetScatteringLength, &I3CLSimMediumProperties::SetScatteringLength)
        .add_property("PhaseRefractiveIndex", &I3CLSimMediumProperties::GetPhaseRefractiveIndex, &I3CLSimMediumProperties::SetPhaseRefractiveIndex)
        .add_property("GroupRefractiveIndexOverride", &I3CLSimMediumProperties::GetGroupRefractiveIndexOverride, &I3CLSimMediumProperties::SetGroupRefractiveIndexOverride)
        .add_property("ScatteringCosAngleDistribution", &I3CLSimMediumProperties::GetScatteringCosAngleDistribution, &I3CLSimMediumProperties::SetScatteringCosAngleDistribution)

        .def("GetMinWavelength", &I3CLSimMediumProperties::GetMinWavelength)
        .def("GetMaxWavelength", &I3CLSimMediumProperties::GetMaxWavelength)
        .add_property("MinWavelength", &I3CLSimMediumProperties::GetMinWavelength)
        .add_property("MaxWavelength", &I3CLSimMediumProperties::GetMaxWavelength)

        .def("GetMediumDensity", &I3CLSimMediumProperties::GetMediumDensity)
        .def("GetLayersNum", &I3CLSimMediumProperties::GetLayersNum)
        .def("GetLayersZStart", &I3CLSimMediumProperties::GetLayersZStart)
        .def("GetLayersHeight", &I3CLSimMediumProperties::GetLayersHeight)
        .add_property("MediumDensity", &I3CLSimMediumProperties::GetMediumDensity)
        .add_property("LayersNum", &I3CLSimMediumProperties::GetLayersNum)
        .add_property("LayersZStart", &I3CLSimMediumProperties::GetLayersZStart)
        .add_property("LayersHeight", &I3CLSimMediumProperties::GetLayersHeight)

        .def("GetRockZCoord", &I3CLSimMediumProperties::GetRockZCoord)
        .def("GetAirZCoord", &I3CLSimMediumProperties::GetAirZCoord)
        .add_property("RockZCoord", &I3CLSimMediumProperties::GetRockZCoord)
        .add_property("AirZCoord", &I3CLSimMediumProperties::GetAirZCoord)

        .def("GetForcedMinWlen", &I3CLSimMediumProperties::GetForcedMinWlen)
        .def("GetForcedMaxWlen", &I3CLSimMediumProperties::GetForcedMaxWlen)
        .def("SetForcedMinWlen", &I3CLSimMediumProperties::SetForcedMinWlen)
        .def("SetForcedMaxWlen", &I3CLSimMediumProperties::SetForcedMaxWlen)
        .add_property("ForcedMinWlen", &I3CLSimMediumProperties::GetForcedMinWlen, &I3CLSimMediumProperties::SetForcedMinWlen)
        .add_property("ForcedMaxWlen", &I3CLSimMediumProperties::GetForcedMaxWlen, &I3CLSimMediumProperties::SetForcedMaxWlen)

        //.def_pickle(boost_serializable_pickle_suite<I3CLSimMediumProperties>())
        ;
    }
    
    bp::implicitly_convertible<shared_ptr<I3CLSimMediumProperties>, shared_ptr<const I3CLSimMediumProperties> >();
    
    
}
