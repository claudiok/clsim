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
 * @file I3CLSimMediumProperties.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/I3CLSimMediumProperties.h>

#include <boost/preprocessor/seq.hpp>

#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;
namespace bp = boost::python;

void register_I3CLSimMediumProperties()
{
    {
        bp::scope I3CLSimMediumProperties_scope = 
        bp::class_<I3CLSimMediumProperties, bases<I3FrameObject>, boost::shared_ptr<I3CLSimMediumProperties> >
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
        .def("GetDirectionalAbsorptionLengthCorrection", &I3CLSimMediumProperties::GetDirectionalAbsorptionLengthCorrection)
        .def("GetPreScatterDirectionTransform", &I3CLSimMediumProperties::GetPreScatterDirectionTransform)
        .def("GetPostScatterDirectionTransform", &I3CLSimMediumProperties::GetPostScatterDirectionTransform)
        .def("GetIceTiltZShift", &I3CLSimMediumProperties::GetIceTiltZShift)
        .def("SetAbsorptionLength", &I3CLSimMediumProperties::SetAbsorptionLength)
        .def("SetScatteringLength", &I3CLSimMediumProperties::SetScatteringLength)
        .def("SetPhaseRefractiveIndex", &I3CLSimMediumProperties::SetPhaseRefractiveIndex)
        .def("SetGroupRefractiveIndexOverride", &I3CLSimMediumProperties::SetGroupRefractiveIndexOverride)
        .def("SetScatteringCosAngleDistribution", &I3CLSimMediumProperties::SetScatteringCosAngleDistribution)
        .def("SetDirectionalAbsorptionLengthCorrection", &I3CLSimMediumProperties::SetDirectionalAbsorptionLengthCorrection)
        .def("SetPreScatterDirectionTransform", &I3CLSimMediumProperties::SetPreScatterDirectionTransform)
        .def("SetPostScatterDirectionTransform", &I3CLSimMediumProperties::SetPostScatterDirectionTransform)
        .def("SetIceTiltZShift", &I3CLSimMediumProperties::SetIceTiltZShift)
        .add_property("AbsorptionLength", &I3CLSimMediumProperties::GetAbsorptionLength, &I3CLSimMediumProperties::SetAbsorptionLength)
        .add_property("ScatteringLength", &I3CLSimMediumProperties::GetScatteringLength, &I3CLSimMediumProperties::SetScatteringLength)
        .add_property("PhaseRefractiveIndex", &I3CLSimMediumProperties::GetPhaseRefractiveIndex, &I3CLSimMediumProperties::SetPhaseRefractiveIndex)
        .add_property("GroupRefractiveIndexOverride", &I3CLSimMediumProperties::GetGroupRefractiveIndexOverride, &I3CLSimMediumProperties::SetGroupRefractiveIndexOverride)
        .add_property("ScatteringCosAngleDistribution", &I3CLSimMediumProperties::GetScatteringCosAngleDistribution, &I3CLSimMediumProperties::SetScatteringCosAngleDistribution)
        .add_property("DirectionalAbsorptionLengthCorrection", &I3CLSimMediumProperties::GetDirectionalAbsorptionLengthCorrection, &I3CLSimMediumProperties::SetDirectionalAbsorptionLengthCorrection)
        .add_property("PreScatterDirectionTransform", &I3CLSimMediumProperties::GetPreScatterDirectionTransform, &I3CLSimMediumProperties::SetPreScatterDirectionTransform)
        .add_property("PostScatterDirectionTransform", &I3CLSimMediumProperties::GetPostScatterDirectionTransform, &I3CLSimMediumProperties::SetPostScatterDirectionTransform)
        .add_property("IceTiltZShift", &I3CLSimMediumProperties::GetIceTiltZShift, &I3CLSimMediumProperties::SetIceTiltZShift)

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

        .def(dataclass_suite<I3CLSimMediumProperties>())
        ;
    }
    
    register_pointer_conversions<I3CLSimMediumProperties>();
}
