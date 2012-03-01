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

#include <clsim/I3CLSimEventStatistics.h>

#include <boost/preprocessor/seq.hpp>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;


void register_I3CLSimEventStatistics()
{
    uint64_t (I3CLSimEventStatistics::* GetNumberOfPhotonsGeneratedForParticle_oneary)(const I3Particle &) const = &I3CLSimEventStatistics::GetNumberOfPhotonsGeneratedForParticle;
    uint64_t (I3CLSimEventStatistics::* GetNumberOfPhotonsGeneratedForParticle_twoary)(uint64_t, int) const = &I3CLSimEventStatistics::GetNumberOfPhotonsGeneratedForParticle;
    double (I3CLSimEventStatistics::* GetSumOfWeightsPhotonsGeneratedForParticle_oneary)(const I3Particle &) const = &I3CLSimEventStatistics::GetSumOfWeightsPhotonsGeneratedForParticle;
    double (I3CLSimEventStatistics::* GetSumOfWeightsPhotonsGeneratedForParticle_twoary)(uint64_t, int) const = &I3CLSimEventStatistics::GetSumOfWeightsPhotonsGeneratedForParticle;

    uint64_t (I3CLSimEventStatistics::* GetNumberOfPhotonsAtDOMsForParticle_oneary)(const I3Particle &) const = &I3CLSimEventStatistics::GetNumberOfPhotonsAtDOMsForParticle;
    uint64_t (I3CLSimEventStatistics::* GetNumberOfPhotonsAtDOMsForParticle_twoary)(uint64_t, int) const = &I3CLSimEventStatistics::GetNumberOfPhotonsAtDOMsForParticle;
    double (I3CLSimEventStatistics::* GetSumOfWeightsPhotonsAtDOMsForParticle_oneary)(const I3Particle &) const = &I3CLSimEventStatistics::GetSumOfWeightsPhotonsAtDOMsForParticle;
    double (I3CLSimEventStatistics::* GetSumOfWeightsPhotonsAtDOMsForParticle_twoary)(uint64_t, int) const = &I3CLSimEventStatistics::GetSumOfWeightsPhotonsAtDOMsForParticle;

    
    void (I3CLSimEventStatistics::* AddNumPhotonsGeneratedWithWeights_oneary)(uint64_t, double, const I3Particle &particle) = &I3CLSimEventStatistics::AddNumPhotonsGeneratedWithWeights;
    void (I3CLSimEventStatistics::* AddNumPhotonsGeneratedWithWeights_twoary)(uint64_t, double, uint64_t, int) = &I3CLSimEventStatistics::AddNumPhotonsGeneratedWithWeights;

    void (I3CLSimEventStatistics::* AddNumPhotonsAtDOMsWithWeights_oneary)(uint64_t, double, const I3Particle &particle) = &I3CLSimEventStatistics::AddNumPhotonsAtDOMsWithWeights;
    void (I3CLSimEventStatistics::* AddNumPhotonsAtDOMsWithWeights_twoary)(uint64_t, double, uint64_t, int) = &I3CLSimEventStatistics::AddNumPhotonsAtDOMsWithWeights;


    
    {
        scope I3CLSimEventStatistics_scope = 
        class_<I3CLSimEventStatistics, bases<I3FrameObject>, boost::shared_ptr<I3CLSimEventStatistics> >("I3CLSimEventStatistics")

        .def("GetNumberOfPhotonsGeneratedForParticle", GetNumberOfPhotonsGeneratedForParticle_oneary, bp::arg("particle"))
        .def("GetNumberOfPhotonsGeneratedForParticle", GetNumberOfPhotonsGeneratedForParticle_twoary, bp::arg("majorID"), bp::arg("minorID"))
        .def("GetTotalNumberOfPhotonsGenerated", &I3CLSimEventStatistics::GetTotalNumberOfPhotonsGenerated)

        .def("GetSumOfWeightsPhotonsGeneratedForParticle", GetSumOfWeightsPhotonsGeneratedForParticle_oneary, bp::arg("particle"))
        .def("GetSumOfWeightsPhotonsGeneratedForParticle", GetSumOfWeightsPhotonsGeneratedForParticle_twoary, bp::arg("majorID"), bp::arg("minorID"))
        .def("GetTotalSumOfWeightsPhotonsGenerated", &I3CLSimEventStatistics::GetTotalSumOfWeightsPhotonsGenerated)

        
        .def("GetNumberOfPhotonsAtDOMsForParticle", GetNumberOfPhotonsAtDOMsForParticle_oneary, bp::arg("particle"))
        .def("GetNumberOfPhotonsAtDOMsForParticle", GetNumberOfPhotonsAtDOMsForParticle_twoary, bp::arg("majorID"), bp::arg("minorID"))
        .def("GetTotalNumberOfPhotonsAtDOMs", &I3CLSimEventStatistics::GetTotalNumberOfPhotonsAtDOMs)
        
        .def("GetSumOfWeightsPhotonsAtDOMsForParticle", GetSumOfWeightsPhotonsAtDOMsForParticle_oneary, bp::arg("particle"))
        .def("GetSumOfWeightsPhotonsAtDOMsForParticle", GetSumOfWeightsPhotonsAtDOMsForParticle_twoary, bp::arg("majorID"), bp::arg("minorID"))
        .def("GetTotalSumOfWeightsPhotonsAtDOMs", &I3CLSimEventStatistics::GetTotalSumOfWeightsPhotonsAtDOMs)

        
        .def("AddNumPhotonsGeneratedWithWeights", AddNumPhotonsGeneratedWithWeights_oneary, bp::arg("numPhotons"), bp::arg("weightsForPhotons"), bp::arg("particle"))
        .def("AddNumPhotonsGeneratedWithWeights", AddNumPhotonsGeneratedWithWeights_twoary, bp::args("numPhotons", "weightsForPhotons", "majorID", "minorID"))

        .def("AddNumPhotonsAtDOMsWithWeights", AddNumPhotonsAtDOMsWithWeights_oneary, bp::arg("numPhotons"), bp::arg("weightsForPhotons"), bp::arg("particle"))
        .def("AddNumPhotonsAtDOMsWithWeights", AddNumPhotonsAtDOMsWithWeights_twoary, bp::args("numPhotons", "weightsForPhotons", "majorID", "minorID"))

        .def(dataclass_suite<I3CLSimEventStatistics>())
        ;
    }

    register_pointer_conversions<I3CLSimEventStatistics>();
}
