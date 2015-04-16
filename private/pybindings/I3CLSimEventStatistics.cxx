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
 * @file I3CLSimEventStatistics.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

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
