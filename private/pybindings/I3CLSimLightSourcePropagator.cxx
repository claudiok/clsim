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
 * @file I3CLSimLightSourceToStepConverter.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <clsim/I3CLSimLightSourcePropagator.h>
#include <clsim/I3CLSimLightSourcePropagatorFromI3PropagatorService.h>
#include <icetray/python/list_indexing_suite.hpp>

I3_POINTER_TYPEDEFS(I3CLSimLightSourcePropagator);
typedef std::vector<I3CLSimLightSourcePropagatorPtr> I3CLSimLightSourcePropagatorSeries;
I3_POINTER_TYPEDEFS(I3CLSimLightSourcePropagatorSeries);

#include "clsim/I3CLSimStep.h"

namespace bp = boost::python;

namespace {

void Convert(I3CLSimLightSourcePropagator &prop, I3CLSimLightSourceConstPtr source, uint32_t identifier, bp::object secondary_callback, bp::object step_callback)
{
    I3CLSimLightSourcePropagator::secondary_callback emitSecondary = [&](I3CLSimLightSourceConstPtr &lightSource, uint32_t lightSourceIdentifier) ->bool
    {
        bp::object ret = secondary_callback(boost::const_pointer_cast<I3CLSimLightSource>(lightSource), lightSourceIdentifier);
        return bp::extract<bool>(ret);
    };
    I3CLSimLightSourcePropagator::step_callback emitStep = [&](const I3CLSimStep &step)
    {
        step_callback(step);
    };
    
    prop.Convert(source, identifier, emitSecondary, emitStep);
}

}

void register_I3CLSimLightSourcePropagator()
{
    bp::class_<I3CLSimLightSourcePropagator, boost::shared_ptr<I3CLSimLightSourcePropagator>, boost::noncopyable>("I3CLSimLightSourcePropagator", bp::no_init)
    .def("SetRandomService", bp::pure_virtual(&I3CLSimLightSourcePropagator::SetRandomService))
    .def("SetWlenBias", bp::pure_virtual(&I3CLSimLightSourcePropagator::SetWlenBias))
    .def("SetMediumProperties", bp::pure_virtual(&I3CLSimLightSourcePropagator::SetMediumProperties))
    .def("Initialize", bp::pure_virtual(&I3CLSimLightSourcePropagator::Initialize))
    .def("IsInitialized", bp::pure_virtual(&I3CLSimLightSourcePropagator::IsInitialized))
    .def("Convert", &Convert)
    ;
    
    bp::class_<I3CLSimLightSourcePropagatorSeries, I3CLSimLightSourcePropagatorSeriesPtr>("I3CLSimLightSourcePropagatorSeries")
    .def(bp::list_indexing_suite<I3CLSimLightSourcePropagatorSeries>())
    ;
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourcePropagator>, boost::shared_ptr<const I3CLSimLightSourcePropagator> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimLightSourcePropagatorSeries>, boost::shared_ptr<const I3CLSimLightSourcePropagatorSeries> >();

    from_python_sequence<I3CLSimLightSourcePropagatorSeries, variable_capacity_policy>();
    
    bp::class_<I3CLSimLightSourcePropagatorFromI3PropagatorService,
               boost::shared_ptr<I3CLSimLightSourcePropagatorFromI3PropagatorService>,
               bp::bases<I3CLSimLightSourcePropagator>,
               boost::noncopyable>("I3CLSimLightSourcePropagatorFromI3PropagatorService",
                                   bp::init<I3ParticleTypePropagatorServiceMapPtr>())
    ;
}
