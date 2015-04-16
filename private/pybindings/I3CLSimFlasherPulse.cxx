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
 * @file I3CLSimFlasherPulse.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <icetray/I3Units.h>

#include <clsim/I3CLSimFlasherPulse.h>
#include <boost/preprocessor/seq.hpp>
     
#include <icetray/python/copy_suite.hpp>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/dataclass_suite.hpp>
     
namespace bp = boost::python;

#define ENUM_DEF(r,data,T) .value(BOOST_PP_STRINGIZE(T), data::T)


void register_I3CLSimFlasherPulse()
{
    {
        bp::scope I3CLSimFlasherPulse_scope = 
        bp::class_<I3CLSimFlasherPulse, boost::shared_ptr<I3CLSimFlasherPulse> >
        ("I3CLSimFlasherPulse")
        
        .def("GetType", &I3CLSimFlasherPulse::GetType)
        .def("SetType", &I3CLSimFlasherPulse::SetType)
        .add_property("type", &I3CLSimFlasherPulse::GetType, &I3CLSimFlasherPulse::SetType)

        .def("GetPos", (const I3Position& (I3CLSimFlasherPulse::*)()) &I3CLSimFlasherPulse::GetPos, bp::return_value_policy<bp::copy_const_reference>() )
        .def("SetPos", (void (I3CLSimFlasherPulse::*)(const I3Position&)) &I3CLSimFlasherPulse::SetPos)
        .def("GetDir", (const I3Direction& (I3CLSimFlasherPulse::*)()) &I3CLSimFlasherPulse::GetDir, bp::return_value_policy<bp::copy_const_reference>() )
        .def("SetDir", (void (I3CLSimFlasherPulse::*)(const I3Direction&)) &I3CLSimFlasherPulse::SetDir)

        .add_property("pos", bp::make_function( (const I3Position& (I3CLSimFlasherPulse::*)()) &I3CLSimFlasherPulse::GetPos, bp::return_value_policy<bp::copy_const_reference>() ),
                      (void (I3CLSimFlasherPulse::*)(const I3Position&)) &I3CLSimFlasherPulse::SetPos ) 
        .add_property("dir", bp::make_function( (const I3Direction& (I3CLSimFlasherPulse::*)()) &I3CLSimFlasherPulse::GetDir, bp::return_value_policy<bp::copy_const_reference>() ),
                      (void (I3CLSimFlasherPulse::*)(const I3Direction&)) &I3CLSimFlasherPulse::SetDir ) 

        .def("GetTime", &I3CLSimFlasherPulse::GetTime)
        .def("SetTime", &I3CLSimFlasherPulse::SetTime)
        .add_property("time", &I3CLSimFlasherPulse::GetTime, &I3CLSimFlasherPulse::SetTime)

        .def("GetNumberOfPhotonsNoBias", &I3CLSimFlasherPulse::GetNumberOfPhotonsNoBias)
        .def("SetNumberOfPhotonsNoBias", &I3CLSimFlasherPulse::SetNumberOfPhotonsNoBias)
        .add_property("numberOfPhotonsNoBias", &I3CLSimFlasherPulse::GetNumberOfPhotonsNoBias, &I3CLSimFlasherPulse::SetNumberOfPhotonsNoBias)

        .def("GetPulseWidth", &I3CLSimFlasherPulse::GetPulseWidth)
        .def("SetPulseWidth", &I3CLSimFlasherPulse::SetPulseWidth)
        .add_property("pulseWidth", &I3CLSimFlasherPulse::GetPulseWidth, &I3CLSimFlasherPulse::SetPulseWidth)

        .def("GetAngularEmissionSigmaPolar", &I3CLSimFlasherPulse::GetAngularEmissionSigmaPolar)
        .def("SetAngularEmissionSigmaPolar", &I3CLSimFlasherPulse::SetAngularEmissionSigmaPolar)
        .add_property("angularEmissionSigmaPolar", &I3CLSimFlasherPulse::GetAngularEmissionSigmaPolar, &I3CLSimFlasherPulse::SetAngularEmissionSigmaPolar)

        .def("GetAngularEmissionSigmaAzimuthal", &I3CLSimFlasherPulse::GetAngularEmissionSigmaAzimuthal)
        .def("SetAngularEmissionSigmaAzimuthal", &I3CLSimFlasherPulse::SetAngularEmissionSigmaAzimuthal)
        .add_property("angularEmissionSigmaAzimuthal", &I3CLSimFlasherPulse::GetAngularEmissionSigmaAzimuthal, &I3CLSimFlasherPulse::SetAngularEmissionSigmaAzimuthal)

        .def_pickle(bp::boost_serializable_pickle_suite<I3CLSimFlasherPulse>())
        .def(bp::copy_suite<I3CLSimFlasherPulse>())
        ;
        
        
        bp::enum_<I3CLSimFlasherPulse::FlasherPulseType>("FlasherPulseType")
        BOOST_PP_SEQ_FOR_EACH(ENUM_DEF,I3CLSimFlasherPulse,I3CLSIMFLASHERPULSE_H_I3CLSimFlasherPulse_FlasherPulseType)
        .export_values()
        ;

    }
    
    bp::class_<I3CLSimFlasherPulseSeries, bp::bases<I3FrameObject>, I3CLSimFlasherPulseSeriesPtr>("I3CLSimFlasherPulseSeries")
    .def(bp::dataclass_suite<I3CLSimFlasherPulseSeries>())
    
    ;

    bp::implicitly_convertible<shared_ptr<I3CLSimFlasherPulse>, shared_ptr<const I3CLSimFlasherPulse> >();
    register_pointer_conversions<I3CLSimFlasherPulseSeries>();
}
