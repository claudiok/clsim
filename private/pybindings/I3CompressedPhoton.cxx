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
 * @file I3Photon.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <icetray/I3Units.h>

#include <clsim/I3CompressedPhoton.h>
#include <boost/preprocessor/seq.hpp>

#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

static std::string 
i3compressedphoton_prettyprint(const I3CompressedPhoton& s)
{
    std::ostringstream oss;
    oss << "[ I3CompressedPhoton    : "  << std::endl
        << "            pos (x,y,z) : [" << s.GetPos().GetX()/I3Units::m << ", " << s.GetPos().GetY()/I3Units::m << ", " << s.GetPos().GetZ()/I3Units::m << "]m" << std::endl
        << "        dir (theta,phi) : [" << s.GetDir().CalcTheta()/I3Units::rad << ", " << s.GetDir().CalcPhi()/I3Units::rad << "]rad, [" << s.GetDir().CalcTheta()/I3Units::deg << ", " << s.GetDir().CalcPhi()/I3Units::deg << "]deg" << std::endl
        << "          dir (zen,azi) : [" << s.GetDir().GetZenith()/I3Units::rad << ", " << s.GetDir().GetAzimuth()/I3Units::rad << "]rad, [" << s.GetDir().GetZenith()/I3Units::deg << ", " << s.GetDir().GetAzimuth()/I3Units::deg << "]deg" << std::endl
        << "            dir (x,y,z) : [" << s.GetDir().GetX() << ", " << s.GetDir().GetY() << ", " << s.GetDir().GetZ() << "]" << std::endl
        << "                   time : " << s.GetTime()/I3Units::ns << "ns" << std::endl

        << "             wavelength : " << s.GetWavelength()/I3Units::nanometer << "nm" << std::endl
        << "                 weight : " << s.GetWeight() << std::endl
        << "        particleMajorID : " << s.GetParticleMajorID() << std::endl
        << "        particleMinorID : " << s.GetParticleMinorID() << std::endl
        << "]" ;
    
    return oss.str();
}


void register_I3CompressedPhoton()
{
    {
        scope i3compressedphoton_scope = 
        class_<I3CompressedPhoton, boost::shared_ptr<I3CompressedPhoton> >("I3CompressedPhoton", init<>() )
        .def(init<const I3Photon &>( (bp::arg("p")) ) )

        .add_property("pos", 
                      make_function(&I3CompressedPhoton::GetPos),
                      make_function(&I3CompressedPhoton::SetPos))
        .def("GetPos", &I3CompressedPhoton::GetPos)
        .def("SetPos", &I3CompressedPhoton::SetPos)

        .add_property("dir", 
                      make_function(&I3CompressedPhoton::GetDir),
                      make_function(&I3CompressedPhoton::SetDir))
        .def("GetDir", &I3CompressedPhoton::GetDir)
        .def("SetDir", &I3CompressedPhoton::SetDir)

        .add_property("wavelength", &I3CompressedPhoton::GetWavelength, &I3CompressedPhoton::SetWavelength)
        .def("GetWavelength", &I3CompressedPhoton::GetWavelength)
        .def("SetWavelength", &I3CompressedPhoton::SetWavelength)

        .add_property("time", &I3CompressedPhoton::GetTime, &I3CompressedPhoton::SetTime)
        .def("GetTime", &I3CompressedPhoton::GetTime)
        .def("SetTime", &I3CompressedPhoton::SetTime)

        .add_property("weight", &I3CompressedPhoton::GetWeight, &I3CompressedPhoton::SetWeight)
        .def("GetWeight", &I3CompressedPhoton::GetWeight)
        .def("SetWeight", &I3CompressedPhoton::SetWeight)

        .add_property("particleMajorID", &I3CompressedPhoton::GetParticleMajorID, &I3CompressedPhoton::SetParticleMajorID)
        .def("GetParticleMajorID", &I3CompressedPhoton::GetParticleMajorID)
        .def("SetParticleMajorID", &I3CompressedPhoton::SetParticleMajorID)

        .add_property("particleMinorID", &I3CompressedPhoton::GetParticleMinorID, &I3CompressedPhoton::SetParticleMinorID)
        .def("GetParticleMinorID", &I3CompressedPhoton::GetParticleMinorID)
        .def("SetParticleMinorID", &I3CompressedPhoton::SetParticleMinorID)

        .def("SetParticleID", &I3CompressedPhoton::SetParticleID)

        .def("__str__", i3compressedphoton_prettyprint)

        .def(dataclass_suite<I3CompressedPhoton>())
        ;
    }
    
    class_<I3CompressedPhotonSeries, bases<I3FrameObject>, I3CompressedPhotonSeriesPtr>("I3CompressedPhotonSeries")
    .def(dataclass_suite<I3CompressedPhotonSeries>())
    ;

    class_<I3CompressedPhotonSeriesMap, bases<I3FrameObject>, I3CompressedPhotonSeriesMapPtr>("I3CompressedPhotonSeriesMap")
    .def(dataclass_suite<I3CompressedPhotonSeriesMap>())
    ;

    // not an I3FrameObject:
    // register_pointer_conversions<I3CompressedPhoton>();
    
    register_pointer_conversions<I3CompressedPhotonSeries>();
    register_pointer_conversions<I3CompressedPhotonSeriesMap>();
}
