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

#include <clsim/I3Photon.h>
#include <boost/preprocessor/seq.hpp>

#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

static std::string 
i3photon_prettyprint(const I3Photon& s)
{
    std::ostringstream oss;
    oss << "[      I3Photon id : "  << s.GetID() << std::endl
        << "       pos (x,y,z) : [" << s.GetPos().GetX()/I3Units::m << ", " << s.GetPos().GetY()/I3Units::m << ", " << s.GetPos().GetZ()/I3Units::m << "]m" << std::endl
        << "   dir (theta,phi) : [" << s.GetDir().CalcTheta()/I3Units::rad << ", " << s.GetDir().CalcPhi()/I3Units::rad << "]rad, [" << s.GetDir().CalcTheta()/I3Units::deg << ", " << s.GetDir().CalcPhi()/I3Units::deg << "]deg" << std::endl
        << "     dir (zen,azi) : [" << s.GetDir().GetZenith()/I3Units::rad << ", " << s.GetDir().GetAzimuth()/I3Units::rad << "]rad, [" << s.GetDir().GetZenith()/I3Units::deg << ", " << s.GetDir().GetAzimuth()/I3Units::deg << "]deg" << std::endl
        << "       dir (x,y,z) : [" << s.GetDir().GetX() << ", " << s.GetDir().GetY() << ", " << s.GetDir().GetZ() << "]" << std::endl
        << "              time : " << s.GetTime()/I3Units::ns << "ns" << std::endl

        << "  startPos (x,y,z) : [" << s.GetStartPos().GetX()/I3Units::m << ", " << s.GetStartPos().GetY()/I3Units::m << ", " << s.GetStartPos().GetZ()/I3Units::m << "]m" << std::endl
        << "strDir (theta,phi) : [" << s.GetStartDir().CalcTheta()/I3Units::rad << ", " << s.GetStartDir().CalcPhi()/I3Units::rad << "]rad, [" << s.GetStartDir().CalcTheta()/I3Units::deg << ", " << s.GetStartDir().CalcPhi()/I3Units::deg << "]deg" << std::endl
        << "startDir (zen,azi) : [" << s.GetStartDir().GetZenith()/I3Units::rad << ", " << s.GetStartDir().GetAzimuth()/I3Units::rad << "]rad, [" << s.GetStartDir().GetZenith()/I3Units::deg << ", " << s.GetStartDir().GetAzimuth()/I3Units::deg << "]deg" << std::endl
        << "  startDir (x,y,z) : [" << s.GetStartDir().GetX() << ", " << s.GetStartDir().GetY() << ", " << s.GetStartDir().GetZ() << "]" << std::endl
        << "         startTime : " << s.GetStartTime()/I3Units::ns << "ns" << std::endl
    
        << "        wavelength : " << s.GetWavelength()/I3Units::nanometer << "nm" << std::endl
        << "    group velocity : " << s.GetGroupVelocity()/(I3Units::meter/I3Units::nanosecond) << "m/ns" << std::endl
        << "      numScattered : " << s.GetNumScattered() << std::endl
        << "     distInAbsLens : " << s.GetDistanceInAbsorptionLengths() << std::endl
        << "            weight : " << s.GetWeight() << std::endl
        << "     cherenkovDist : " << s.GetCherenkovDist()/I3Units::m << "m" << std::endl
        << "     cherenkovTime : " << s.GetCherenkovTime()/I3Units::ns << "ns" << std::endl
        << "   particleMajorID : " << s.GetParticleMajorID() << std::endl
        << "   particleMinorID : " << s.GetParticleMinorID() << std::endl
        << " len(positionList) : " << s.GetNumPositionListEntries() << std::endl;

    for (uint32_t i=0;i<s.GetNumPositionListEntries();++i)
    {
        I3PositionConstPtr pos = s.GetPositionListEntry(i);
        double distInAbsLensAtPos = s.GetDistanceInAbsorptionLengthsAtPositionListEntry(i);
        if (i==0) {
            oss << "         (initial) : ";
        } else if (i==s.GetNumPositionListEntries()-1) {
            oss << "           (final) : ";
        } else {
            oss << "                   : ";
        }

        if (!pos) {
            oss << "(not saved)" << std::endl;
        } else {
            oss << "(" << pos->GetX()/I3Units::m << "," << pos->GetY()/I3Units::m << "," << pos->GetZ()/I3Units::m << ")m; absLens=" << distInAbsLensAtPos << std::endl; 
        }
    }
    
    oss << "]" ;
    
    return oss.str();
}


namespace I3Photon_python_helper
{
    boost::python::list GetPositionList(const I3Photon &photon)
    {
        boost::python::list pylist;
        
        for (uint32_t i=0;i<photon.GetNumPositionListEntries();++i)
        {
            I3PositionConstPtr posListPhoton = photon.GetPositionListEntry(i);
            if (posListPhoton) {
                pylist.append( I3PositionPtr(new I3Position(*posListPhoton)) );
            } else {
                pylist.append( I3PositionPtr() );
            }
        }
        
        return pylist;
    }
}

void register_I3Photon()
{
    {
        scope i3photon_scope = 
        class_<I3Photon, boost::shared_ptr<I3Photon> >("I3Photon")
        .add_property("pos", 
                      make_function(&I3Photon::GetPos, return_value_policy<copy_const_reference>()),
                      make_function(&I3Photon::SetPos))
        .def("GetPos", &I3Photon::GetPos, return_value_policy<copy_const_reference>())
        .def("SetPos", &I3Photon::SetPos)

        .add_property("dir", 
                      make_function(&I3Photon::GetDir, return_value_policy<copy_const_reference>()),
                      make_function(&I3Photon::SetDir))
        .def("GetDir", &I3Photon::GetDir, return_value_policy<copy_const_reference>())
        .def("SetDir", &I3Photon::SetDir)

        .add_property("startPos", 
                      make_function(&I3Photon::GetStartPos, return_value_policy<copy_const_reference>()),
                      make_function(&I3Photon::SetStartPos))
        .def("GetStartPos", &I3Photon::GetStartPos, return_value_policy<copy_const_reference>())
        .def("SetStartPos", &I3Photon::SetStartPos)
        
        .add_property("startDir", 
                      make_function(&I3Photon::GetStartDir, return_value_policy<copy_const_reference>()),
                      make_function(&I3Photon::SetStartDir))
        .def("GetStartDir", &I3Photon::GetStartDir, return_value_policy<copy_const_reference>())
        .def("SetStartDir", &I3Photon::SetStartDir)

        .add_property("numScattered", &I3Photon::GetNumScattered, &I3Photon::SetNumScattered)
        .def("GetNumScattered", &I3Photon::GetNumScattered)
        .def("SetNumScattered", &I3Photon::SetNumScattered)

        .add_property("distanceInAbsorptionLengths", &I3Photon::GetDistanceInAbsorptionLengths, &I3Photon::SetDistanceInAbsorptionLengths)
        .def("GetDistanceInAbsorptionLengths", &I3Photon::GetDistanceInAbsorptionLengths)
        .def("SetDistanceInAbsorptionLengths", &I3Photon::SetDistanceInAbsorptionLengths)

        .add_property("wavelength", &I3Photon::GetWavelength, &I3Photon::SetWavelength)
        .def("GetWavelength", &I3Photon::GetWavelength)
        .def("SetWavelength", &I3Photon::SetWavelength)

        .add_property("groupVelocity", &I3Photon::GetGroupVelocity, &I3Photon::SetGroupVelocity)
        .def("GetGroupVelocity", &I3Photon::GetGroupVelocity)
        .def("SetGroupVelocity", &I3Photon::SetGroupVelocity)

        .add_property("time", &I3Photon::GetTime, &I3Photon::SetTime)
        .def("GetTime", &I3Photon::GetTime)
        .def("SetTime", &I3Photon::SetTime)

        .add_property("startTime", &I3Photon::GetStartTime, &I3Photon::SetStartTime)
        .def("GetStartTime", &I3Photon::GetStartTime)
        .def("SetStartTime", &I3Photon::SetStartTime)

        .add_property("cherenkovDist", &I3Photon::GetCherenkovDist, &I3Photon::SetCherenkovDist)
        .def("GetCherenkovDist", &I3Photon::GetCherenkovDist)
        .def("SetCherenkovDist", &I3Photon::SetCherenkovDist)

        .add_property("cherenkovTime", &I3Photon::GetCherenkovTime)
        .def("GetCherenkovTime", &I3Photon::GetCherenkovTime)
        
        .add_property("ID", &I3Photon::GetID, &I3Photon::SetID)
        .def("GetID", &I3Photon::GetID)
        .def("SetID", &I3Photon::SetID)

        .add_property("weight", &I3Photon::GetWeight, &I3Photon::SetWeight)
        .def("GetWeight", &I3Photon::GetWeight)
        .def("SetWeight", &I3Photon::SetWeight)

        .add_property("particleMajorID", &I3Photon::GetParticleMajorID, &I3Photon::SetParticleMajorID)
        .def("GetParticleMajorID", &I3Photon::GetParticleMajorID)
        .def("SetParticleMajorID", &I3Photon::SetParticleMajorID)

        .add_property("particleMinorID", &I3Photon::GetParticleMinorID, &I3Photon::SetParticleMinorID)
        .def("GetParticleMinorID", &I3Photon::GetParticleMinorID)
        .def("SetParticleMinorID", &I3Photon::SetParticleMinorID)

        .def("SetParticleID", &I3Photon::SetParticleID)

        .add_property("numPositionListEntries", &I3Photon::GetNumPositionListEntries)
        .def("GetNumPositionListEntries", &I3Photon::GetNumPositionListEntries)
        .def("GetPositionList", &I3Photon_python_helper::GetPositionList)
        .add_property("positionList", &I3Photon_python_helper::GetPositionList)
        .def("GetPositionListEntry", &I3Photon::GetPositionListEntry)
        .def("GetDistanceInAbsorptionLengthsAtPositionListEntry", &I3Photon::GetDistanceInAbsorptionLengthsAtPositionListEntry)
        .def("AppendToIntermediatePositionList", &I3Photon::AppendToIntermediatePositionList)

        .def("__str__", i3photon_prettyprint)

        .def(dataclass_suite<I3Photon>())
        ;
    }
    
    class_<I3PhotonSeries, bases<I3FrameObject>, I3PhotonSeriesPtr>("I3PhotonSeries")
    .def(dataclass_suite<I3PhotonSeries>())
    ;

    class_<I3PhotonSeriesMap, bases<I3FrameObject>, I3PhotonSeriesMapPtr>("I3PhotonSeriesMap")
    .def(dataclass_suite<I3PhotonSeriesMap>())
    ;

    register_pointer_conversions<I3Photon>();
    
    register_pointer_conversions<I3PhotonSeries>();
    register_pointer_conversions<I3PhotonSeriesMap>();
}
