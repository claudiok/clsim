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

#include <clsim/I3Photon.h>
#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>
#include <icetray/python/std_map_indexing_suite.hpp>

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
        << "    group velocity : " << s.GetGroupVelocity()/I3Units::nanometer << "nm" << std::endl
        << "      numScattered : " << s.GetNumScattered() << std::endl
        << "            weight : " << s.GetWeight() << std::endl
        << "     cherenkovDist : " << s.GetCherenkovDist()/I3Units::m << "m" << std::endl
        << "     cherenkovTime : " << s.GetCherenkovTime()/I3Units::ns << "ns" << std::endl
        << "   particleMajorID : " << s.GetParticleMajorID() << std::endl
        << "   particleMinorID : " << s.GetParticleMinorID() << std::endl
        << "]" ;
    
    return oss.str();
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

        .add_property("particleMajorID", &I3Photon::GetParticleMajorID, &I3Photon::SetParticleMajorID)
        .def("GetParticleMajorID", &I3Photon::GetParticleMajorID)
        .def("SetParticleMajorID", &I3Photon::SetParticleMajorID)

        .add_property("particleMinorID", &I3Photon::GetParticleMinorID, &I3Photon::SetParticleMinorID)
        .def("GetParticleMinorID", &I3Photon::GetParticleMinorID)
        .def("SetParticleMinorID", &I3Photon::SetParticleMinorID)

        .def("SetParticleID", &I3Photon::SetParticleID)

        .def("__str__", i3photon_prettyprint)
        ;
    }
    
    class_<I3PhotonSeries, bases<I3FrameObject>, I3PhotonSeriesPtr>("I3PhotonSeries")
    .def(std_vector_indexing_suite<I3PhotonSeries>())
    ;

    class_<I3PhotonSeriesMap, bases<I3FrameObject>, I3PhotonSeriesMapPtr>("I3PhotonSeriesMap")
    .def(std_map_indexing_suite<I3PhotonSeriesMap>())
    ;

    register_pointer_conversions<I3Photon>();
    
    register_pointer_conversions<I3PhotonSeries>();
    register_pointer_conversions<I3PhotonSeriesMap>();
}
