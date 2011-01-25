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

#include <icetray/I3Units.h>

#include <clsim/I3CLSimPhoton.h>
#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>

using namespace boost::python;

static std::string 
i3clsimphoton_prettyprint(const I3CLSimPhoton& s)
{
    I3DirectionPtr dir = s.GetDir();
    if (!dir) dir=I3DirectionPtr(new I3Direction());
    
    std::ostringstream oss;
    oss << "[ I3CLSimPhoton id : " << s.GetID() << std::endl
        << "       pos (x,y,z) : [" << s.GetPosX()/I3Units::m << ", " << s.GetPosY()/I3Units::m << ", " << s.GetPosZ()/I3Units::m << "]m" << std::endl
        << "   dir (theta,phi) : [" << s.GetDirTheta()/I3Units::rad << ", " << s.GetDirPhi()/I3Units::rad << "]rad, [" << s.GetDirTheta()/I3Units::deg << ", " << s.GetDirPhi()/I3Units::deg << "]deg" << std::endl
        << "     dir (zen,azi) : [" << dir->GetZenith()/I3Units::rad << ", " << dir->GetAzimuth()/I3Units::rad << "]rad, [" << dir->GetZenith()/I3Units::deg << ", " << dir->GetAzimuth()/I3Units::deg << "]deg" << std::endl
        << "       dir (x,y,z) : [" << dir->GetX() << ", " << dir->GetY() << ", " << dir->GetZ() << "]" << std::endl
        << "              time : " << s.GetTime()/I3Units::ns << "ns" << std::endl
        << "        wavelength : " << s.GetWavelength()/I3Units::nanometer << "nm" << std::endl
        << "       numScatters : " << s.GetNumScatters() << std::endl
        << "            weight : " << s.GetWeight() << std::endl
        << "     cherenkovDist : " << s.GetCherenkovDist()/I3Units::m << "m" << std::endl
        << "]" ;
    
    return oss.str();
}


void register_I3CLSimPhoton()
{
    {
        void (I3CLSimPhoton::* SetDir_oneary)(const I3Direction&) = &I3CLSimPhoton::SetDir; 
        void (I3CLSimPhoton::* SetDir_threeary)(const double &x, const double &y, const double &z) = &I3CLSimPhoton::SetDir;
        
        scope clsimstep_scope = 
        class_<I3CLSimPhoton, boost::shared_ptr<I3CLSimPhoton> >("I3CLSimPhoton")
        .add_property("x", &I3CLSimPhoton::GetPosX, &I3CLSimPhoton::SetPosX)
        .add_property("y", &I3CLSimPhoton::GetPosY, &I3CLSimPhoton::SetPosY)
        .add_property("z", &I3CLSimPhoton::GetPosZ, &I3CLSimPhoton::SetPosZ)
        .add_property("time", &I3CLSimPhoton::GetTime, &I3CLSimPhoton::SetTime)

        .add_property("theta", &I3CLSimPhoton::GetDirTheta, &I3CLSimPhoton::SetDirTheta)
        .add_property("phi", &I3CLSimPhoton::GetDirPhi, &I3CLSimPhoton::SetDirPhi)
        .add_property("wavelength", &I3CLSimPhoton::GetWavelength, &I3CLSimPhoton::SetWavelength)
        .add_property("cherenkovDist", &I3CLSimPhoton::GetCherenkovDist, &I3CLSimPhoton::SetCherenkovDist)

        .add_property("numScatters", &I3CLSimPhoton::GetNumScatters, &I3CLSimPhoton::SetNumScatters)
        .add_property("weight", &I3CLSimPhoton::GetWeight, &I3CLSimPhoton::SetWeight)
        .add_property("id", &I3CLSimPhoton::GetID, &I3CLSimPhoton::SetID)

        .add_property("pos", &I3CLSimPhoton::GetPos, &I3CLSimPhoton::SetPos)
        .add_property("dir", &I3CLSimPhoton::GetDir, SetDir_oneary)
        
        .def("SetDirXYZ", SetDir_threeary)
        
        .def("__str__", i3clsimphoton_prettyprint)
        ;
    }
    
    class_<I3CLSimPhotonSeries, bases<I3FrameObject>, I3CLSimPhotonSeriesPtr>("I3CLSimPhotonSeries")
    .def(std_vector_indexing_suite<I3CLSimPhotonSeries>())
    ;
    
    // does not base on I3FrameObject, so register only the shared_ptr<T>-to-shared_ptr<const T> conversion
    //register_pointer_conversions<I3CLSimPhoton>();
    boost::python::implicitly_convertible<shared_ptr<I3CLSimPhoton>, shared_ptr<const I3CLSimPhoton> >();
    
    register_pointer_conversions<I3CLSimPhotonSeries>();
}
