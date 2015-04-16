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
 * @file I3CLSimPhoton.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <icetray/I3Units.h>

#include <clsim/I3CLSimPhoton.h>
#include <boost/preprocessor/seq.hpp>

#include <icetray/python/list_indexing_suite.hpp>
#include <icetray/python/std_map_indexing_suite.hpp>
#include <icetray/python/copy_suite.hpp>
#include <icetray/python/boost_serializable_pickle_suite.hpp>

using namespace boost::python;

static std::string 
i3clsimphoton_prettyprint(const I3CLSimPhoton& s)
{
    I3DirectionPtr dir = s.GetDir();
    if (!dir) dir=I3DirectionPtr(new I3Direction());

    I3DirectionPtr startDir = s.GetStartDir();
    if (!startDir) startDir=I3DirectionPtr(new I3Direction());

    std::ostringstream oss;
    oss << "[ I3CLSimPhoton id : " << s.GetID() << std::endl
        << "       pos (x,y,z) : [" << s.GetPosX()/I3Units::m << ", " << s.GetPosY()/I3Units::m << ", " << s.GetPosZ()/I3Units::m << "]m" << std::endl
        << "   dir (theta,phi) : [" << s.GetDirTheta()/I3Units::rad << ", " << s.GetDirPhi()/I3Units::rad << "]rad, [" << s.GetDirTheta()/I3Units::deg << ", " << s.GetDirPhi()/I3Units::deg << "]deg" << std::endl
        << "     dir (zen,azi) : [" << dir->GetZenith()/I3Units::rad << ", " << dir->GetAzimuth()/I3Units::rad << "]rad, [" << dir->GetZenith()/I3Units::deg << ", " << dir->GetAzimuth()/I3Units::deg << "]deg" << std::endl
        << "       dir (x,y,z) : [" << dir->GetX() << ", " << dir->GetY() << ", " << dir->GetZ() << "]" << std::endl
        << "              time : " << s.GetTime()/I3Units::ns << "ns" << std::endl

        << "  startPos (x,y,z) : [" << s.GetStartPosX()/I3Units::m << ", " << s.GetStartPosY()/I3Units::m << ", " << s.GetStartPosZ()/I3Units::m << "]m" << std::endl
        << "strDir (theta,phi) : [" << s.GetStartDirTheta()/I3Units::rad << ", " << s.GetStartDirPhi()/I3Units::rad << "]rad, [" << s.GetStartDirTheta()/I3Units::deg << ", " << s.GetStartDirPhi()/I3Units::deg << "]deg" << std::endl
        << "startDir (zen,azi) : [" << startDir->GetZenith()/I3Units::rad << ", " << startDir->GetAzimuth()/I3Units::rad << "]rad, [" << startDir->GetZenith()/I3Units::deg << ", " << startDir->GetAzimuth()/I3Units::deg << "]deg" << std::endl
        << "  startDir (x,y,z) : [" << startDir->GetX() << ", " << startDir->GetY() << ", " << startDir->GetZ() << "]" << std::endl
        << "         startTime : " << s.GetStartTime()/I3Units::ns << "ns" << std::endl

        << "        wavelength : " << s.GetWavelength()/I3Units::nanometer << "nm" << std::endl
        << "       numScatters : " << s.GetNumScatters() << std::endl
        << "            weight : " << s.GetWeight() << std::endl
        << "     cherenkovDist : " << s.GetCherenkovDist()/I3Units::m << "m" << std::endl
        << "     groupVelocity : " << s.GetGroupVelocity()/(I3Units::m/I3Units::ns) << "m/ns" << std::endl
        << "     distInAbsLens : " << s.GetDistInAbsLens() << " absorption lengths" << std::endl
    
        << "          stringID : " << s.GetStringID() << std::endl
        << "              omID : " << s.GetOMID() << std::endl
        << "]" ;
    
    return oss.str();
}


void register_I3CLSimPhoton()
{
    {
        void (I3CLSimPhoton::* SetDir_oneary)(const I3Direction&) = &I3CLSimPhoton::SetDir; 
        void (I3CLSimPhoton::* SetDir_threeary)(const double &x, const double &y, const double &z) = &I3CLSimPhoton::SetDir;

        void (I3CLSimPhoton::* SetStartDir_oneary)(const I3Direction&) = &I3CLSimPhoton::SetStartDir; 
        void (I3CLSimPhoton::* SetStartDir_threeary)(const double &x, const double &y, const double &z) = &I3CLSimPhoton::SetStartDir;

        scope clsimstep_scope = 
        class_<I3CLSimPhoton, boost::shared_ptr<I3CLSimPhoton> >("I3CLSimPhoton")
        .add_property("x", &I3CLSimPhoton::GetPosX, &I3CLSimPhoton::SetPosX)
        .add_property("y", &I3CLSimPhoton::GetPosY, &I3CLSimPhoton::SetPosY)
        .add_property("z", &I3CLSimPhoton::GetPosZ, &I3CLSimPhoton::SetPosZ)
        .add_property("time", &I3CLSimPhoton::GetTime, &I3CLSimPhoton::SetTime)

        .add_property("startX", &I3CLSimPhoton::GetStartPosX, &I3CLSimPhoton::SetStartPosX)
        .add_property("startY", &I3CLSimPhoton::GetStartPosY, &I3CLSimPhoton::SetStartPosY)
        .add_property("startZ", &I3CLSimPhoton::GetStartPosZ, &I3CLSimPhoton::SetStartPosZ)
        .add_property("startTime", &I3CLSimPhoton::GetStartTime, &I3CLSimPhoton::SetStartTime)

        .add_property("theta", &I3CLSimPhoton::GetDirTheta, &I3CLSimPhoton::SetDirTheta)
        .add_property("phi", &I3CLSimPhoton::GetDirPhi, &I3CLSimPhoton::SetDirPhi)

        .add_property("startTheta", &I3CLSimPhoton::GetStartDirTheta, &I3CLSimPhoton::SetStartDirTheta)
        .add_property("startPhi", &I3CLSimPhoton::GetStartDirPhi, &I3CLSimPhoton::SetStartDirPhi)
        
        .add_property("wavelength", &I3CLSimPhoton::GetWavelength, &I3CLSimPhoton::SetWavelength)
        .add_property("cherenkovDist", &I3CLSimPhoton::GetCherenkovDist, &I3CLSimPhoton::SetCherenkovDist)
        .add_property("groupVelocity", &I3CLSimPhoton::GetGroupVelocity, &I3CLSimPhoton::SetGroupVelocity)
        .add_property("distInAbsLens", &I3CLSimPhoton::GetDistInAbsLens, &I3CLSimPhoton::SetDistInAbsLens)

        .add_property("numScatters", &I3CLSimPhoton::GetNumScatters, &I3CLSimPhoton::SetNumScatters)
        .add_property("weight", &I3CLSimPhoton::GetWeight, &I3CLSimPhoton::SetWeight)
        .add_property("id", &I3CLSimPhoton::GetID, &I3CLSimPhoton::SetID)
        .add_property("stringID", &I3CLSimPhoton::GetStringID, &I3CLSimPhoton::SetStringID)
        .add_property("omID", &I3CLSimPhoton::GetOMID, &I3CLSimPhoton::SetOMID)

        .add_property("pos", &I3CLSimPhoton::GetPos, &I3CLSimPhoton::SetPos)
        .add_property("dir", &I3CLSimPhoton::GetDir, SetDir_oneary)

        .add_property("startPos", &I3CLSimPhoton::GetStartPos, &I3CLSimPhoton::SetStartPos)
        .add_property("startDir", &I3CLSimPhoton::GetStartDir, SetStartDir_oneary)

        .def("SetDirXYZ", SetDir_threeary)
        .def("SetStartDirXYZ", SetStartDir_threeary)
        
        .def("__str__", i3clsimphoton_prettyprint)

        .def(bp::copy_suite<I3CLSimPhoton>())
        .def_pickle(bp::boost_serializable_pickle_suite<I3CLSimPhoton>())
        ;
    }
    
    class_<I3CLSimPhotonSeries, bases<I3FrameObject>, I3CLSimPhotonSeriesPtr>("I3CLSimPhotonSeries")
    .def(list_indexing_suite<I3CLSimPhotonSeries>())
    .def_pickle(bp::boost_serializable_pickle_suite<I3CLSimPhotonSeries>())
    ;

    class_<I3CLSimPhotonSeriesMap, bases<I3FrameObject>, I3CLSimPhotonSeriesMapPtr>("I3CLSimPhotonSeriesMap")
    .def(std_map_indexing_suite<I3CLSimPhotonSeriesMap>())
    .def_pickle(bp::boost_serializable_pickle_suite<I3CLSimPhotonSeriesMap>())
    ;

    // does not base on I3FrameObject, so register only the shared_ptr<T>-to-shared_ptr<const T> conversion
    //register_pointer_conversions<I3CLSimPhoton>();
    boost::python::implicitly_convertible<shared_ptr<I3CLSimPhoton>, shared_ptr<const I3CLSimPhoton> >();
    
    register_pointer_conversions<I3CLSimPhotonSeries>();
    register_pointer_conversions<I3CLSimPhotonSeriesMap>();
}
