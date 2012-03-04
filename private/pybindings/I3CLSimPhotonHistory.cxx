//
//   Copyright (c) 2012  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module clsim.
//
//   clsim is free software; you can redistribute it and/or modify
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

#include <clsim/I3CLSimPhotonHistory.h>
#include <boost/preprocessor/seq.hpp>
#include <icetray/python/std_vector_indexing_suite.hpp>
#include <icetray/python/copy_suite.hpp>
#include <icetray/python/boost_serializable_pickle_suite.hpp>

using namespace boost::python;

static std::string 
i3clsimphotonhistory_prettyprint(const I3CLSimPhotonHistory& s)
{

    std::ostringstream oss;
    oss << "[ I3CLSimPhotonHistory  len=" << s.size() << std::endl;
    oss << "[ (most recent point first)" << std::endl;
    
    for (std::size_t i=0;i<s.size();++i)
    {
        oss << "  index " << i << " : pos=(" << s.GetX(i)/I3Units::m << "," << s.GetY(i)/I3Units::m << "," << s.GetZ(i)/I3Units::m << ")m" << std::endl;
    }

    oss << "]" << std::endl;
    return oss.str();
}


void register_I3CLSimPhotonHistory()
{
    {
        void (I3CLSimPhotonHistory::* append_oneary)(const I3Position &) = &I3CLSimPhotonHistory::push_back;
        void (I3CLSimPhotonHistory::* append_threeary)(float x, float y, float z) = &I3CLSimPhotonHistory::push_back;

        scope clsimphotonhistory_scope = 
        class_<I3CLSimPhotonHistory, boost::shared_ptr<I3CLSimPhotonHistory> >("I3CLSimPhotonHistory")
        .def("append", append_oneary)
        .def("append", append_threeary)
        .def("__len__", &I3CLSimPhotonHistory::size)
        .def("__getitem__", &I3CLSimPhotonHistory::at)
        
        .def("__str__", i3clsimphotonhistory_prettyprint)

        .def(bp::copy_suite<I3CLSimPhotonHistory>())
        .def_pickle(bp::boost_serializable_pickle_suite<I3CLSimPhotonHistory>())
        ;
    }
    
    class_<I3CLSimPhotonHistorySeries, bases<I3FrameObject>, I3CLSimPhotonHistorySeriesPtr>("I3CLSimPhotonHistorySeries")
    .def(std_vector_indexing_suite<I3CLSimPhotonHistorySeries>())
    .def_pickle(bp::boost_serializable_pickle_suite<I3CLSimPhotonHistorySeries>())
    ;

    // does not base on I3FrameObject, so register only the shared_ptr<T>-to-shared_ptr<const T> conversion
    //register_pointer_conversions<I3CLSimPhoton>();
    boost::python::implicitly_convertible<shared_ptr<I3CLSimPhotonHistory>, shared_ptr<const I3CLSimPhotonHistory> >();
    
    register_pointer_conversions<I3CLSimPhotonHistorySeries>();
}
