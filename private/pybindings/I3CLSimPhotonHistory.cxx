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
 * @file I3CLSimPhotonHistory.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <icetray/I3Units.h>

#include <clsim/I3CLSimPhotonHistory.h>
#include <boost/preprocessor/seq.hpp>

#include <icetray/python/list_indexing_suite.hpp>
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
        oss << "  index " << i << " : pos=(" << s.GetX(i)/I3Units::m << "," << s.GetY(i)/I3Units::m << "," << s.GetZ(i)/I3Units::m << ")m; abslens=" << s.GetDistanceInAbsorptionLengths(i) << std::endl;
    }

    oss << "]" << std::endl;
    return oss.str();
}


void register_I3CLSimPhotonHistory()
{
    {
        void (I3CLSimPhotonHistory::* append_twoary)(const I3Position &, float) = &I3CLSimPhotonHistory::push_back;
        void (I3CLSimPhotonHistory::* append_fourary)(float, float, float, float) = &I3CLSimPhotonHistory::push_back;

        scope clsimphotonhistory_scope = 
        class_<I3CLSimPhotonHistory, boost::shared_ptr<I3CLSimPhotonHistory> >("I3CLSimPhotonHistory")
        .def("append", append_twoary)
        .def("append", append_fourary)
        .def("__len__", &I3CLSimPhotonHistory::size)
        .def("__getitem__", &I3CLSimPhotonHistory::at)
        
        .def("__str__", i3clsimphotonhistory_prettyprint)

        .def(bp::copy_suite<I3CLSimPhotonHistory>())
        .def_pickle(bp::boost_serializable_pickle_suite<I3CLSimPhotonHistory>())
        ;
    }
    
    class_<I3CLSimPhotonHistorySeries, bases<I3FrameObject>, I3CLSimPhotonHistorySeriesPtr>("I3CLSimPhotonHistorySeries")
    .def(list_indexing_suite<I3CLSimPhotonHistorySeries>())
    .def_pickle(bp::boost_serializable_pickle_suite<I3CLSimPhotonHistorySeries>())
    ;

    // does not base on I3FrameObject, so register only the shared_ptr<T>-to-shared_ptr<const T> conversion
    //register_pointer_conversions<I3CLSimPhoton>();
    boost::python::implicitly_convertible<shared_ptr<I3CLSimPhotonHistory>, shared_ptr<const I3CLSimPhotonHistory> >();
    
    register_pointer_conversions<I3CLSimPhotonHistorySeries>();
}
