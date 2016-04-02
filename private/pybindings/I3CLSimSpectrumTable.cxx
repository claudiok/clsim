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
 * @file I3CLSimSpectrumTable.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/I3CLSimSpectrumTable.h>
#include <boost/foreach.hpp>

namespace bp = boost::python;

namespace {
    bp::list I3CLSimSpectrumTable_GetSpectra_python(const I3CLSimSpectrumTable &self)
    {
        const std::vector<I3CLSimFunctionConstPtr> &spectra = self.GetSpectra();

        bp::list t;
        
        BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, spectra)
        {
            t.append(ptr);
        }
        
        return t;
    }
}

void register_I3CLSimSpectrumTable()
{
    {
        bp::scope i3clsimspectrumtable_scope = 
        bp::class_<I3CLSimSpectrumTable, boost::shared_ptr<I3CLSimSpectrumTable> >
        ("I3CLSimSpectrumTable")

        .add_property("spectra", I3CLSimSpectrumTable_GetSpectra_python)
        .def("GetSpectra", I3CLSimSpectrumTable_GetSpectra_python)

        .def("append", &I3CLSimSpectrumTable::append)

        .def("__len__", &I3CLSimSpectrumTable::size)
        .def("__getitem__", &I3CLSimSpectrumTable::operator[])
        ;
        
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimSpectrumTable>, boost::shared_ptr<const I3CLSimSpectrumTable> >();
}
