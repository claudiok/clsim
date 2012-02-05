//
//   Copyright (c) 2012  Claudio Kopper
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

#include <clsim/I3CLSimSpectrumTable.h>
#include <boost/foreach.hpp>

namespace bp = boost::python;

namespace {
    bp::list I3CLSimSpectrumTable_GetSpectra_python(const I3CLSimSpectrumTable &self)
    {
        const std::vector<I3CLSimWlenDependentValueConstPtr> &spectra = self.GetSpectra();

        bp::list t;
        
        BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, spectra)
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
    
    bp::implicitly_convertible<shared_ptr<I3CLSimSpectrumTable>, shared_ptr<const I3CLSimSpectrumTable> >();
}
