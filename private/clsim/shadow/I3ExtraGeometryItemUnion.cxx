/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ExtraGeometryItemUnion.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 *
 *
 *  This file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *  
 */

#include <icetray/serialization.h>
#include <clsim/shadow/I3ExtraGeometryItemUnion.h>

#include <boost/foreach.hpp>

I3ExtraGeometryItemUnion::I3ExtraGeometryItemUnion() {;}

I3ExtraGeometryItemUnion::
I3ExtraGeometryItemUnion(const std::vector<I3ExtraGeometryItemConstPtr> &elements)
:
elements_(elements)
{
    
}

I3ExtraGeometryItemUnion::~I3ExtraGeometryItemUnion() { }


bool
I3ExtraGeometryItemUnion::DoesLineIntersect
(const I3Position &lineStart,
 const I3Position &lineEnd) const
{
    
    
    return false;
}

std::pair<I3Position, I3Position>
I3ExtraGeometryItemUnion::GetBoundingBox
() const
{
    
    
    return std::make_pair(I3Position(), I3Position());
}



template <class Archive>
void I3ExtraGeometryItemUnion::load(Archive &ar, unsigned version)
{
    if (version > i3extrageometryitemunion_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3ExtraGeometryItem class.",version,i3extrageometryitemunion_version_);
    
    ar >> make_nvp("I3ExtraGeometryItem", base_object<I3ExtraGeometryItem>(*this));
    
    {
        // work-around to allow de-serialization into vector<const ptr>
        std::vector<I3ExtraGeometryItemPtr> elements_non_const;
        ar >> make_nvp("elements", elements_non_const);
        BOOST_FOREACH(const I3ExtraGeometryItemPtr &item, elements_non_const) {
            elements_.push_back(item);
        }
    }
}


template <class Archive>
void I3ExtraGeometryItemUnion::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("I3ExtraGeometryItem", base_object<I3ExtraGeometryItem>(*this));
    ar << make_nvp("elements", elements_);
}


std::ostream& I3ExtraGeometryItemUnion::operator<<(std::ostream& oss) const
{
    oss << "[ I3ExtraGeometryItemUnion :" << std::endl;

    BOOST_FOREACH(const I3ExtraGeometryItemConstPtr &ptr, elements_)
    {
        oss << "  -> subitem: (not impl yet)" << std::endl;
    }
    
    oss << "]" << std::endl;
    return oss;
}


I3_SPLIT_SERIALIZABLE(I3ExtraGeometryItemUnion);
