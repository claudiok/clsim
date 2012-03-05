/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ExtraGeometryItem.cxx
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
#include <clsim/shadow/I3ExtraGeometryItem.h>

I3ExtraGeometryItem::~I3ExtraGeometryItem() { }


template <class Archive>
void I3ExtraGeometryItem::serialize (Archive &ar, unsigned version)
{
    if (version > i3extrageometryitem_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3ExtraGeometryItem class.",version,i3extrageometryitem_version_);
    
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    
        
}


std::ostream& I3ExtraGeometryItem::operator<<(std::ostream& oss) const
{
    oss << "[ I3ExtraGeometryItem ]" << std::endl;
    return oss;
}

I3_SERIALIZABLE(I3ExtraGeometryItem);
