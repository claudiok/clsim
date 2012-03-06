/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ExtraGeometryItemCylinder.cxx
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

#include <limits>

#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/shadow/I3ExtraGeometryItemCylinder.h>

#include <boost/foreach.hpp>

I3ExtraGeometryItemCylinder::I3ExtraGeometryItemCylinder()
:
radius_(NAN),
boundingBoxCalculated_(false)
{;}

I3ExtraGeometryItemCylinder::
I3ExtraGeometryItemCylinder(const I3Position &from, const I3Position &to, double radius)
:
from_(from),
to_(to),
radius_(radius),
boundingBoxCalculated_(false)
{
}

I3ExtraGeometryItemCylinder::~I3ExtraGeometryItemCylinder() { }

void I3ExtraGeometryItemCylinder::CalculateBoundingBox() const
{
    if (boundingBoxCalculated_) return;
    
    
    boundingBoxLower_.SetPos(NAN,NAN,NAN);
    boundingBoxUpper_.SetPos(NAN,NAN,NAN);
    
    boundingBoxCalculated_=true;
}


bool
I3ExtraGeometryItemCylinder::DoesLineIntersect
(const I3Position &lineStart,
 const I3Position &lineEnd) const
{
    if (isnan(radius_)) return false;
    
    
    
    return false;
}

std::pair<I3Position, I3Position>
I3ExtraGeometryItemCylinder::GetBoundingBox
() const
{
    CalculateBoundingBox();
    return std::make_pair(boundingBoxLower_, boundingBoxUpper_);
}



template <class Archive>
void I3ExtraGeometryItemCylinder::load(Archive &ar, unsigned version)
{
    if (version > i3extrageometryitemcylinder_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3ExtraGeometryItem class.",version,i3extrageometryitemcylinder_version_);
    
    ar >> make_nvp("I3ExtraGeometryItem", base_object<I3ExtraGeometryItem>(*this));

    ar >> make_nvp("from", from_);
    ar >> make_nvp("to", to_);
    ar >> make_nvp("radius", radius_);
}


template <class Archive>
void I3ExtraGeometryItemCylinder::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("I3ExtraGeometryItem", base_object<I3ExtraGeometryItem>(*this));

    ar << make_nvp("from", from_);
    ar << make_nvp("to", to_);
    ar << make_nvp("radius", radius_);
}


std::ostream& I3ExtraGeometryItemCylinder::operator<<(std::ostream& oss) const
{
    oss << "[ I3ExtraGeometryItemCylinder :" << std::endl;
    oss << "                         from : (" << from_.GetX()/I3Units::m << "," << from_.GetY()/I3Units::m << "," << from_.GetZ()/I3Units::m << ")m" << std::endl;
    oss << "                           to : (" << to_.GetX()/I3Units::m << "," << to_.GetY()/I3Units::m << "," << to_.GetZ()/I3Units::m << ")m" << std::endl;
    oss << "                       radius : " << radius_/I3Units::m << "m" << std::endl;
    oss << "]" << std::endl;
    return oss;
}


I3_SPLIT_SERIALIZABLE(I3ExtraGeometryItemCylinder);
