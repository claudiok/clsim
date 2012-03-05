/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ExtraGeometryItemMove.cxx
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
#include <clsim/shadow/I3ExtraGeometryItemMove.h>

#include <boost/foreach.hpp>

I3ExtraGeometryItemMove::I3ExtraGeometryItemMove()
:
boundingBoxCalculated_(false)
{;}

I3ExtraGeometryItemMove::
I3ExtraGeometryItemMove(I3ExtraGeometryItemConstPtr element, const I3Position &offset)
:
element_(element),
offset_(offset),
boundingBoxCalculated_(false)
{
}

I3ExtraGeometryItemMove::~I3ExtraGeometryItemMove() { }

void I3ExtraGeometryItemMove::CalculateBoundingBox() const
{
    if (boundingBoxCalculated_) return;
    
    if (!element_) {
        boundingBoxLower_.SetPos(NAN, NAN, NAN);
        boundingBoxUpper_.SetPos(NAN, NAN, NAN);
        
        boundingBoxCalculated_=true;
        return;
    }
    
    const std::pair<I3Position, I3Position> box =
    element_->GetBoundingBox();
    
    boundingBoxLower_.SetPos(box.first.GetX() +offset_.GetX(), box.first.GetY() +offset_.GetY(), box.first.GetZ() +offset_.GetZ());
    boundingBoxUpper_.SetPos(box.second.GetX()+offset_.GetX(), box.second.GetY()+offset_.GetY(), box.second.GetZ()+offset_.GetZ());
    
    boundingBoxCalculated_=true;
}


bool
I3ExtraGeometryItemMove::DoesLineIntersect
(const I3Position &lineStart,
 const I3Position &lineEnd) const
{
    if (!element_) return false;

    return element_->DoesLineIntersect(I3Position(lineStart.GetX()-offset_.GetX(),
                                                  lineStart.GetY()-offset_.GetY(),
                                                  lineStart.GetZ()-offset_.GetZ()),
                                       I3Position(lineEnd.GetX()  -offset_.GetX(),
                                                  lineEnd.GetY()  -offset_.GetY(),
                                                  lineEnd.GetZ()  -offset_.GetZ()));
}

std::pair<I3Position, I3Position>
I3ExtraGeometryItemMove::GetBoundingBox
() const
{
    CalculateBoundingBox();
    return std::make_pair(boundingBoxLower_, boundingBoxUpper_);
}



template <class Archive>
void I3ExtraGeometryItemMove::load(Archive &ar, unsigned version)
{
    if (version > i3extrageometryitemmove_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3ExtraGeometryItem class.",version,i3extrageometryitemmove_version_);
    
    ar >> make_nvp("I3ExtraGeometryItem", base_object<I3ExtraGeometryItem>(*this));
    
    {
        I3ExtraGeometryItemPtr element_nonconst;
        ar >> make_nvp("element", element_nonconst);
        element_ = element_nonconst;
    }
}


template <class Archive>
void I3ExtraGeometryItemMove::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("I3ExtraGeometryItem", base_object<I3ExtraGeometryItem>(*this));
    ar << make_nvp("element", element_);
}


std::ostream& I3ExtraGeometryItemMove::operator<<(std::ostream& oss) const
{
    oss << "[ I3ExtraGeometryItemMove :" << std::endl;
    oss << "                   offset : (" << offset_.GetX()/I3Units::m << "," << offset_.GetY()/I3Units::m << "," << offset_.GetZ()/I3Units::m << ")m" << std::endl;

    oss << "  -> subitem: (not impl yet)" << std::endl;
    
    oss << "]" << std::endl;
    return oss;
}


I3_SPLIT_SERIALIZABLE(I3ExtraGeometryItemMove);
