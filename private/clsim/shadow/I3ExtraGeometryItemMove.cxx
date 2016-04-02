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
 * @file I3ExtraGeometryItemMove.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
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
        boundingBoxLower_=I3Position();
        boundingBoxUpper_=I3Position();
        
        boundingBoxCalculated_=true;
        return;
    }
    
    const std::pair<I3Position, I3Position> box =
    element_->GetBoundingBox();
    
    boundingBoxLower_=box.first+offset_;
    boundingBoxUpper_=box.second+offset_;
    
    boundingBoxCalculated_=true;
}


bool
I3ExtraGeometryItemMove::DoesLineIntersect
(const I3Position &lineStart,
 const I3Position &lineEnd) const
{
    if (!element_) return false;

    return element_->DoesLineIntersect(lineStart-offset_,
                                       lineEnd-offset_);
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

// TODO: the following line needs to go away once a more recent version of icetray makes it into a release:
#include "split_serializable_backport.h"

I3_SPLIT_SERIALIZABLE(I3ExtraGeometryItemMove);
