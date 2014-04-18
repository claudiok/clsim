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
 * @file I3ExtraGeometryItemUnion.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <limits>

#include <icetray/serialization.h>
#include <clsim/shadow/I3ExtraGeometryItemUnion.h>

#include <boost/foreach.hpp>

I3ExtraGeometryItemUnion::I3ExtraGeometryItemUnion()
:
boundingBoxCalculated_(false)
{;}

I3ExtraGeometryItemUnion::
I3ExtraGeometryItemUnion(const std::vector<I3ExtraGeometryItemConstPtr> &elements)
:
elements_(elements),
boundingBoxCalculated_(false)
{
}

I3ExtraGeometryItemUnion::~I3ExtraGeometryItemUnion() { }

void I3ExtraGeometryItemUnion::CalculateBoundingBox() const
{
    if (boundingBoxCalculated_) return;
    
    double lowX=std::numeric_limits<double>::infinity();
    double lowY=std::numeric_limits<double>::infinity();
    double lowZ=std::numeric_limits<double>::infinity();
    double highX=-std::numeric_limits<double>::infinity();
    double highY=-std::numeric_limits<double>::infinity();
    double highZ=-std::numeric_limits<double>::infinity();
    
    BOOST_FOREACH(const I3ExtraGeometryItemConstPtr &ptr, elements_)
    {
        if (!ptr) continue;
        
        const std::pair<I3Position, I3Position> box =
        ptr->GetBoundingBox();
        
        if (box.first.GetX()  < lowX)  lowX  = box.first.GetX();
        if (box.second.GetX() < lowX)  lowX  = box.second.GetX();
        if (box.first.GetX()  > highX) highX = box.first.GetX();
        if (box.second.GetX() > highX) highX = box.second.GetX();

        if (box.first.GetY()  < lowY)  lowY  = box.first.GetY();
        if (box.second.GetY() < lowY)  lowY  = box.second.GetY();
        if (box.first.GetY()  > highY) highY = box.first.GetY();
        if (box.second.GetY() > highY) highY = box.second.GetY();

        if (box.first.GetZ()  < lowZ)  lowZ  = box.first.GetZ();
        if (box.second.GetZ() < lowZ)  lowZ  = box.second.GetZ();
        if (box.first.GetZ()  > highZ) highZ = box.first.GetZ();
        if (box.second.GetZ() > highZ) highZ = box.second.GetZ();
    }
    
    boundingBoxLower_=I3Position(lowX, lowY, lowZ);
    boundingBoxUpper_=I3Position(highX, highY, highZ);
    
    boundingBoxCalculated_=true;
}


bool
I3ExtraGeometryItemUnion::DoesLineIntersect
(const I3Position &lineStart,
 const I3Position &lineEnd) const
{
    BOOST_FOREACH(const I3ExtraGeometryItemConstPtr &ptr, elements_)
    {
        if (!ptr) continue;
        if (ptr->DoesLineIntersect(lineStart, lineEnd)) return true;
    }
    
    return false;
}

std::pair<I3Position, I3Position>
I3ExtraGeometryItemUnion::GetBoundingBox
() const
{
    CalculateBoundingBox();
    
    return std::make_pair(boundingBoxLower_, boundingBoxUpper_);
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

// TODO: the following line needs to go away once a more recent version of icetray makes it into a release:
#include "split_serializable_backport.h"

I3_SPLIT_SERIALIZABLE(I3ExtraGeometryItemUnion);
