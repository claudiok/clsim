/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ExtraGeometryItem.h
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

#ifndef I3EXTRAGEOMETRYITEM_H_INCLUDED
#define I3EXTRAGEOMETRYITEM_H_INCLUDED


#include <cstring>

#include "icetray/I3FrameObject.h"

#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Position.h"


/**
 * @brief Abstract base class for a "extra" geometry item
 *  (i.e. anything that is not a DOM). This includes cables
 * and anything else that could cast a shadow on DOMs.
 */
static const unsigned i3extrageometryitem_version_ = 0;

struct I3ExtraGeometryItem : public I3FrameObject
{
public:
    I3ExtraGeometryItem() {;}
    virtual ~I3ExtraGeometryItem();

    /**
     * Should return true if the line intersects the item.
     */
    virtual bool DoesLineIntersect(const I3Position &lineStart,
                           const I3Position &lineEnd) const = 0;

    /**
     * Returns two corners of the object's bounding box.
     */
    virtual std::pair<I3Position, I3Position> GetBoundingBox() const = 0;
    
    virtual std::ostream& operator<<(std::ostream& oss) const;

private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

inline std::ostream& operator<<(std::ostream& oss, const I3ExtraGeometryItem &item)
{
    return item.operator<<(oss);
}

BOOST_CLASS_VERSION(I3ExtraGeometryItem, i3extrageometryitem_version_);

I3_POINTER_TYPEDEFS(I3ExtraGeometryItem);

#endif //I3EXTRAGEOMETRYITEM_H_INCLUDED
