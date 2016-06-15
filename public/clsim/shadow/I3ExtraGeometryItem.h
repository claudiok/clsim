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
 * @file I3ExtraGeometryItem.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
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
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

inline std::ostream& operator<<(std::ostream& oss, const I3ExtraGeometryItem &item)
{
    return item.operator<<(oss);
}

I3_CLASS_VERSION(I3ExtraGeometryItem, i3extrageometryitem_version_);

I3_POINTER_TYPEDEFS(I3ExtraGeometryItem);

#endif //I3EXTRAGEOMETRYITEM_H_INCLUDED
