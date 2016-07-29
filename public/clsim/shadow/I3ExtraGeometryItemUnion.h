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
 * @file I3ExtraGeometryItemUnion.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3EXTRAGEOMETRYITEMUNION_H_INCLUDED
#define I3EXTRAGEOMETRYITEMUNION_H_INCLUDED

#include <vector>

#include "clsim/shadow/I3ExtraGeometryItem.h"

/**
 * @brief Describes a union af a list of items.
 */
static const unsigned i3extrageometryitemunion_version_ = 0;

struct I3ExtraGeometryItemUnion : public I3ExtraGeometryItem
{
public:
    virtual ~I3ExtraGeometryItemUnion();

    I3ExtraGeometryItemUnion(const std::vector<I3ExtraGeometryItemConstPtr> &elements);
    I3ExtraGeometryItemUnion();
    
    virtual bool DoesLineIntersect(const I3Position &lineStart,
                           const I3Position &lineEnd) const;
    virtual std::pair<I3Position, I3Position> GetBoundingBox() const;

    virtual std::ostream& operator<<(std::ostream& oss) const;

private:
    void CalculateBoundingBox() const;

    std::vector<I3ExtraGeometryItemConstPtr> elements_;

    mutable bool boundingBoxCalculated_;
    mutable I3Position boundingBoxLower_;
    mutable I3Position boundingBoxUpper_;
    
    friend class icecube::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;

    I3_SERIALIZATION_SPLIT_MEMBER();
};

I3_CLASS_VERSION(I3ExtraGeometryItemUnion, i3extrageometryitemunion_version_);

I3_POINTER_TYPEDEFS(I3ExtraGeometryItemUnion);

#endif //I3EXTRAGEOMETRYITEMUNION_H_INCLUDED
