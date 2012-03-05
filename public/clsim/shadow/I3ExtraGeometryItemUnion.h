/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ExtraGeometryItemUnion.h
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
    
    std::vector<I3ExtraGeometryItemConstPtr> elements_;
    
    friend class boost::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;

    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

BOOST_CLASS_VERSION(I3ExtraGeometryItemUnion, i3extrageometryitemunion_version_);

I3_POINTER_TYPEDEFS(I3ExtraGeometryItemUnion);

#endif //I3EXTRAGEOMETRYITEMUNION_H_INCLUDED
