/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ExtraGeometryItemMove.h
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

#ifndef I3EXTRAGEOMETRYITEMMOVE_H_INCLUDED
#define I3EXTRAGEOMETRYITEMMOVE_H_INCLUDED

#include <vector>

#include "clsim/shadow/I3ExtraGeometryItem.h"

/**
 * @brief Describes a union af a list of items.
 */
static const unsigned i3extrageometryitemmove_version_ = 0;

struct I3ExtraGeometryItemMove : public I3ExtraGeometryItem
{
public:
    virtual ~I3ExtraGeometryItemMove();

    I3ExtraGeometryItemMove(I3ExtraGeometryItemConstPtr element, const I3Position &offset);
    I3ExtraGeometryItemMove();
    
    virtual bool DoesLineIntersect(const I3Position &lineStart,
                           const I3Position &lineEnd) const;
    virtual std::pair<I3Position, I3Position> GetBoundingBox() const;

    virtual std::ostream& operator<<(std::ostream& oss) const;

private:

    I3ExtraGeometryItemConstPtr element_;
    I3Position offset_;
    
    void CalculateBoundingBox() const;
    mutable bool boundingBoxCalculated_;
    mutable I3Position boundingBoxLower_;
    mutable I3Position boundingBoxUpper_;
    
    friend class boost::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;

    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

BOOST_CLASS_VERSION(I3ExtraGeometryItemMove, i3extrageometryitemmove_version_);

I3_POINTER_TYPEDEFS(I3ExtraGeometryItemMove);

#endif //I3ExtraGeometryItemMove_H_INCLUDED
