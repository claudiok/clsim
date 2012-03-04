/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3CLSimPhotonHistory.h
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

#ifndef I3CLSIMPHOTONHISTORY_H_INCLUDED
#define I3CLSIMPHOTONHISTORY_H_INCLUDED

#include <cstring>

#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Position.h"


/**
 * @brief Stores the history of a photon as the last
 * points of scatter before reaching a DOM.
 */
static const unsigned i3clsimphotonhistory_version_ = 0;

struct I3CLSimPhotonHistory
{
    
public:
    
    I3CLSimPhotonHistory() {;}
    
    ~I3CLSimPhotonHistory();

    inline std::size_t size() const {return posX_.size();}
    
    inline float GetX(std::size_t index) const {return posX_[index];}
    inline float GetY(std::size_t index) const {return posY_[index];}
    inline float GetZ(std::size_t index) const {return posZ_[index];}
    
    inline I3PositionPtr operator[](std::size_t index) const {return I3PositionPtr(new I3Position(posX_[index], posY_[index], posZ_[index]));}
    inline I3PositionPtr at(std::size_t index) const {return I3PositionPtr(new I3Position(posX_.at(index), posY_.at(index), posZ_.at(index)));}
    
    inline void push_back(float x, float y, float z) {posX_.push_back(x); posY_.push_back(y); posZ_.push_back(z);} 
    inline void push_back(const I3Position &pos) {push_back(pos.GetX(), pos.GetY(), pos.GetZ());} 
    
private:
    std::vector<float> posX_;
    std::vector<float> posY_;
    std::vector<float> posZ_;
    
private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

inline bool operator==(const I3CLSimPhotonHistory &a, const I3CLSimPhotonHistory &b) 
{
    if (a.size() != b.size()) return false;
    
    for (std::size_t i=0;i<a.size();++i)
    {
        if (a.GetX(i) != b.GetX(i)) return false;
        if (a.GetY(i) != b.GetY(i)) return false;
        if (a.GetZ(i) != b.GetZ(i)) return false;
    }
    
    return true;
}

BOOST_CLASS_VERSION(I3CLSimPhotonHistory, i3clsimphotonhistory_version_);

typedef I3Vector<I3CLSimPhotonHistory> I3CLSimPhotonHistorySeries;

I3_POINTER_TYPEDEFS(I3CLSimPhotonHistory);
I3_POINTER_TYPEDEFS(I3CLSimPhotonHistorySeries);

#endif //I3CLSIMPHOTONHISTORY_H_INCLUDED
