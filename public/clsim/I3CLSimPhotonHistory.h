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
 * @file I3CLSimPhotonHistory.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
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
static const unsigned i3clsimphotonhistory_version_ = 1;

struct I3CLSimPhotonHistory
{
    
public:
    
    I3CLSimPhotonHistory() {;}
    
    ~I3CLSimPhotonHistory();

    inline std::size_t size() const {return posX_.size();}
    
    inline float GetX(std::size_t index) const {return posX_[index];}
    inline float GetY(std::size_t index) const {return posY_[index];}
    inline float GetZ(std::size_t index) const {return posZ_[index];}
    inline float GetDistanceInAbsorptionLengths(std::size_t index) const {return distanceInAbsorptionLengths_[index];}
    
    inline I3PositionPtr operator[](std::size_t index) const {return I3PositionPtr(new I3Position(posX_[index], posY_[index], posZ_[index]));}
    inline I3PositionPtr at(std::size_t index) const {return I3PositionPtr(new I3Position(posX_.at(index), posY_.at(index), posZ_.at(index)));}
    
    inline void push_back(float x, float y, float z, float abslens) {posX_.push_back(x); posY_.push_back(y); posZ_.push_back(z); distanceInAbsorptionLengths_.push_back(abslens);} 
    inline void push_back(const I3Position &pos, float abslens) {push_back(pos.GetX(), pos.GetY(), pos.GetZ(), abslens);} 
    
private:
    std::vector<float> posX_;
    std::vector<float> posY_;
    std::vector<float> posZ_;
    std::vector<float> distanceInAbsorptionLengths_;
    
private:
    friend class icecube::serialization::access;
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
        if (a.GetDistanceInAbsorptionLengths(i) != b.GetDistanceInAbsorptionLengths(i)) return false;
    }
    
    return true;
}

I3_CLASS_VERSION(I3CLSimPhotonHistory, i3clsimphotonhistory_version_);

std::ostream& operator<<(std::ostream&, const I3CLSimPhotonHistory&);

typedef I3Vector<I3CLSimPhotonHistory> I3CLSimPhotonHistorySeries;

I3_POINTER_TYPEDEFS(I3CLSimPhotonHistory);
I3_POINTER_TYPEDEFS(I3CLSimPhotonHistorySeries);

#endif //I3CLSIMPHOTONHISTORY_H_INCLUDED
