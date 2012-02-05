//
//   Copyright (c) 2012  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module clsim
//
//   g4sim-intrface is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef I3CLSIMFLASHERPULSE_H_INCLUDED
#define I3CLSIMFLASHERPULSE_H_INCLUDED

#include "icetray/I3TrayHeaders.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/I3Position.h"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#define I3CLSIMFLASHERPULSE_H_I3CLSimFlasherPulse_FlasherPulseType \
    (Unknown)(LED340nm)(LED370nm)(LED405nm)(LED450nm)(LED505nm)    \
    (SC1)(SC2)

static const unsigned i3clsimflasherpulse_version_ = 0;

/**
 * @brief Describes a single flasher LED (or SC) pulse.
 * I3FlasherInfo is not really suited for this, so we keep
 * our own structure here. One I3FlasherInfo object can
 * describe multiple LEDs at the same time (using a bitmask),
 * which is something we don't want here. (These LEDs could
 * have different wavelengths, ...)
 */
class I3CLSimFlasherPulse
{
public:
    ~I3CLSimFlasherPulse();
    I3CLSimFlasherPulse();

    enum FlasherPulseType {
        Unknown  = 0,
        LED340nm = 1,
        LED370nm = 2,
        LED405nm = 3,
        LED450nm = 4,
        LED505nm = 5,
        SC1      = 6,
        SC2      = 7
    };

    inline FlasherPulseType GetType() const {return flasherPulseType_;}
    inline void SetType(FlasherPulseType newType) {flasherPulseType_ = newType;}
    
    
    inline const I3Position &GetPos() const {return pos_;}
    inline const I3Direction &GetDir() const {return dir_;}

    inline void SetPos(const I3Position &newPos) {pos_=newPos;}
    inline void SetDir(const I3Direction &newDir) {dir_=newDir;}


    inline double GetTime() const {return time_;}
    inline void SetTime(double newTime) {time_=newTime;}

    inline double GetNumberOfPhotonsNoBias() const {return numberOfPhotonsNoBias_;}
    inline void SetNumberOfPhotonsNoBias(double newNumber) {numberOfPhotonsNoBias_=newNumber;}

    inline double GetPulseWidth() const {return pulseWidth_;}
    inline void SetPulseWidth(double newWidth) {pulseWidth_=newWidth;}
    
    
private:
    FlasherPulseType flasherPulseType_;

    I3Direction dir_;
    I3Position pos_;
    double time_;
    
    double numberOfPhotonsNoBias_;
    double pulseWidth_; // FWHM

private: // static stuff
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);

    friend bool operator==(const I3CLSimFlasherPulse &, const I3CLSimFlasherPulse &);
};
bool operator==(const I3CLSimFlasherPulse &a, const I3CLSimFlasherPulse &b);

BOOST_CLASS_VERSION(I3CLSimFlasherPulse, i3clsimflasherpulse_version_);

typedef I3Vector<I3CLSimFlasherPulse> I3CLSimFlasherPulseSeries;

I3_POINTER_TYPEDEFS(I3CLSimFlasherPulse);
I3_POINTER_TYPEDEFS(I3CLSimFlasherPulseSeries);

#endif //I3CLSIMFLASHERPULSE_H_INCLUDED
