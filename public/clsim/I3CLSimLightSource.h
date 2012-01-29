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

#ifndef I3CLSIMLIGHTSOURCE_H_INCLUDED
#define I3CLSIMLIGHTSOURCE_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3FlasherInfo.h>

#define I3CLSIMLIGHTSOURCE_H_I3CLSimLightSource_LightSourceType \
    (Unknown)(Particle)(Flasher)

/**
 * @brief Describes a light source. This is currently
 * either a I3Particle or a flasher described by
 * I3FlasherInfo and some additional variables.
 */
class I3CLSimLightSource
{
public:
    ~I3CLSimLightSource();

    I3CLSimLightSource(const I3Particle &particle);
    I3CLSimLightSource(const I3FlasherInfo &flasher);
    
    enum LightSourceType {
        Unknown = 0,
        Particle = 1,
        Flasher = 2
    };

    inline LightSourceType GetType() const {return lightSourceType_;}
    
    const I3Particle &GetParticle() const;
    const I3FlasherInfo &GetFlasherInfo() const;

private:
    I3CLSimLightSource(); // no default construction
    
    LightSourceType lightSourceType_;

    I3Particle particle_;
    I3FlasherInfo flasher_;

private: // static stuff
    friend bool operator==(const I3CLSimLightSource &, const I3CLSimLightSource &);
};
bool operator==(const I3CLSimLightSource &a, const I3CLSimLightSource &b);

typedef std::vector<I3CLSimLightSource> I3CLSimLightSourceSeries;

I3_POINTER_TYPEDEFS(I3CLSimLightSource);
I3_POINTER_TYPEDEFS(I3CLSimLightSourceSeries);

#endif //I3CLSIMLIGHTSOURCE_H_INCLUDED
