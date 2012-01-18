/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3HoleIceSimulator.h
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

#ifndef I3HOLEICESIMULATOR_H_INCLUDED
#define I3HOLEICESIMULATOR_H_INCLUDED

#include "phys-services/I3RandomService.h"

#include "dataclasses/I3Position.h"

#include "clsim/I3Photon.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimWlenDependentValue.h"

#include <string>


/**
 * @brief This module propagates photons on the surface of
 * an oversized DOM to the actual DOM surface. It can take
 * hole ice and cable shadow effects into account.
 *
 */
class I3HoleIceSimulator
{
public:
    I3HoleIceSimulator(I3RandomServicePtr random,
                       double DOMRadius,
                       double holeRadius,
                       I3CLSimMediumPropertiesConstPtr mediumProperties,
                       I3CLSimWlenDependentValueConstPtr holeIceAbsorptionLength,
                       I3CLSimWlenDependentValueConstPtr holeIceScatteringLength);
    ~I3HoleIceSimulator();
    
    
    /**
     * Tracks an I3Photon in place
     */
    bool TrackPhoton(I3Photon &photon,
                     const I3Position &om_pos);
    

    
private:
    // parameters
    I3RandomServicePtr random_;
    double DOMRadius_;
    double holeRadius_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    I3CLSimWlenDependentValueConstPtr holeIceAbsorptionLength_;
    I3CLSimWlenDependentValueConstPtr holeIceScatteringLength_;
    
private:
    // default, assignment, and copy constructor declared private
    I3HoleIceSimulator();
    I3HoleIceSimulator(const I3HoleIceSimulator&);
    I3HoleIceSimulator& operator=(const I3HoleIceSimulator&);
    
    SET_LOGGER("I3HoleIceSimulator");
};

I3_POINTER_TYPEDEFS(I3HoleIceSimulator);

#endif //I3HOLEICESIMULATOR_H_INCLUDED
