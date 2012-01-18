/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3HoleIceSimulatorModule.h
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

#ifndef I3HOLEICESIMULATORMODULE_H_INCLUDED
#define I3HOLEICESIMULATORMODULE_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "dataclasses/geometry/I3Geometry.h"

#include "phys-services/I3RandomService.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimWlenDependentValue.h"

#include "clsim/I3HoleIceSimulator.h"

#include <string>


/**
 * @brief This module propagates photons on the surface of
 * an oversized DOM to the actual DOM surface. It can take
 * hole ice and cable shadow effects into account.
 *
 */
class I3HoleIceSimulatorModule : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3HoleIceSimulatorModule(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    ~I3HoleIceSimulatorModule();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
    /**
     * The module needs to process Physics frames
     */
#ifdef IS_Q_FRAME_ENABLED
    void DAQ(I3FramePtr frame);
#else
    void Physics(I3FramePtr frame);
#endif

    
private:
    // parameters
    
    /// Parameter: A random number generating service (derived from I3RandomService).
    I3RandomServicePtr randomService_;

    /// Parameter: Name of the input I3PhotonSeriesMap frame object. 
    std::string inputPhotonSeriesMapName_;

    /// Parameter: Name of the output I3PhotonSeriesMap frame object. 
    std::string outputPhotonSeriesMapName_;

    /// Parameter: Specifiy the DOM radius. Do not include oversize factors here.
    double DOMRadiusWithoutOversize_;

    /// Parameter: Specifiy the hole radius.
    double holeRadius_;

    /// Parameter: An instance of I3CLSimMediumProperties describing the ice/water properties.
    I3CLSimMediumPropertiesConstPtr mediumProperties_;

    /// Parameter: An instance of I3CLSimWlenDependentValue describing the hole ice absorption lengths.
    I3CLSimWlenDependentValueConstPtr holeIceAbsorptionLength_;

    /// Parameter: An instance of I3CLSimWlenDependentValue describing the hole ice scattering lengths.
    I3CLSimWlenDependentValueConstPtr holeIceScatteringLength_;

    
    I3HoleIceSimulatorPtr holeIceSimulator_;
    
private:
    // default, assignment, and copy constructor declared private
    I3HoleIceSimulatorModule();
    I3HoleIceSimulatorModule(const I3HoleIceSimulatorModule&);
    I3HoleIceSimulatorModule& operator=(const I3HoleIceSimulatorModule&);

    SET_LOGGER("I3HoleIceSimulatorModule");
};

#endif //I3HOLEICESIMULATORMODULE_H_INCLUDED
