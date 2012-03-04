/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ShadowedPhotonRemoverModule.h
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

#ifndef I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED
#define I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "dataclasses/geometry/I3Geometry.h"

#include "clsim/I3ShadowedPhotonRemover.h"

#include <string>


/**
 * @brief This module removes photons that have paths intersecting 
 *   with any shadowing part of the detecor (such as cables).
 *
 */
class I3ShadowedPhotonRemoverModule : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3ShadowedPhotonRemoverModule(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    ~I3ShadowedPhotonRemoverModule();
    
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
    
    /// Parameter: Name of the input I3PhotonSeriesMap frame object. 
    std::string inputPhotonSeriesMapName_;

    /// Parameter: Name of the output I3PhotonSeriesMap frame object. 
    std::string outputPhotonSeriesMapName_;

    I3ShadowedPhotonRemoverPtr shadowedPhotonRemover_;
    
private:
    // default, assignment, and copy constructor declared private
    I3ShadowedPhotonRemoverModule();
    I3ShadowedPhotonRemoverModule(const I3ShadowedPhotonRemoverModule&);
    I3ShadowedPhotonRemoverModule& operator=(const I3ShadowedPhotonRemoverModule&);

    SET_LOGGER("I3ShadowedPhotonRemoverModule");
};

#endif //I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED
