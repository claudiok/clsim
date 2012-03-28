/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ShadowedPhotonRemover.h
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

#ifndef I3SHADOWEDPHOTONREMOVER_H_INCLUDED
#define I3SHADOWEDPHOTONREMOVER_H_INCLUDED

#include "clsim/I3Photon.h"

#include <string>


/**
 * @brief This class checks if a photon path intersects with 
 *   any shadowing part of the detecor (such as cables).
 *   This code is NOT functional at the moment.
 */
class I3ShadowedPhotonRemover
{
public:
    I3ShadowedPhotonRemover();
    ~I3ShadowedPhotonRemover();
    
    
    /**
     * returns true if the photon hits any of the extra geometry
     */
    bool IsPhotonShadowed(const I3Photon &photon) const;
    

    
private:
    // parameters
    
private:
    // assignment and copy constructor declared private
    //I3ShadowedPhotonRemover();
    I3ShadowedPhotonRemover(const I3ShadowedPhotonRemover&);
    I3ShadowedPhotonRemover& operator=(const I3ShadowedPhotonRemover&);
    
    SET_LOGGER("I3ShadowedPhotonRemover");
};

I3_POINTER_TYPEDEFS(I3ShadowedPhotonRemover);

#endif //I3SHADOWEDPHOTONREMOVER_H_INCLUDED
