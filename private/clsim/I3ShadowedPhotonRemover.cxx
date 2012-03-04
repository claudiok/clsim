/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3ShadowedPhotonRemover.cxx
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

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif 
#include <inttypes.h>

#include <limits>

#include "clsim/I3ShadowedPhotonRemover.h"

#include "dataclasses/I3Constants.h"

I3ShadowedPhotonRemover::I3ShadowedPhotonRemover() //(I3ExtraGeometry extraGeometry) 
//:
//extraGeometry_(extraGeometry)
{
    log_trace("%s", __PRETTY_FUNCTION__);

}

I3ShadowedPhotonRemover::~I3ShadowedPhotonRemover()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}



bool I3ShadowedPhotonRemover::IsPhotonShadowed(const I3Photon &photon) const
{

    
    return false;
}
