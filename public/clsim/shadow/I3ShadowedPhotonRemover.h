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
 * @file I3ShadowedPhotonRemover.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3SHADOWEDPHOTONREMOVER_H_INCLUDED
#define I3SHADOWEDPHOTONREMOVER_H_INCLUDED

#include "simclasses/I3CompressedPhoton.h"

#include <string>

#include "simclasses/I3ExtraGeometryItem.h"

#include "dataclasses/geometry/I3OMGeo.h"

#include "simclasses/I3CylinderMap.h"

/**
 * @brief This class checks if a photon path intersects with 
 *   any shadowing part of the detecor (such as cables).
 *   This code is NOT functional at the moment.
 */
class I3ShadowedPhotonRemover
{
public:
  I3ShadowedPhotonRemover(const I3CylinderMap &cylinder_map , const double &distance);
    ~I3ShadowedPhotonRemover();
    
    
    /**
     * returns true if the photon hits any of the extra geometry
     */
    bool IsPhotonShadowed(const I3CompressedPhoton &photon);
    double direction_azimuth;
    double direction_zenith;
    double dx;
    double dy;
    double dz;
    I3Position start_position;

private:
    const I3CylinderMap& cylinder_map_;
    const double& distance_;
    // parameters

    SET_LOGGER("I3ShadowedPhotonRemover");
};

I3_POINTER_TYPEDEFS(I3ShadowedPhotonRemover);

#endif //I3SHADOWEDPHOTONREMOVER_H_INCLUDED
