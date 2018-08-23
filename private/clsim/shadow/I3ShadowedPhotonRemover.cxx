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
 * @file I3ShadowedPhotonRemover.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif 
#include <inttypes.h>

#include <limits>

#include "simclasses/I3Photon.h"

#include "clsim/shadow/I3ShadowedPhotonRemover.h"

#include "dataclasses/I3Constants.h"

#include "clsim/shadow/I3ExtraGeometryItem.h"

#include <cmath>

I3ShadowedPhotonRemover::I3ShadowedPhotonRemover(const I3ExtraGeometryItem &cylinder , const double &distance) : cylinder_(cylinder) , distance_(distance) //(I3ExtraGeometry extraGeometry) 
//:
//extraGeometry_(extraGeometry)
{
  log_trace("%s", __PRETTY_FUNCTION__);

}

I3ShadowedPhotonRemover::~I3ShadowedPhotonRemover()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


//This is a boolean that will return if the photon hits the cable or cylinder
bool I3ShadowedPhotonRemover::IsPhotonShadowed(const I3Photon &photon) 
{
  direction_azimuth = photon.GetDir().GetAzimuth();
  direction_zenith = photon.GetDir().GetZenith();
  dx = distance_ * sin ( direction_zenith ) * cos ( direction_azimuth );
  dy = distance_ * sin ( direction_zenith ) * cos ( direction_azimuth );
  dz = distance_ * cos ( direction_zenith );
  start_position = photon.GetPos() + I3Position(dx , dy , dz);
  return cylinder_.DoesLineIntersect( start_position , photon.GetPos());
}
