/**
 * Copyright (c) 2015
 * Jakob van Santen <jvansanten@icecube.wisc.edu>
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
 * @file cylindrical_coordinates.c.cl
 * @version $LastChangedRevision$
 * @date $Date$
 * @author Jakob van Santen
 */

#ifdef TABULATE_IMPACT_ANGLE
typedef float8 coordinate_t;
#else
typedef float4 coordinate_t;
#endif

inline floating_t magnitude(floating4_t vec)
{
    return my_sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

// Coordinate transform for infinite muon tracks. Here the source depth is
// degenerate with time, so we have 1 fewer dimension
inline coordinate_t
getCoordinates(const floating4_t absPos, floating4_t dirAndWlen,
    const struct I3CLSimReferenceParticle *source, RNG_ARGS)
{
    coordinate_t coords;
    
    // NB: the reference vectors sourcePos, sourceDir, and perpDir are
    //     defined as static variables at compile time
    floating4_t pos = absPos - source->posAndTime;
    floating_t l = dot(pos, source->dir);
    floating4_t rho = pos - l*source->dir;
    
    // perpendicular distance
    coords.s0 = magnitude(rho);
    // azimuth (in radians, because that's how the first implementation worked)
    coords.s1 = (coords.s0 > 0) ?
        acos(dot(rho,source->perpDir)/coords.s0) : 0;
    // depth of closest approach
    coords.s2 = source->posAndTime.z + l*source->dir.z;
    // delay time
    coords.s3 = pos.w - (l + coords.s0*tan_thetaC)*recip_speedOfLight;
#ifdef TABULATE_IMPACT_ANGLE
    // s4 is the cosine of the opening angle between a vector connecting
    // the DOM's center to the photon impact point and a vector connecting
    // the center to the nominal Cherenkov emission point. Here we average over possible DOM positions
    // by randomizing the impact position in the cross-sectional area of the
    // DOM. Note that because the impact parameter is expressed as a 
    // rotation across the surface of the DOM, it is independent of the DOM
    // radius.
    floating_t sina = my_sqrt(RNG_CALL_UNIFORM_CO);
    scatterDirectionByAngle(my_sqrt(1-sina*sina), sina, &dirAndWlen, RNG_CALL_UNIFORM_CO);
    // find the vector connecting the center of the DOM to the nominal Cherenkov
    // emission point (i.e. perpendicular to the Cherenkov light front)
    floating4_t cpos = absPos - (source->posAndTime + (l-rho*my_recip(tan_thetaC))*source->dir);
    floating_t cdist = magnitude(cpos);
    coords.s4 = (cdist > 0) ? my_divide(dot(dirAndWlen, cpos), cdist) : 1;
#endif
    
    return coords;
}
