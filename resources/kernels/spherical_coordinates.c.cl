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
 * @file spherical_coordinates.c.cl
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
    floating_t n_rho = magnitude(rho);
    
    // radius
    coords.s0 = magnitude(pos);
    // azimuth
    floating_t azimuth = (n_rho > 0) ? acos(dot(rho, source->perpDir)/n_rho)/(PI/180) : 0;
#ifdef HAS_FULL_AZIMUTH_EXTENSION
        // need to determine direction of rho in case table has an azimuth extension up to 360 deg
        floating4_t azisignvec = cross(rho, source->perpDir);
        floating_t azisign = dot(azisignvec, source->dir);
        coords.s1 = (azisign > 0) ? 360.-azimuth : azimuth;
#else
        coords.s1 = azimuth;
#endif
    // cos(polar angle)
    coords.s2 = (coords.s0 > 0) ? my_divide(l, coords.s0) : 0;
    // delay time
    coords.s3 = pos.w - coords.s0*min_invGroupVel;
#ifdef TABULATE_IMPACT_ANGLE
    // s4 is the cosine of the opening angle between a vector connecting
    // the DOM's center to the photon impact point and a vector connecting
    // the center to the emitter. Here we average over possible DOM positions
    // by randomizing the impact position in the cross-sectional area of the
    // DOM. Note that because the impact parameter is expressed as a 
    // rotation across the surface of the DOM, it is independent of the DOM
    // radius.
    floating_t sina = my_sqrt(RNG_CALL_UNIFORM_CO);
    scatterDirectionByAngle(my_sqrt(1-sina*sina), sina, &dirAndWlen, RNG_CALL_UNIFORM_CO);
    coords.s4 = (coords.s0 > 0) ? my_divide(dot(dirAndWlen, pos), coords.s0) : 1;
#endif
    
    //dbg_printf("     %4.1f %4.1f %4.2f %6.2f\n", coords.s0, coords.s1, coords.s2, coords.s3);
    
    return coords;
}
