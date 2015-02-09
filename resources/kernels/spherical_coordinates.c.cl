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

typedef floating4_t coordinate_t;

inline floating_t magnitude(floating4_t vec)
{
    return my_sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

inline coordinate_t
getCoordinates(const floating4_t absPos, const struct I3CLSimReferenceParticle *source)
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
    coords.s1 = (n_rho > 0) ?
        acos(-dot(rho,source->perpDir)/n_rho)/(PI/180) : 0;
    // cos(polar angle)
    coords.s2 = (coords.s0 > 0) ? my_divide(l, coords.s0) : 0;
    // delay time
    coords.s3 = pos.w - coords.s0*min_invPhaseVel;
    
    dbg_printf("     %4.1f %4.1f %4.2f %6.2f\n", coords.s0, coords.s1, coords.s2, coords.s3);
    
    return coords;
}
