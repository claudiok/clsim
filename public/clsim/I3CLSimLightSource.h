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
 * @file I3CLSimLightSource.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMLIGHTSOURCE_H_INCLUDED
#define I3CLSIMLIGHTSOURCE_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include <dataclasses/physics/I3Particle.h>
#include <clsim/I3CLSimFlasherPulse.h>

#define I3CLSIMLIGHTSOURCE_H_I3CLSimLightSource_LightSourceType \
    (Unknown)(Particle)(Flasher)

/**
 * @brief Describes a light source. This is currently
 * either a I3Particle or a flasher described by
 * I3CLSimFlasherPulse and some additional variables.
 */
class I3CLSimLightSource
{
public:
    ~I3CLSimLightSource();

    // copy construction:
    I3CLSimLightSource(const I3CLSimLightSource &lightSource);

    I3CLSimLightSource(const I3Particle &particle);
    I3CLSimLightSource(const I3CLSimFlasherPulse &flasher);
    
    enum LightSourceType {
        Unknown = 0,
        Particle = 1,
        Flasher = 2
    };

    inline LightSourceType GetType() const {return lightSourceType_;}
    
    const I3Particle &GetParticle() const;
    const I3CLSimFlasherPulse &GetFlasherPulse() const;

private:
    I3CLSimLightSource(); // no default construction
    
    LightSourceType lightSourceType_;

    I3Particle particle_;
    I3CLSimFlasherPulse flasher_;

private: // static stuff
    friend bool operator==(const I3CLSimLightSource &, const I3CLSimLightSource &);
};
bool operator==(const I3CLSimLightSource &a, const I3CLSimLightSource &b);

typedef std::vector<I3CLSimLightSource> I3CLSimLightSourceSeries;

I3_POINTER_TYPEDEFS(I3CLSimLightSource);
I3_POINTER_TYPEDEFS(I3CLSimLightSourceSeries);

#endif //I3CLSIMLIGHTSOURCE_H_INCLUDED
