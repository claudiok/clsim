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

#include <clsim/I3CLSimLightSource.h>

#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>



// construction & destruction
I3CLSimLightSource::I3CLSimLightSource(const I3Particle &particle)
:
lightSourceType_(I3CLSimLightSource::Particle),
particle_(particle)
{
    
}

I3CLSimLightSource::I3CLSimLightSource(const I3CLSimFlasherPulse &flasher)
:
lightSourceType_(I3CLSimLightSource::Flasher),
flasher_(flasher)
{

}

I3CLSimLightSource::I3CLSimLightSource()
:
lightSourceType_(I3CLSimLightSource::Unknown)
{
    
}

I3CLSimLightSource::I3CLSimLightSource(const I3CLSimLightSource &lightSource)
:
lightSourceType_(lightSource.lightSourceType_)
{
    // copy either particle or flasher info depending on type
    switch(lightSourceType_)
    {
        case Unknown:
            break;
        case Particle:
            particle_ = lightSource.particle_;
            break;
        case Flasher:
            flasher_ = lightSource.flasher_;
            break;
        default:
            break;
    }
}


I3CLSimLightSource::~I3CLSimLightSource() 
{ 
    
}


// comparison
bool operator==(const I3CLSimLightSource &a, const I3CLSimLightSource &b)
{
    if (a.lightSourceType_ != b.lightSourceType_) return false;
    
    if (a.lightSourceType_ == I3CLSimLightSource::Unknown) return true; // nothing else to compare

    if (a.lightSourceType_ == I3CLSimLightSource::Particle) {
        return (a.particle_==b.particle_);
    } else if (a.lightSourceType_ == I3CLSimLightSource::Flasher) {
        return (a.flasher_==b.flasher_);
    } else {
        throw std::runtime_error("I3CLSimLightSource has an undefined type!");
    }

    return true;
}



const I3Particle &I3CLSimLightSource::GetParticle() const
{
    if (lightSourceType_ != Particle)
        throw std::runtime_error("this light source is not an I3Particle");

    return particle_;
}


const I3CLSimFlasherPulse &I3CLSimLightSource::GetFlasherPulse() const
{
    if (lightSourceType_ != Flasher)
        throw std::runtime_error("this light source is not an I3CLSimFlasherPulse");
    
    return flasher_;
}

