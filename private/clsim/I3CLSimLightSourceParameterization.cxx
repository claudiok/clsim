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
 * @file I3CLSimLightSourceParameterization.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/I3CLSimLightSourceParameterization.h>

#include <limits>
#include <cmath>

const I3CLSimLightSourceParameterization::AllParticles_t I3CLSimLightSourceParameterization::AllParticles = I3CLSimLightSourceParameterization::AllParticles_t();

I3CLSimLightSourceParameterization::I3CLSimLightSourceParameterization()
:
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
forPdgEncoding(0),
#else
forParticleType(I3Particle::unknown),
#endif
fromEnergy(0.),
toEnergy(std::numeric_limits<double>::infinity()),
needsLength(false),
catchAll(false),
flasherMode(false),
forFlasherPulseType(I3CLSimFlasherPulse::Unknown)
{ 
    
}

I3CLSimLightSourceParameterization::~I3CLSimLightSourceParameterization() 
{ 

}

I3CLSimLightSourceParameterization::I3CLSimLightSourceParameterization
(I3CLSimLightSourceToStepConverterPtr converter_,
 I3Particle::ParticleType forParticleType_,
 double fromEnergy_, double toEnergy_,
 bool needsLength_
 )
:
converter(converter_),
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
forPdgEncoding(I3Particle(I3Particle::Null,forParticleType_).GetPdgEncoding()),
#else
forParticleType(forParticleType_),
#endif
fromEnergy(fromEnergy_),
toEnergy(toEnergy_),
needsLength(needsLength_),
catchAll(false),
flasherMode(false),
forFlasherPulseType(I3CLSimFlasherPulse::Unknown)
{
    
}

I3CLSimLightSourceParameterization::I3CLSimLightSourceParameterization
(I3CLSimLightSourceToStepConverterPtr converter_,
 const I3Particle &forParticleType_,
 double fromEnergy_, double toEnergy_,
 bool needsLength_
 )
:
converter(converter_),
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
forPdgEncoding(forParticleType_.GetPdgEncoding()),
#else
forParticleType(forParticleType_.GetType()),
#endif
fromEnergy(fromEnergy_),
toEnergy(toEnergy_),
needsLength(needsLength_),
catchAll(false),
flasherMode(false),
forFlasherPulseType(I3CLSimFlasherPulse::Unknown)
{
    
}

I3CLSimLightSourceParameterization::I3CLSimLightSourceParameterization
(I3CLSimLightSourceToStepConverterPtr converter_,
 const I3CLSimLightSourceParameterization::AllParticles_t &dummy,
 double fromEnergy_, double toEnergy_,
 bool needsLength_
 )
:
converter(converter_),
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
forPdgEncoding(0),
#else
forParticleType(I3Particle::unknown),
#endif
fromEnergy(fromEnergy_),
toEnergy(toEnergy_),
needsLength(needsLength_),
catchAll(true),
flasherMode(false),
forFlasherPulseType(I3CLSimFlasherPulse::Unknown)
{
    
}

I3CLSimLightSourceParameterization::I3CLSimLightSourceParameterization
(I3CLSimLightSourceToStepConverterPtr converter_,
 I3CLSimFlasherPulse::FlasherPulseType forFlasherPulseType_
)
:
converter(converter_),
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
forPdgEncoding(0),
#else
forParticleType(I3Particle::unknown),
#endif
fromEnergy(NAN),
toEnergy(NAN),
needsLength(false),
catchAll(false),
flasherMode(true),
forFlasherPulseType(forFlasherPulseType_)
{
    
}

bool I3CLSimLightSourceParameterization::IsValidForParticle(const I3Particle &particle) const
{
    if (flasherMode) return false; // flasher params. cannot be valid for particles
    
    if (std::isnan(particle.GetEnergy())) {
        log_warn("I3CLSimLightSourceParameterization::IsValid() called with particle with NaN energy. Parameterization is NOT valid.");
        return false;
    }

    return IsValid(particle.GetType(), particle.GetEnergy(), particle.GetLength());
}

bool I3CLSimLightSourceParameterization::IsValid(I3Particle::ParticleType type, double energy, double length) const
{
    if (flasherMode) return false; // flasher params. cannot be valid for particles

    if (std::isnan(energy)) return false;
    if (!catchAll) {
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        const int32_t encoding = I3Particle(I3Particle::Null,type).GetPdgEncoding();
        if (encoding==0) return false;
        if (encoding != forPdgEncoding) return false;
#else
        if (type != forParticleType) return false;
#endif
    }
    if ((energy < fromEnergy) || (energy > toEnergy)) return false;
    if ((needsLength) && (std::isnan(length))) return false;
    
    return true;
}

#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
bool I3CLSimLightSourceParameterization::IsValidForPdgEncoding(int32_t encoding, double energy, double length) const
{
    if (flasherMode) return false; // flasher params. cannot be valid for particles

    if (std::isnan(energy)) return false;
    if (!catchAll) {
        if (encoding != forPdgEncoding) return false;
    }
    if ((energy < fromEnergy) || (energy > toEnergy)) return false;
    if ((needsLength) && (std::isnan(length))) return false;
    
    return true;
}
#endif

bool I3CLSimLightSourceParameterization::IsValidForLightSource(const I3CLSimLightSource &lightSource) const
{
    if (lightSource.GetType() == I3CLSimLightSource::Particle) {
        return IsValidForParticle(lightSource.GetParticle());
    } else if (lightSource.GetType() == I3CLSimLightSource::Flasher) {
        if (!flasherMode) return false;

        // return true if the flasher pulse has the same type as the configured type
        return (lightSource.GetFlasherPulse().GetType() == forFlasherPulseType);
    } else {
        log_error("Parameterization got light source with invalid or unknown type.");
        return false;
    }
}
