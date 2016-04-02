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
 * @file I3CLSimLightSourceParameterization.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMLIGHTSOURCEPARAMETERIZATION_H_INCLUDED
#define I3CLSIMLIGHTSOURCEPARAMETERIZATION_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include "dataclasses/physics/I3Particle.h"
#include "clsim/I3CLSimLightSource.h"

#include <vector>

// forward declarations
struct I3CLSimLightSourceToStepConverter;
I3_POINTER_TYPEDEFS(I3CLSimLightSourceToStepConverter);

/**
 * @brief Defines a particle parameterization for fast simulation,
 * bypassing Geant4
 */
struct I3CLSimLightSourceParameterization 
{
public:
    I3CLSimLightSourceParameterization();
    ~I3CLSimLightSourceParameterization();

    // for I3Particle
    I3CLSimLightSourceParameterization(I3CLSimLightSourceToStepConverterPtr converter_,
                                       I3Particle::ParticleType forParticleType_,
                                       double fromEnergy_, double toEnergy_,
                                       bool needsLength_=false);

    I3CLSimLightSourceParameterization(I3CLSimLightSourceToStepConverterPtr converter_,
                                       const I3Particle &forParticleType_,
                                       double fromEnergy_, double toEnergy_,
                                       bool needsLength_=false);

    struct AllParticles_t {};
    static const AllParticles_t AllParticles;
    I3CLSimLightSourceParameterization(I3CLSimLightSourceToStepConverterPtr converter_,
                                       const AllParticles_t &,
                                       double fromEnergy_, double toEnergy_,
                                       bool needsLength_=false);

    // for I3CLSimFlasherPulse
    I3CLSimLightSourceParameterization(I3CLSimLightSourceToStepConverterPtr converter_,
                                       I3CLSimFlasherPulse::FlasherPulseType forFlasherPulseType_);

    
    I3CLSimLightSourceToStepConverterPtr converter;
#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
    I3Particle::ParticleType forParticleType;
#else
    int32_t forPdgEncoding;
#endif
    double fromEnergy, toEnergy;
    bool needsLength;
    bool catchAll;

    bool flasherMode;
    I3CLSimFlasherPulse::FlasherPulseType forFlasherPulseType;
    
    bool IsValidForParticle(const I3Particle &particle) const;
    bool IsValid(I3Particle::ParticleType type, double energy, double length=NAN) const;
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
    bool IsValidForPdgEncoding(int32_t encoding, double energy, double length=NAN) const;
#endif
    
    bool IsValidForLightSource(const I3CLSimLightSource &lightSource) const;
    
private:
};

inline bool operator==(const I3CLSimLightSourceParameterization &a, const I3CLSimLightSourceParameterization &b)
{
    if (a.converter != b.converter) return false;
    if (a.flasherMode != b.flasherMode) return false;

    if (a.flasherMode) {
        if (a.fromEnergy != b.fromEnergy) return false;
        if (a.toEnergy != b.toEnergy) return false;
#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        if (a.forParticleType != b.forParticleType) return false;
#else
        if (a.forPdgEncoding != b.forPdgEncoding) return false;
#endif
        if (a.needsLength != b.needsLength) return false;
        if (a.catchAll != b.catchAll) return false;
    } else {
        if (a.forFlasherPulseType != b.forFlasherPulseType) return false;
    }
    
    return true;
}

typedef std::vector<I3CLSimLightSourceParameterization> I3CLSimLightSourceParameterizationSeries;

I3_POINTER_TYPEDEFS(I3CLSimLightSourceParameterization);
I3_POINTER_TYPEDEFS(I3CLSimLightSourceParameterizationSeries);

#endif //I3CLSIMLIGHTSOURCEPARAMETERIZATION_H_INCLUDED
