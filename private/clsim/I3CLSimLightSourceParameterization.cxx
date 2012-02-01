#include <icetray/serialization.h>
#include <clsim/I3CLSimLightSourceParameterization.h>

#include <limits>

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
catchFlashers(false)
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
forPdgEncoding(I3Particle::ConvertToPdgEncoding(forParticleType_)),
#else
forParticleType(forParticleType_),
#endif
fromEnergy(fromEnergy_),
toEnergy(toEnergy_),
needsLength(needsLength_),
catchAll(false),
catchFlashers(false)
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
catchFlashers(false)
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
catchFlashers(false)
{
    
}

I3CLSimLightSourceParameterization::I3CLSimLightSourceParameterization
(I3CLSimLightSourceToStepConverterPtr converter_,
 const I3CLSimLightSourceParameterization::AllFlashers_t &dummy
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
catchFlashers(true)
{
    
}

bool I3CLSimLightSourceParameterization::IsValidForParticle(const I3Particle &particle) const
{
    if (catchFlashers) return false; // flasher params. cannot be valid for particles
    
    if (isnan(particle.GetEnergy())) {
        log_warn("I3CLSimLightSourceParameterization::IsValid() called with particle with NaN energy. Parameterization is NOT valid.");
        return false;
    }

    return IsValid(particle.GetType(), particle.GetEnergy(), particle.GetLength());
}

bool I3CLSimLightSourceParameterization::IsValid(I3Particle::ParticleType type, double energy, double length) const
{
    if (catchFlashers) return false; // flasher params. cannot be valid for particles

    if (isnan(energy)) return false;
    if (!catchAll) {
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        const int32_t encoding = I3Particle::ConvertToPdgEncoding(type);
        if (encoding==0) return false;
        if (encoding != forPdgEncoding) return false;
#else
        if (type != forParticleType) return false;
#endif
    }
    if ((energy < fromEnergy) || (energy > toEnergy)) return false;
    if ((needsLength) && (isnan(length))) return false;
    
    return true;
}

#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
bool I3CLSimLightSourceParameterization::IsValidForPdgEncoding(int32_t encoding, double energy, double length) const
{
    if (catchFlashers) return false; // flasher params. cannot be valid for particles

    if (isnan(energy)) return false;
    if (!catchAll) {
        if (encoding != forPdgEncoding) return false;
    }
    if ((energy < fromEnergy) || (energy > toEnergy)) return false;
    if ((needsLength) && (isnan(length))) return false;
    
    return true;
}
#endif

bool I3CLSimLightSourceParameterization::IsValidForLightSource(const I3CLSimLightSource &lightSource) const
{
    if (lightSource.GetType() == I3CLSimLightSource::Particle) {
        return IsValidForParticle(lightSource.GetParticle());
    } else if (lightSource.GetType() == I3CLSimLightSource::Flasher) {
        if (catchFlashers) return true;
        return false;
    } else {
        log_error("Parameterization got light source with invalid or unknown type.");
        return false;
    }
}
