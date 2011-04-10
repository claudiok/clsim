#include <icetray/serialization.h>
#include <clsim/I3CLSimParticleParameterization.h>

#include <limits>

I3CLSimParticleParameterization::I3CLSimParticleParameterization()
:
forParticleType(I3Particle::unknown),
fromEnergy(0.),
toEnergy(std::numeric_limits<double>::infinity())
{ 
    
}

I3CLSimParticleParameterization::~I3CLSimParticleParameterization() 
{ 

}

I3CLSimParticleParameterization::I3CLSimParticleParameterization
(I3CLSimParticleToStepConverterPtr converter_,
 I3Particle::ParticleType forParticleType_,
 double fromEnergy_, double toEnergy_,
 bool needsLength_
 )
:
converter(converter_),
forParticleType(forParticleType_),
fromEnergy(fromEnergy_),
toEnergy(toEnergy_),
needsLength(needsLength_)
{
    
}

bool I3CLSimParticleParameterization::IsValidForParticle(const I3Particle &particle) const
{
    if (isnan(particle.GetEnergy())) {
        log_warn("I3CLSimParticleParameterization::IsValid() called with particle with NaN energy. Parameterization is NOT valid.");
        return false;
    }
    
    return IsValid(particle.GetType(), particle.GetEnergy(), particle.GetLength());
}

bool I3CLSimParticleParameterization::IsValid(I3Particle::ParticleType type, double energy, double length) const
{
    if (isnan(energy)) return false;
    if (type != forParticleType) return false;
    if ((energy < fromEnergy) || (energy > toEnergy)) return false;
    if ((needsLength) && (isnan(length))) return false;
    
    return true;
}
