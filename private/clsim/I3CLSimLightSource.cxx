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

I3CLSimLightSource::I3CLSimLightSource(const I3FlasherInfo &flasher)
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


const I3FlasherInfo &I3CLSimLightSource::GetFlasherInfo() const
{
    if (lightSourceType_ != Flasher)
        throw std::runtime_error("this light source is not an I3FlasherInfo");
    
    return flasher_;
}

