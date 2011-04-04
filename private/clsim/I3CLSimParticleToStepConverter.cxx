#include <icetray/serialization.h>
#include <clsim/I3CLSimParticleToStepConverter.h>

I3CLSimParticleToStepConverter::I3CLSimParticleToStepConverter() {;}
I3CLSimParticleToStepConverter::~I3CLSimParticleToStepConverter() {;}

void I3CLSimParticleToStepConverter::SetParticleParameterizationSeries
(const I3CLSimParticleParameterizationSeries &parameterizationSeries_)
{
    if (IsInitialized())
        throw I3CLSimParticleToStepConverter_exception("SetParticleParameterizationSeries() called after Initialize().");
        
    parameterizationSeries = parameterizationSeries_;
}

const I3CLSimParticleParameterizationSeries &
I3CLSimParticleToStepConverter::GetParticleParameterizationSeries() const
{ 
    return parameterizationSeries;
}

I3CLSimStepSeriesConstPtr I3CLSimParticleToStepConverter::GetConversionResult(double timeout)
{
    bool dummy;
    return GetConversionResultWithBarrierInfo(dummy, timeout);
}
