#include <icetray/serialization.h>
#include <clsim/I3CLSimLightSourceToStepConverter.h>

I3CLSimLightSourceToStepConverter::I3CLSimLightSourceToStepConverter() {;}
I3CLSimLightSourceToStepConverter::~I3CLSimLightSourceToStepConverter() {;}

void I3CLSimLightSourceToStepConverter::SetParticleParameterizationSeries
(const I3CLSimLightSourceParameterizationSeries &parameterizationSeries_)
{
    if (IsInitialized())
        throw I3CLSimLightSourceToStepConverter_exception("SetParticleParameterizationSeries() called after Initialize().");
        
    parameterizationSeries = parameterizationSeries_;
}

const I3CLSimLightSourceParameterizationSeries &
I3CLSimLightSourceToStepConverter::GetParticleParameterizationSeries() const
{ 
    return parameterizationSeries;
}

I3CLSimStepSeriesConstPtr I3CLSimLightSourceToStepConverter::GetConversionResult(double timeout)
{
    bool dummy;
    return GetConversionResultWithBarrierInfo(dummy, timeout);
}
