#ifndef I3CLSimParticleToStepConverterPPC_H_INCLUDED
#define I3CLSimParticleToStepConverterPPC_H_INCLUDED

#include "clsim/I3CLSimParticleToStepConverter.h"
#include "dataclasses/physics/I3Particle.h"

#include "phys-services/I3RandomService.h"

#include "clsim/I3CLSimQueue.h"

#include <map>
#include <string>

/**
 * @brief A particle-to-step converter for cascades
 * using pre-defined parameterizations.
 */
struct I3CLSimParticleToStepConverterPPC : public I3CLSimParticleToStepConverter
{
public:
    static const uint32_t default_photonsPerStep;
    
    I3CLSimParticleToStepConverterPPC(I3RandomServicePtr randomService,
                                                          uint32_t photonsPerStep=default_photonsPerStep);
    virtual ~I3CLSimParticleToStepConverterPPC();

    // inherited:
    
    virtual void SetBunchSizeGranularity(uint64_t num);

    virtual void SetMaxBunchSize(uint64_t num);

    virtual void SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias);

    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties);
    
    virtual void Initialize();

    virtual bool IsInitialized() const;
    
    virtual void EnqueueParticle(const I3Particle &particle, uint32_t identifier);
    
    virtual void EnqueueBarrier();
    
    virtual bool BarrierActive() const;
    
    virtual bool MoreStepsAvailable() const;

    virtual I3CLSimStepSeriesConstPtr GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout=NAN);
    
private:
    I3CLSimStepSeriesPtr currentStepSeries_;
    
    I3RandomServicePtr randomService_;
    
    bool initialized_;
    bool barrier_is_enqueued_;
    uint64_t bunchSizeGranularity_;
    uint64_t maxBunchSize_;
    uint32_t photonsPerStep_;

    I3CLSimWlenDependentValueConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    
    std::vector<double> meanPhotonsPerMeterInLayer_;
};

I3_POINTER_TYPEDEFS(I3CLSimParticleToStepConverterPPC);

#endif //I3CLSimParticleToStepConverterPPC_H_INCLUDED
