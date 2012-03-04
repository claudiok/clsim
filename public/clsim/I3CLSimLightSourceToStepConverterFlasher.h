#ifndef I3CLSIMLIGHTSOURCETOSTEPCONVERTERFLASHER_H_INCLUDED
#define I3CLSIMLIGHTSOURCETOSTEPCONVERTERFLASHER_H_INCLUDED

#include "clsim/I3CLSimLightSourceToStepConverter.h"
#include "clsim/I3CLSimLightSource.h"
#include "clsim/function/I3CLSimFunction.h"
#include "clsim/I3CLSimSpectrumTable.h"

#include <string>
#include <vector>
#include <deque>



/**
 * @brief A particle-to-step converter for IceCube LED flashers.
 */
struct I3CLSimLightSourceToStepConverterFlasher : public I3CLSimLightSourceToStepConverter
{
public:
    static const uint32_t default_photonsPerStep;

    I3CLSimLightSourceToStepConverterFlasher(I3CLSimFunctionConstPtr flasherSpectrumNoBias,
                                             I3CLSimSpectrumTablePtr spectrumTable,
                                             uint32_t photonsPerStep=default_photonsPerStep);

    virtual ~I3CLSimLightSourceToStepConverterFlasher();

    // inherited:
    
    virtual void SetBunchSizeGranularity(uint64_t num);

    virtual void SetMaxBunchSize(uint64_t num);

    virtual void SetRandomService(I3RandomServicePtr random);

    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias);

    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties);
    
    virtual void Initialize();

    virtual bool IsInitialized() const;
    
    virtual void EnqueueLightSource(const I3CLSimLightSource &lightSource, uint32_t identifier);
    
    virtual void EnqueueBarrier();
    
    virtual bool BarrierActive() const;
    
    virtual bool MoreStepsAvailable() const;

    virtual I3CLSimStepSeriesConstPtr GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout=NAN);
    
private:
    // this function performs the actual conversion
    I3CLSimStepSeriesConstPtr MakeSteps(bool &barrierWasReset);

    void FillStep(I3CLSimStep &step,
                  uint32_t numberOfPhotons,
                  const I3CLSimFlasherPulse &flasherPulse,
                  uint32_t identifier);

    ///////////////
    // definitions used in the internal queue
    
    struct LightSourceData_t {
        bool isBarrier;
        I3CLSimFlasherPulse flasherPulse;
        uint32_t identifier;
        
        uint64_t numPhotonsWithBias;
    };
    std::deque<LightSourceData_t> inputQueue_;

    
    
    I3RandomServicePtr randomService_;
    
    bool initialized_;
    bool barrier_is_enqueued_;
    uint64_t bunchSizeGranularity_;
    uint64_t maxBunchSize_;
    uint32_t photonsPerStep_;
    uint8_t spectrumSourceTypeIndex_;
    
    I3CLSimFunctionConstPtr flasherSpectrumNoBias_;
    I3CLSimFunctionConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    
    double photonNumberCorrectionFactorForBias_;
    

};

I3_POINTER_TYPEDEFS(I3CLSimLightSourceToStepConverterFlasher);

#endif //I3CLSIMLIGHTSOURCETOSTEPCONVERTERFLASHER_H_INCLUDED
