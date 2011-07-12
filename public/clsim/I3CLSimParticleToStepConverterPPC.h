#ifndef I3CLSimParticleToStepConverterPPC_H_INCLUDED
#define I3CLSimParticleToStepConverterPPC_H_INCLUDED

#include "clsim/I3CLSimParticleToStepConverter.h"
#include "dataclasses/physics/I3Particle.h"

#include "phys-services/I3RandomService.h"

#include "clsim/I3CLSimQueue.h"

#include <map>
#include <string>
#include <deque>

#include <boost/variant.hpp>

/**
 * @brief A particle-to-step converter for cascades
 * using pre-defined parameterizations.
 */
struct I3CLSimParticleToStepConverterPPC : public I3CLSimParticleToStepConverter
{
public:
    static const uint32_t default_photonsPerStep;
    static const uint32_t default_highPhotonsPerStep;
    static const double default_useHighPhotonsPerStepStartingFromNumPhotons;

    I3CLSimParticleToStepConverterPPC(I3RandomServicePtr randomService,
                                      uint32_t photonsPerStep=default_photonsPerStep,
                                      uint32_t highPhotonsPerStep=default_highPhotonsPerStep,
                                      double useHighPhotonsPerStepStartingFromNumPhotons=default_useHighPhotonsPerStepStartingFromNumPhotons);
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
    ///////////////
    // definitions used in the internal queue
    
    struct CascadeStepData_t {
        I3Particle particle;
        uint32_t particleIdentifier;
        uint64_t photonsPerStep;
        uint64_t numSteps;
        uint64_t numPhotonsInLastStep;
        
        double pa,pb;
    };
    struct MuonStepData_t {
        I3Particle particle;
        uint32_t particleIdentifier;
        uint64_t photonsPerStep;
        uint64_t numSteps;
        uint64_t numPhotonsInLastStep;

        bool stepIsCascadeLike;
        double length;
    };
    struct BarrierData_t {
    };
    typedef boost::variant<CascadeStepData_t, MuonStepData_t, BarrierData_t> StepData_t;
    
    std::deque<StepData_t> stepGenerationQueue_;
    
    I3CLSimStepSeriesConstPtr MakeSteps(bool &barrierWasReset);
    class MakeSteps_visitor : public boost::static_visitor<std::pair<I3CLSimStepSeriesConstPtr, bool> >
    {
    public:
        MakeSteps_visitor(I3RandomService &randomService, uint64_t maxNumStepsPerStepSeries);
        template <typename T>
        std::pair<I3CLSimStepSeriesConstPtr, bool> operator()(T &data) const;
        
    private:
        void FillStep(I3CLSimParticleToStepConverterPPC::CascadeStepData_t &data, I3CLSimStep &newStep, uint64_t photonsPerStep) const;
        void FillStep(I3CLSimParticleToStepConverterPPC::MuonStepData_t &data, I3CLSimStep &newStep, uint64_t photonsPerStep) const;
        
        I3RandomService &randomService_;
        uint64_t maxNumStepsPerStepSeries_;
    };
    //////////////////
    
    I3RandomServicePtr randomService_;
    
    bool initialized_;
    bool barrier_is_enqueued_;
    uint64_t bunchSizeGranularity_;
    uint64_t maxBunchSize_;
    uint32_t photonsPerStep_;
    uint32_t highPhotonsPerStep_;
    double useHighPhotonsPerStepStartingFromNumPhotons_;
    
    I3CLSimWlenDependentValueConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    
    std::vector<double> meanPhotonsPerMeterInLayer_;
};

I3_POINTER_TYPEDEFS(I3CLSimParticleToStepConverterPPC);

#endif //I3CLSimParticleToStepConverterPPC_H_INCLUDED
