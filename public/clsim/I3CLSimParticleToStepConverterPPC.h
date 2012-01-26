#ifndef I3CLSimParticleToStepConverterPPC_H_INCLUDED
#define I3CLSimParticleToStepConverterPPC_H_INCLUDED

#include "clsim/I3CLSimParticleToStepConverter.h"
#include "dataclasses/physics/I3Particle.h"

#include "clsim/I3CLSimQueue.h"

#include <map>
#include <string>
#include <vector>
#include <deque>

#include <boost/variant.hpp>

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>


// forward decl
namespace I3CLSimParticleToStepConverterUtils {
    class GenerateStepPreCalculator;
}

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

    I3CLSimParticleToStepConverterPPC(uint32_t photonsPerStep=default_photonsPerStep,
                                      uint32_t highPhotonsPerStep=default_highPhotonsPerStep,
                                      double useHighPhotonsPerStepStartingFromNumPhotons=default_useHighPhotonsPerStepStartingFromNumPhotons);
    virtual ~I3CLSimParticleToStepConverterPPC();

    // inherited:
    
    virtual void SetBunchSizeGranularity(uint64_t num);

    virtual void SetMaxBunchSize(uint64_t num);

    virtual void SetRandomService(I3RandomServicePtr random);

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
    
    // forward declaration
    class GenerateStepPreCalculator;

    
    I3CLSimStepSeriesConstPtr MakeSteps(bool &barrierWasReset);
    class MakeSteps_visitor : public boost::static_visitor<std::pair<I3CLSimStepSeriesConstPtr, bool> >
    {
    public:
        MakeSteps_visitor(uint64_t &rngState, uint32_t rngA,
                          uint64_t maxNumStepsPerStepSeries,
                          GenerateStepPreCalculator &preCalc);
        template <typename T>
        std::pair<I3CLSimStepSeriesConstPtr, bool> operator()(T &data) const;
        
    private:
        void FillStep(I3CLSimParticleToStepConverterPPC::CascadeStepData_t &data, I3CLSimStep &newStep, uint64_t photonsPerStep) const;
        void FillStep(I3CLSimParticleToStepConverterPPC::MuonStepData_t &data, I3CLSimStep &newStep, uint64_t photonsPerStep) const;
        
        uint64_t &rngState_;
        uint32_t rngA_;
        //I3RandomService &randomService_;
        uint64_t maxNumStepsPerStepSeries_;
        GenerateStepPreCalculator &preCalc_;
    };
    //////////////////
    
    I3RandomServicePtr randomService_;
    uint64_t rngState_;
    uint32_t rngA_;
    
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
    
    shared_ptr<GenerateStepPreCalculator> preCalc_;
    
    
    
    ////////////////////
    // HELPERS
    ////////////////////
    
    class GenerateStepPreCalculator
    {
    public:
        GenerateStepPreCalculator(I3RandomServicePtr randomService,
                                  double angularDist_a=0.39,
                                  double angularDist_b=2.61,
                                  std::size_t numberOfValues=102400);
        ~GenerateStepPreCalculator();
        
        inline void GetAngularCosSinValue(double &angular_cos, double &angular_sin, double &random_value)
        {
            if (index_ >= numberOfValues_) RegenerateValues();
            
            const std::pair<std::pair<double, double>, double> &currentPair = (*currentVector_)[index_];
            
            angular_sin = currentPair.first.first;
            angular_cos = currentPair.first.second;
            random_value = currentPair.second;
            
            ++index_;
        }
        
    private:
        double angularDist_a_;
        double one_over_angularDist_a_;
        double angularDist_b_;
        double angularDist_I_;
        
        std::size_t numberOfValues_;
        std::size_t index_;

        typedef std::vector<std::pair<std::pair<double, double>, double> > queueVector_t;
        shared_ptr<queueVector_t> currentVector_;
        
        I3CLSimQueue<shared_ptr<queueVector_t> > queueFromFeederThreads_;
        std::vector<shared_ptr<boost::thread> > feederThreads_;
        
        void FeederThread(unsigned int threadId, uint64_t initialRngState, uint32_t rngA);
        void RegenerateValues();
    };

    
    static void GenerateStep(I3CLSimStep &newStep,
                             const I3Particle &p,
                             uint32_t identifier,
                             uint32_t photonsPerStep,
                             const double &longitudinalPos,
                             GenerateStepPreCalculator &preCalc);

    static void GenerateStepForMuon(I3CLSimStep &newStep,
                                    const I3Particle &p,
                                    uint32_t identifier,
                                    uint32_t photonsPerStep,
                                    double length);

};

// forward-declare template specialization
template <>
std::pair<I3CLSimStepSeriesConstPtr, bool>
I3CLSimParticleToStepConverterPPC::MakeSteps_visitor::operator()
(I3CLSimParticleToStepConverterPPC::BarrierData_t &data) const;

I3_POINTER_TYPEDEFS(I3CLSimParticleToStepConverterPPC);

#endif //I3CLSimParticleToStepConverterPPC_H_INCLUDED
