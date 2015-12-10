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
 * @file I3CLSimLightSourceToStepConverterPPC.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMLIGHTSOURCETOSTEPCONVERTERPPC_H_INCLUDED
#define I3CLSIMLIGHTSOURCETOSTEPCONVERTERPPC_H_INCLUDED

#include "clsim/I3CLSimLightSourceToStepConverter.h"
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
namespace I3CLSimLightSourceToStepConverterUtils {
    class GenerateStepPreCalculator;
}

/**
 * @brief A particle-to-step converter for cascades
 * using pre-defined parameterizations.
 */
struct I3CLSimLightSourceToStepConverterPPC : public I3CLSimLightSourceToStepConverter
{
public:
    static const uint32_t default_photonsPerStep;
    static const uint32_t default_highPhotonsPerStep;
    static const double default_useHighPhotonsPerStepStartingFromNumPhotons;

    I3CLSimLightSourceToStepConverterPPC(uint32_t photonsPerStep=default_photonsPerStep,
                                      uint32_t highPhotonsPerStep=default_highPhotonsPerStep,
                                      double useHighPhotonsPerStepStartingFromNumPhotons=default_useHighPhotonsPerStepStartingFromNumPhotons);
    virtual ~I3CLSimLightSourceToStepConverterPPC();

    void SetUseCascadeExtension(bool v) { useCascadeExtension_ = v; };

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
        void FillStep(I3CLSimLightSourceToStepConverterPPC::CascadeStepData_t &data,
                      I3CLSimStep &newStep,
                      uint64_t photonsPerStep,
                      double particleDir_x, double particleDir_y, double particleDir_z
                     ) const;
        void FillStep(I3CLSimLightSourceToStepConverterPPC::MuonStepData_t &data,
                      I3CLSimStep &newStep,
                      uint64_t photonsPerStep,
                      double particleDir_x, double particleDir_y, double particleDir_z
                     ) const;
        
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
    bool useCascadeExtension_;
    
    I3CLSimFunctionConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    
    std::vector<double> meanPhotonsPerMeterInLayer_;
    
    boost::shared_ptr<GenerateStepPreCalculator> preCalc_;
    
    
    
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
        double one_over_angularDist_a_;
        double angularDist_b_;
        double angularDist_I_;
        
        std::size_t numberOfValues_;
        std::size_t index_;

        typedef std::vector<std::pair<std::pair<double, double>, double> > queueVector_t;
        boost::shared_ptr<queueVector_t> currentVector_;
        
        I3CLSimQueue<boost::shared_ptr<queueVector_t> > queueFromFeederThreads_;
        std::vector<boost::shared_ptr<boost::thread> > feederThreads_;
        
        void FeederThread(unsigned int threadId, uint64_t initialRngState, uint32_t rngA);
        void RegenerateValues();
    };

    
    static void GenerateStep(I3CLSimStep &newStep,
                             const I3Particle &p,
                             double particleDir_x, double particleDir_y, double particleDir_z,
                             uint32_t identifier,
                             uint32_t photonsPerStep,
                             const double &longitudinalPos,
                             GenerateStepPreCalculator &preCalc);

    static void GenerateStepForMuon(I3CLSimStep &newStep,
                                    const I3Particle &p,
                                    double particleDir_x, double particleDir_y, double particleDir_z,
                                    uint32_t identifier,
                                    uint32_t photonsPerStep,
                                    double length);

};

// forward-declare template specialization
template <>
std::pair<I3CLSimStepSeriesConstPtr, bool>
I3CLSimLightSourceToStepConverterPPC::MakeSteps_visitor::operator()
(I3CLSimLightSourceToStepConverterPPC::BarrierData_t &data) const;

I3_POINTER_TYPEDEFS(I3CLSimLightSourceToStepConverterPPC);

#endif //I3CLSIMLIGHTSOURCETOSTEPCONVERTERPPC_H_INCLUDED
