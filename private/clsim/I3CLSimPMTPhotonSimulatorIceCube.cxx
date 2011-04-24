#include "clsim/I3CLSimPMTPhotonSimulatorIceCube.h"

#include "hit-maker/AfterPulseGenerator.h"
#include "hit-maker/LatePulseGenerator.h"

#include <boost/foreach.hpp>

const double I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_jitter = 2.*I3Units::ns;
const double I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_pre_pulse_probability=0.007;
const double I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_late_pulse_probability=0.035;
const double I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_after_pulse_probability=0.0593;


I3CLSimPMTPhotonSimulatorIceCube::I3CLSimPMTPhotonSimulatorIceCube(I3RandomServicePtr randomService,
                                                                   double jitter,
                                                                   double pre_pulse_probability,
                                                                   double late_pulse_probability,
                                                                   double after_pulse_probability)
:
randomService_(randomService),
jitter_(jitter),
pre_pulse_probability_(pre_pulse_probability),
late_pulse_probability_(late_pulse_probability),
after_pulse_probability_(after_pulse_probability)
{ 
    afterPulseGenerator_ = AfterPulseGeneratorPtr(new AfterPulseGenerator(randomService));;
    latePulseGenerator_ = LatePulseGeneratorPtr(new LatePulseGenerator(randomService));;
}

I3CLSimPMTPhotonSimulatorIceCube::~I3CLSimPMTPhotonSimulatorIceCube() 
{ 
    
}

void I3CLSimPMTPhotonSimulatorIceCube::SetCalibration(I3CalibrationConstPtr calibration)
{
    calibration_=calibration;
}

void I3CLSimPMTPhotonSimulatorIceCube::SetDetectorStatus(I3DetectorStatusConstPtr status)
{
    status_=status;
}

namespace {
    inline double jitter(I3RandomServicePtr randomService)
    {
        for(;;) 
        {
            const double time = randomService->Gaus(0.,2.*I3Units::ns);
            
            if ((time >= 5.*I3Units::ns) && (time <= 25.*I3Units::ns))
                return time;
        }
    }
    
    inline double pts(double v)
    {
        return -31.8*I3Units::ns*sqrt(1345*I3Units::V/v);
    }

}

void I3CLSimPMTPhotonSimulatorIceCube::
ApplyAfterPulseLatePulseAndJitterSim(const OMKey &key,
                                     const I3MCHit &input_hit,
                                     std::vector<I3MCHit> &output_vector) const
{
    if (!calibration_) log_fatal("no calibration has been set");
    if (!status_) log_fatal("no detector status has been set");

    if (isnan(input_hit.GetTime())) { // do nothing for 
        log_debug("Hit with NaN time. Adding it unchanged.");
        output_vector.push_back(input_hit);
        return;
    }

    // retrieve the HV for this DOM
    std::map<OMKey, I3DOMStatus>::const_iterator om_stat = status_->domStatus.find(key);
    if (om_stat==status_->domStatus.end())
        log_fatal("No DOMStatus entry found for OMKey(%i,%u)",
                  key.GetString(), key.GetOM());
    const double pmtHV = om_stat->second.pmtHV;
    
    
    // add the input hit to the output vector
    output_vector.push_back(input_hit);
    I3MCHit &output_hit = output_vector.back();
    
    
    // taken from PPC by D. Chirkin
    const double rnd = randomService_->Uniform();

    double hit_time = output_hit.GetTime();
    
    if(rnd < pre_pulse_probability_)
    {
        output_hit.SetHitSource(I3MCHit::PRE_PULSE);
        hit_time += pts(pmtHV) + jitter(randomService_);
    }
    else if(rnd < pre_pulse_probability_+late_pulse_probability_)
    {
        output_hit.SetHitSource(I3MCHit::ELASTIC_LATE_PULSE);
        hit_time = latePulseGenerator_->GenerateSingleLatePulse(hit_time, pmtHV);
    }
    else
    {
        output_hit.SetHitSource(I3MCHit::SPE);
        hit_time += jitter(randomService_);
    }
    
    output_hit.SetTime(hit_time);

    
    if (randomService_->Uniform() < after_pulse_probability_)
    {
        std::vector<I3MCHit> afterPulses;
        afterPulseGenerator_->GenerateAfterPulses(afterPulses, output_hit);
        
        BOOST_FOREACH(const I3MCHit &afterPulse, afterPulses)
        {
            output_vector.push_back(afterPulse);
        }
    }
    // end - PPC code
    
    
}
