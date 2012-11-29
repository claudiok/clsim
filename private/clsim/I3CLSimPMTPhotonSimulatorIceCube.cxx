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
 * @file I3CLSimPMTPhotonSimulatorIceCube.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "clsim/I3CLSimPMTPhotonSimulatorIceCube.h"

#include "hit-maker/AfterPulseGenerator.h"
#include "hit-maker/LatePulseGenerator.h"

#include <boost/foreach.hpp>

const double I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_jitter = 2.*I3Units::ns;
const double I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_pre_pulse_probability=0.007;
const double I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_late_pulse_probability=0.035;
const double I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_after_pulse_probability=0.0593;
const bool I3CLSimPMTPhotonSimulatorIceCube::DEFAULT_ppc_jitter_mode=false;


I3CLSimPMTPhotonSimulatorIceCube::I3CLSimPMTPhotonSimulatorIceCube(double jitter,
                                                                   double pre_pulse_probability,
                                                                   double late_pulse_probability,
                                                                   double after_pulse_probability,
                                                                   bool ppc_jitter_mode)
:
//randomService_(randomService),
jitter_(jitter),
pre_pulse_probability_(pre_pulse_probability),
late_pulse_probability_(late_pulse_probability),
after_pulse_probability_(after_pulse_probability),
ppc_jitter_mode_(ppc_jitter_mode)
{ 
    log_trace("pre_pulse_probability=%g, late_pulse_probability=%g, after_pulse_probability=%g",
              pre_pulse_probability_, late_pulse_probability_, after_pulse_probability_);

    if ((!ppc_jitter_mode) && (jitter != DEFAULT_jitter))
        log_fatal("You need to disable ppc_jitter_mode in order to use non-default jitter values!");

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

void I3CLSimPMTPhotonSimulatorIceCube::SetRandomService(I3RandomServicePtr random)
{
    randomService_=random;
}

double I3CLSimPMTPhotonSimulatorIceCube::GenerateJitter() const
{
    if (ppc_jitter_mode_)
    {
        // PPC samples from a cut-off gaussian distribution
        for(;;) 
        {
            const double time = randomService_->Gaus(0.,jitter_);
            
            if ((time >= -5.*I3Units::ns) && (time <= 25.*I3Units::ns))
                return time;
        }
    }
    else
    {
        // this is taken from hit-maker:

        // Quoting hit-maker: 
        // """
        // Fisher-Tippett variables for well behaved time
        // distribtutions.  This replaces the truncated gaussian
        // These values were obtained by fits to Bricktop running
        // at 1345 during DFL studies by C. Wendt.  The fits were
        // performed by R.Porrata.
        // """

        const double FT_BETA(1.9184 * I3Units::nanosecond);                
        const double FT_MU(0.15304 * I3Units::nanosecond);                 
                                                                   
        const double FT_LN_TIME_UPPER_BOUND(0.999998); // Upper bound also needed for FT. tmax ~ 25 ns.
        const double FT_LN_TIME_LOWER_BOUND(1e-7);     // Lower bound ~ 5 ns.
        //const double EARLIEST_TIME(FT_MU - FT_BETA * std::log(-std::log(FT_LN_TIME_LOWER_BOUND)));

        return FT_MU - FT_BETA * log(-log(randomService_->Uniform(FT_LN_TIME_LOWER_BOUND,FT_LN_TIME_UPPER_BOUND)));
    }
}

namespace {
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
    if (!randomService_) log_fatal("no random number generator service has been set");

    if (isnan(input_hit.GetTime())) { // do nothing for 
        log_debug("Hit with NaN time. Adding it unchanged.");
        output_vector.push_back(input_hit);
        return;
    }
    
    bool jitterOnly = ((pre_pulse_probability_<=0.) &&
                       (late_pulse_probability_<=0.) &&
                       (after_pulse_probability_<=0.));

    if (!jitterOnly)
    {
        if (!afterPulseGenerator_)
            afterPulseGenerator_ = AfterPulseGeneratorPtr(new AfterPulseGenerator(randomService_));
        if (!latePulseGenerator_)
            latePulseGenerator_ = LatePulseGeneratorPtr(new LatePulseGenerator(randomService_));
    }

    double pmtHV=NAN;
    
    if (!jitterOnly)
    {
        // retrieve the HV for this DOM
        std::map<OMKey, I3DOMStatus>::const_iterator om_stat = status_->domStatus.find(key);
        if (om_stat==status_->domStatus.end())
            log_fatal("No DOMStatus entry found for OMKey(%i,%u)",
                      key.GetString(), key.GetOM());
        pmtHV = om_stat->second.pmtHV;

        if (isnan(pmtHV))
            log_fatal("OMKey(%i,%u): pmtHV is NaN",
                      key.GetString(), key.GetOM());

        if (pmtHV<0.)
            log_fatal("OMKey(%i,%u): pmtHV<0. (value=%gV)",
                      key.GetString(), key.GetOM(), pmtHV/I3Units::V);

        // ignore hits on DOMs with pmtHV==0
        if (pmtHV==0.) return;
    }
    
    // add the input hit to the output vector
    output_vector.push_back(input_hit);
    I3MCHit &output_hit = output_vector.back();
    
    // taken from PPC by D. Chirkin
    const double rnd = randomService_->Uniform();

    double hit_time = input_hit.GetTime();
    
    if (jitterOnly) {
        hit_time += GenerateJitter();
        output_hit.SetTime(hit_time);
    } else {
        if(rnd < pre_pulse_probability_)
        {
            output_hit.SetHitSource(I3MCHit::PRE_PULSE);
            hit_time += pts(pmtHV) + GenerateJitter();
        }
        else if(rnd < pre_pulse_probability_+late_pulse_probability_)
        {
            output_hit.SetHitSource(I3MCHit::ELASTIC_LATE_PULSE);
            hit_time = latePulseGenerator_->GenerateSingleLatePulse(hit_time, pmtHV);
        }
        else
        {
            output_hit.SetHitSource(I3MCHit::SPE);
            hit_time += GenerateJitter();
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
    
}
