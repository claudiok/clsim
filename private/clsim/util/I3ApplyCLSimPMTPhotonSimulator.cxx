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
 * @file I3ApplyCLSimPMTPhotonSimulator.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "phys-services/I3RandomService.h"

#include "clsim/function/I3CLSimFunction.h"
#include "clsim/I3CLSimPMTPhotonSimulator.h"

#include <boost/foreach.hpp>

#include "dataclasses/physics/I3MCHit.h"

/**
 * Applies jitter, after-/late- and re-pulses as described
 * by a I3CLSimPMTPhotonSimulator object.
 */
class I3ApplyCLSimPMTPhotonSimulator : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3ApplyCLSimPMTPhotonSimulator(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    virtual ~I3ApplyCLSimPMTPhotonSimulator();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
	virtual void Calibration(I3FramePtr frame);
	
	virtual void DetectorStatus(I3FramePtr frame);
	
    /**
     * The module needs to process Physics frames
     */
#ifdef IS_Q_FRAME_ENABLED
    virtual void DAQ(I3FramePtr frame);
#else
    virtual void Physics(I3FramePtr frame);
#endif
    
private:
    // parameters

    /// Parameter: The PMT simulation worker object    
    I3CLSimPMTPhotonSimulatorPtr pmtPhotonSimulator_;

    /// Parameter: Name of the input I3MCHitSeriesMap frame object.
    std::string inputMCHitSeriesMapName_;

    /// Parameter: Name of the output I3MCHitSeriesMap frame object.
    std::string outputMCHitSeriesMapName_;
    

    /// a random service
    I3RandomServicePtr randomService_;
    
private:
    // default, assignment, and copy constructor declared private
    I3ApplyCLSimPMTPhotonSimulator();
    I3ApplyCLSimPMTPhotonSimulator(const I3ApplyCLSimPMTPhotonSimulator&);
    I3ApplyCLSimPMTPhotonSimulator& operator=(const I3ApplyCLSimPMTPhotonSimulator&);
    
    SET_LOGGER("I3ApplyCLSimPMTPhotonSimulator");
};



// The module
I3_MODULE(I3ApplyCLSimPMTPhotonSimulator);

I3ApplyCLSimPMTPhotonSimulator::I3ApplyCLSimPMTPhotonSimulator(const I3Context& context) 
: I3ConditionalModule(context)
{
    AddParameter("PMTPhotonSimulator",
                 "Optional after-pulse, late-pulse and jitter simulator object,\n"
                 "an instance of I3CLSimPMTPhotonSimulator",
                 pmtPhotonSimulator_);

    inputMCHitSeriesMapName_="I3MCTree";
    AddParameter("InputMCHitSeriesMapName",
                 "Name of the input I3MCHitSeriesMap frame object.",
                 inputMCHitSeriesMapName_);

    outputMCHitSeriesMapName_="I3MCTree";
    AddParameter("OutputMCHitSeriesMapName",
                 "Name of the output I3MCHitSeriesMap frame object.",
                 outputMCHitSeriesMapName_);

    AddParameter("RandomService",
                 "An instance of I3RandomService.",
                 randomService_);

    // add an outbox
    AddOutBox("OutBox");

}

I3ApplyCLSimPMTPhotonSimulator::~I3ApplyCLSimPMTPhotonSimulator()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3ApplyCLSimPMTPhotonSimulator::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("PMTPhotonSimulator", pmtPhotonSimulator_);
    GetParameter("InputMCHitSeriesMapName", inputMCHitSeriesMapName_);
    GetParameter("OutputMCHitSeriesMapName", outputMCHitSeriesMapName_);
    randomService_ = context_.Get<I3RandomServicePtr>("I3RandomService");

    if (inputMCHitSeriesMapName_=="")
        log_fatal("The \"InputMCHitSeriesMapName\" parameter must not be empty.");

    if (outputMCHitSeriesMapName_=="")
        log_fatal("The \"OutputMCHitSeriesMapName\" parameter must not be empty.");

    if (!pmtPhotonSimulator_)
        log_fatal("The \"OutputMCHitSeriesMapName\" parameter must not be empty.");
	
	pmtPhotonSimulator_->SetRandomService(randomService_);
}


void I3ApplyCLSimPMTPhotonSimulator::Calibration(I3FramePtr frame){
	pmtPhotonSimulator_->SetCalibration(frame->Get<boost::shared_ptr<const I3Calibration> >());
	PushFrame(frame);
}

void I3ApplyCLSimPMTPhotonSimulator::DetectorStatus(I3FramePtr frame){
	pmtPhotonSimulator_->SetDetectorStatus(frame->Get<boost::shared_ptr<const I3DetectorStatus> >());
	PushFrame(frame);
}

#ifdef IS_Q_FRAME_ENABLED
void I3ApplyCLSimPMTPhotonSimulator::DAQ(I3FramePtr frame)
#else
void I3ApplyCLSimPMTPhotonSimulator::Physics(I3FramePtr frame)
#endif
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3MCHitSeriesMapConstPtr inputHits = frame->Get<I3MCHitSeriesMapConstPtr>(inputMCHitSeriesMapName_);
    if (!inputHits) {
        log_debug("Frame does not contain an I3MCHitSeriesMap named \"%s\".",
                  inputMCHitSeriesMapName_.c_str());
        PushFrame(frame);
        return;
    }

    // allocate the output I3MCHitSeriesMap
    I3MCHitSeriesMapPtr outputHits(new I3MCHitSeriesMap());

    for (I3MCHitSeriesMap::const_iterator it=inputHits->begin();
         it!=inputHits->end(); ++it)
    {
        const OMKey &key = it->first;
        const I3MCHitSeries &inputVector = it->second;

        std::vector<I3MCHit> &outputVector = 
            outputHits->insert(std::make_pair(key, I3MCHitSeries() )).first->second;

        BOOST_FOREACH(const I3MCHit &input_hit, inputVector)
        {
#ifdef I3MCHIT_WEIGHT_IS_DEPRECATED
            const std::size_t numberOfHits = input_hit.GetNPE();
#else
            const std::size_t numberOfHits = input_hit.GetWeight();
            if (static_cast<double>(numberOfHits) != input_hit.GetWeight())
            {
                log_fatal("only integer number of pulses are supported. You have a hit with weight %f.",
                    input_hit.GetWeight());
            }

#endif

            if (numberOfHits==0)
                log_warn("You have a hit with NPE==0");

            if ((input_hit.GetHitSource() != I3MCHit::SPE) && 
                (input_hit.GetHitSource() != I3MCHit::RANDOM))
                log_error("You have a hit where HitSource is neither SPE nor RANDOM.");

            // insert N copies of the original hit, each with NPE=1
            for (std::size_t i=0;i<numberOfHits;++i)
            {
                // let pmtPhotonSimulator add the hit and all afterpulses
                I3MCHit hit = input_hit;

                // fill in all information
                hit.SetTime(input_hit.GetTime());
                hit.SetHitID(input_hit.GetHitID());
#ifdef I3MCHIT_WEIGHT_IS_DEPRECATED
                hit.SetNPE(1);
#else
                hit.SetWeight(1.0);
#endif
                hit.SetHitSource(I3MCHit::SPE); // SPE for now, may be changed by afterpulse simulation

                pmtPhotonSimulator_->ApplyAfterPulseLatePulseAndJitterSim
                    (key, hit, outputVector);
            }

        }

    }



    frame->Put(outputMCHitSeriesMapName_, outputHits);
    
    // that's it!
    PushFrame(frame);
}
