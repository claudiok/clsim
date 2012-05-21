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
 * @file I3CLSimLightSourceToStepConverterFlasher.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMLIGHTSOURCETOSTEPCONVERTERFLASHER_H_INCLUDED
#define I3CLSIMLIGHTSOURCETOSTEPCONVERTERFLASHER_H_INCLUDED

#include "clsim/I3CLSimLightSourceToStepConverter.h"
#include "clsim/I3CLSimLightSource.h"
#include "clsim/function/I3CLSimFunction.h"
#include "clsim/random_value/I3CLSimRandomValue.h"
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
    static const bool default_interpretAngularDistributionsInPolarCoordinates;

    /**
     * Initializes a new converter object for a specific flasher type.
     *
     * flasherSpectrumNoBias: the flasher wavelength spectrum without any
     *   bias factors for DOM acceptance.
     *
     * spectrumTable: the global spectrum table used for OpenCL code
     *   generation.
     *
     * angularProfileDistributionPolar: a random number distribution of
     *   polar angles w.r.t. to the flasher orientation used for smearing
     *   of initial photon directions. The distribution is assumed to take
     *   a single run-time parameter: the "width" as specified in I3CLSimFlasherPulse.
     *
     * angularProfileDistributionAzimuthal: a random number distribution of
     *   azimuthal angles w.r.t. to the flasher orientation used for smearing
     *   of initial photon directions. The distribution is assumed to take
     *   a single run-time parameter: the "width" as specified in I3CLSimFlasherPulse.
     *
     * timeDelayDistribution: a distribution of time delays w.r.t. the nominal
     *   flasher time. Values from this distribution will be added to the nominal
     *   flasher time to generate photon starting times. The distribution is assumed to take
     *   a single run-time parameter: the "width" as specified in I3CLSimFlasherPulse.
     *
     * interpretAngularDistributionsInPolarCoordinates: interpret the
     *   angularProfileDistributionPolar and angularProfileDistributionAzimuthal
     *   distributions in polar coordinates. Instead of a shift in the polar
     *   and azimuthal directions, the polar shift is a shift perpendicular from the
     *   current direction. The azimuthal shift is a rotation around the old axis.
     *   Since the origin of this rotation is undefined, it is best to use a uniform
     *   range of angles from 0deg to 360deg for the azimuthal shift in this case.
     *
     * photonsPerStep: instead of creating photons directly, this class generates
     *   so-called steps (bunches of photons described by a single object).
     *   This controls how many photons there will be per step (some steps may have
     *   less photons than specified here to accomodate all possible numbers of
     *   photons).
     *   All random distributions are only sampled once per step, not once per photon,
     *   thus having too many photons per step may introduce too much grnaularity.
     *   This should not be a problem for flashers since they generate a large number
     *   of photons anyway.
     * 
     */
    I3CLSimLightSourceToStepConverterFlasher(I3CLSimFunctionConstPtr flasherSpectrumNoBias,
                                             I3CLSimSpectrumTablePtr spectrumTable,
                                             I3CLSimRandomValueConstPtr angularProfileDistributionPolar,
                                             I3CLSimRandomValueConstPtr angularProfileDistributionAzimuthal,
                                             I3CLSimRandomValueConstPtr timeDelayDistribution,
                                             bool interpretAngularDistributionsInPolarCoordinates=default_interpretAngularDistributionsInPolarCoordinates,
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
    
    I3CLSimRandomValueConstPtr angularProfileDistributionPolar_;
    I3CLSimRandomValueConstPtr angularProfileDistributionAzimuthal_;
    I3CLSimRandomValueConstPtr timeDelayDistribution_;

    bool interpretAngularDistributionsInPolarCoordinates_;
    
};

I3_POINTER_TYPEDEFS(I3CLSimLightSourceToStepConverterFlasher);

#endif //I3CLSIMLIGHTSOURCETOSTEPCONVERTERFLASHER_H_INCLUDED
