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
 * @file I3CLSimPMTPhotonSimulator.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMPMTPHOTONSIMULATOR_H_INCLUDED
#define I3CLSIMPMTPHOTONSIMULATOR_H_INCLUDED

#include "icetray/I3TrayHeaders.h"
#include "dataclasses/physics/I3MCHit.h"

#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/status/I3DetectorStatus.h"

#include "phys-services/I3RandomService.h"

#include <vector>

/**
 * @brief Base class. Classes derived from this are used
 * for late-pulse, after-pulse and jitter simulations.
 *
 * They provide a single function, receiving one I3MCHit
 * and filling a vector of hits, contining the original
 * hit with jitter and (possibly) some after-pulses.
 *
 */
struct I3CLSimPMTPhotonSimulator 
{
public:
    
    I3CLSimPMTPhotonSimulator();
    virtual ~I3CLSimPMTPhotonSimulator();

    /**
     * Set the current calibration
     */
    virtual void SetCalibration(I3CalibrationConstPtr calibration) = 0;

    /**
     * Set the current status
     */
    virtual void SetDetectorStatus(I3DetectorStatusConstPtr status) = 0;

    /**
     * Sets the random number generator service.
     * This should be an instance of I3RandomService.
     */
    virtual void SetRandomService(I3RandomServicePtr random) = 0;

    /**
     * This function does not clear the vector before filling it. New hits are added to it.
     */
    virtual void ApplyAfterPulseLatePulseAndJitterSim(const OMKey &key,
                                                      const I3MCHit &input_hit,
                                                      std::vector<I3MCHit> &output_vector) const = 0;
    
private:
};


I3_POINTER_TYPEDEFS(I3CLSimPMTPhotonSimulator);

#endif //I3CLSIMPMTPHOTONSIMULATOR_H_INCLUDED
