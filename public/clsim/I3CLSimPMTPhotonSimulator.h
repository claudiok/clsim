#ifndef I3CLSIMPMTPHOTONSIMULATOR_H_INCLUDED
#define I3CLSIMPMTPHOTONSIMULATOR_H_INCLUDED

#include "icetray/I3TrayHeaders.h"
#include "dataclasses/physics/I3MCHit.h"

#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/status/I3DetectorStatus.h"

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
     * This function does not clear the vector before filling it. New hits are added to it.
     */
    virtual void ApplyAfterPulseLatePulseAndJitterSim(const OMKey &key,
                                                      const I3MCHit &input_hit,
                                                      std::vector<I3MCHit> &output_vector) const = 0;
    
private:
};


I3_POINTER_TYPEDEFS(I3CLSimPMTPhotonSimulator);

#endif //I3CLSIMPMTPHOTONSIMULATOR_H_INCLUDED
