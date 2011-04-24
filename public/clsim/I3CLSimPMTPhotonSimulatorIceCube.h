#ifndef I3CLSIMPMTPHOTONSIMULATORICECUBE_H_INCLUDED
#define I3CLSIMPMTPHOTONSIMULATORICECUBE_H_INCLUDED

#include "clsim/I3CLSimPMTPhotonSimulator.h"

#include "icetray/I3Units.h"
#include "phys-services/I3RandomService.h"

// froward declarations
class AfterPulseGenerator;
class LatePulseGenerator;

/**
 * @brief Simulate PMT jitter, late-pulses and
 * after-pulses by using code from hit-maker.
 *
 */
struct I3CLSimPMTPhotonSimulatorIceCube : public I3CLSimPMTPhotonSimulator
{
public:
    static const double DEFAULT_jitter;
    static const double DEFAULT_pre_pulse_probability;
    static const double DEFAULT_late_pulse_probability;
    static const double DEFAULT_after_pulse_probability;
    
    I3CLSimPMTPhotonSimulatorIceCube(I3RandomServicePtr randomService,
                                     double jitter=DEFAULT_jitter,
                                     double pre_pulse_probability=DEFAULT_pre_pulse_probability,
                                     double late_pulse_probability=DEFAULT_late_pulse_probability,
                                     double after_pulse_probability=DEFAULT_after_pulse_probability);
    virtual ~I3CLSimPMTPhotonSimulatorIceCube();

    /**
     * Set the current calibration
     */
    virtual void SetCalibration(I3CalibrationConstPtr calibration);
    
    /**
     * Set the current status
     */
    virtual void SetDetectorStatus(I3DetectorStatusConstPtr status);

    /**
     * This function does not clear the vector before filling it. New hits are added to it.
     */
    virtual void ApplyAfterPulseLatePulseAndJitterSim(const OMKey &key,
                                                      const I3MCHit &input_hit,
                                                      std::vector<I3MCHit> &output_vector) const;
    
private:
    I3RandomServicePtr randomService_;
    double jitter_;
    double pre_pulse_probability_;
    double late_pulse_probability_;
    double after_pulse_probability_;
    
    shared_ptr<AfterPulseGenerator> afterPulseGenerator_;
    shared_ptr<LatePulseGenerator> latePulseGenerator_;

    I3CalibrationConstPtr calibration_;
    I3DetectorStatusConstPtr status_;
};


I3_POINTER_TYPEDEFS(I3CLSimPMTPhotonSimulatorIceCube);

#endif //I3CLSIMPMTPHOTONSIMULATORICECUBE_H_INCLUDED
