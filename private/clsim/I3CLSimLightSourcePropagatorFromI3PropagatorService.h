
#ifndef CLSIM_I3CLSIMLIGHTSOURCEPROPAGATORFROMI3PROPAGATORSERVICE_H_INCLUDED
#define CLSIM_I3CLSIMLIGHTSOURCEPROPAGATORFROMI3PROPAGATORSERVICE_H_INCLUDED

#include <clsim/I3CLSimLightSourcePropagator.h>
#include <sim-services/I3PropagatorService.h>

class I3CLSimLightSourcePropagatorFromI3PropagatorService : public I3CLSimLightSourcePropagator
{
public:
    
    I3CLSimLightSourcePropagatorFromI3PropagatorService(I3ParticleTypePropagatorServiceMapPtr);
    virtual ~I3CLSimLightSourcePropagatorFromI3PropagatorService();

    // inherited:

    virtual void SetRandomService(I3RandomServicePtr random);

    /**
     * Sets the wavelength bias. Set this to a constant value
     * of 1 if you do not need biased photon generation.
     * The Cherenkov spectrum will be multiplied by this
     * value at each wavelength.
     * This will influence the number of photons produced.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias) {};

    /**
     * Sets the medium properties.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) {};
    
    /**
     * Initializes the simulation.
     * Will throw if already initialized.
     */
    virtual void Initialize() { initialized_=true; };

    /**
     * Returns true if initialized.
     * Never throws.
     */
    virtual bool IsInitialized() const { return initialized_; };
    
    virtual bool IsValidForLightSource(const I3CLSimLightSource &source);
    virtual void Convert(I3CLSimLightSourceConstPtr &, uint32_t, secondary_callback, step_callback);
    
private:
    I3ParticleTypePropagatorServiceMapPtr particleToPropagatorServiceMap_;

    bool initialized_;
    
    SET_LOGGER("I3CLSimLightSourcePropagatorFromI3PropagatorService");
};
 

#endif // CLSIM_I3CLSIMLIGHTSOURCEPROPAGATORFROMI3PROPAGATORSERVICE_H_INCLUDED