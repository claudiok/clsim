
#ifndef CLSIM_I3CLSIMLIGHTSOURCEPROPAGATOR_H_INCLUDED
#define CLSIM_I3CLSIMLIGHTSOURCEPROPAGATOR_H_INCLUDED

#include "icetray/I3PointerTypedefs.h"

struct I3CLSimStep;
I3_FORWARD_DECLARATION(I3CLSimLightSource);
I3_FORWARD_DECLARATION(I3CLSimMediumProperties);
struct I3CLSimFunction;
typedef boost::shared_ptr<const I3CLSimFunction> I3CLSimFunctionConstPtr;
I3_FORWARD_DECLARATION(I3RandomService);

class I3CLSimLightSourcePropagator {
public:
    
    /// Function to call when a secondary light source is produced. If the
    /// callback returns true, the light source was accepted by another
    /// propagator, and this propagator should stop processing it. Otherwise,
    /// this propagator needs to break the secondary down into steps itself.
    typedef std::function<bool(I3CLSimLightSourceConstPtr &, uint32_t)> secondary_callback;
    
    /// Function to call when a step is produced
    typedef std::function<void(const I3CLSimStep&)> step_callback;
    
    virtual bool IsValidForLightSource(const I3CLSimLightSource &) = 0;
    virtual void Convert(I3CLSimLightSourceConstPtr &, uint32_t, secondary_callback, step_callback) = 0;
    
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
    
    virtual void SetRandomService(I3RandomServicePtr random) = 0;
    
    virtual void Initialize() = 0;
    virtual bool IsInitialized() const = 0;
};

#endif // #ifndef CLSIM_I3CLSIMLIGHTSOURCEPROPAGATOR_H_INCLUDED
