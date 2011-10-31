#ifndef I3CLSIMSTEPTOPHOTONCONVERTER_H_INCLUDED
#define I3CLSIMSTEPTOPHOTONCONVERTER_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include "clsim/I3CLSimSimpleGeometry.h"

#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimPhoton.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimRandomValue.h"
#include "clsim/I3CLSimWlenDependentValue.h"

#include <boost/noncopyable.hpp>

#include <map>
#include <string>
#include <stdexcept>

/**
 * @brief Base class for objects that receive a vector of
 * I3CLSimSteps and convert them into I3CLSimPhotons.
 *
 * Depending on the implementation, these photons may
 * either be stored directly at their emission point or
 * after propagation to a target (e.g. a DOM)
 */

class I3CLSimStepToPhotonConverter_exception : public std::runtime_error
{
public:
    I3CLSimStepToPhotonConverter_exception(const std::string &msg)
    :std::runtime_error(msg)
    {;}
};

struct I3CLSimStepToPhotonConverter : private boost::noncopyable
{
public:
    typedef std::pair<uint32_t, I3CLSimPhotonSeriesPtr> ConversionResult_t;
    
    //virtual ~I3CLSimStepToPhotonConverter();

    /**
     * Sets the wavelength generator. By default it should
     * return wavelengths according to a Cherenkov spectrum.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenGenerator(I3CLSimRandomValueConstPtr wlenGenerator) = 0;

    /**
     * Sets the wavelength bias. Set this to a constant value
     * of 1 if you do not need biased photon generation.
     * Photons will be assigned a weight equal to 1/bias.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias) = 0;

    /**
     * Sets the medium properties.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) = 0;

    /**
     * Sets the geometry.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetGeometry(I3CLSimSimpleGeometryConstPtr geometry) = 0;
    
    /**
     * Initializes the simulation.
     * Will throw if already initialized.
     */
    virtual void Initialize() = 0;

    /**
     * Returns true if initialized.
     * Never throws.
     */
    virtual bool IsInitialized() const = 0;
    
    /**
     * Adds a new I3CLSimStepSeries to the queue.
     * The resulting I3CLSimPhotonSeries can be retrieved from the
     * I3CLSimStepToPhotonConverter after some processing time.
     *
     * Enqueuing a vector after calling EnqueueBarrier 
     * will throw if not all photons have been retrieved.
     *
     * Will throw if not initialized.
     */
    virtual void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier) = 0;

    /**
     * Reports the current queue size. The queue works asynchronously,
     * so this value will probably have changed once you use it.
     *
     * Will throw if not initialized.
     */
    virtual std::size_t QueueSize() const = 0; 
    
    /**
     * Returns true if more photons are available.
     * If the return value is false, the current simulation is finished
     * and a new step vector may be set.
     * 
     * Will throw if not initialized.
     */
    virtual bool MorePhotonsAvailable() const = 0;

    /**
     * Returns a bunch of photons stored in a vector<I3CLSimPhoton>.
     *
     * The return value is a pair<uint, vector<I3CLSimPhoton> >.
     * The integer is the same identifier as specified in the call
     * to EnqueueSteps().
     *
     * Might block if no photons are available.
     * 
     * Will throw if not initialized.
     */
    virtual ConversionResult_t GetConversionResult() = 0;
    
protected:
};

I3_POINTER_TYPEDEFS(I3CLSimStepToPhotonConverter);

#endif //I3CLSIMSTEPTOPHOTONCONVERTER_H_INCLUDED
