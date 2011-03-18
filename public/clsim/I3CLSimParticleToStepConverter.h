#ifndef I3CLSIMPARTICLETOSTEPCONVERTER_H_INCLUDED
#define I3CLSIMPARTICLETOSTEPCONVERTER_H_INCLUDED

#include "icetray/I3TrayHeaders.h"
#include "dataclasses/physics/I3Particle.h"

#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimMediumProperties.h"

#include <boost/noncopyable.hpp>
#include <boost/variant.hpp>

#include <string>
#include <stdexcept>

/**
 * @brief Base class for objects that get a single I3Particle
 * and convert it into individual I3CLSimStep objects stored
 * as bunches in an I3CLSimStepSeries.
 */

class I3CLSimParticleToStepConverter_exception : public std::runtime_error
{
public:
    I3CLSimParticleToStepConverter_exception(const std::string &msg)
    :std::runtime_error(msg)
    {;}
};

struct I3CLSimParticleToStepConverter : private boost::noncopyable
{
public:
    typedef boost::variant<I3CLSimStepSeriesConstPtr, std::pair<uint32_t, I3ParticleConstPtr> > ConversionResult_t;
    
    //virtual ~I3CLSimParticleToStepConverter();

    /**
     * Sets the granularity of the bunch size for the
     * return vectors. The number of entries in a vector
     * will always be a multiple of this number.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetBunchSizeGranularity(uint64_t num) = 0;

    /**
     * Sets the maximum bunch size.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMaxBunchSize(uint64_t num) = 0;

    /**
     * Sets the medium properties.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) = 0;
    
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
     * Adds a new I3Particle to the queue for use as a primary in tracking.
     * The resulting I3CLSimSteps can be retrieved from the
     * I3CLSimParticleToStepConverter after some processing time.
     *
     * Enqueuing a particle after calling EnqueueBarrier 
     * will throw if not all steps have been retrieved.
     *
     * Each step produced by this particles will be tagged
     * with the id set by "identifier".
     * 
     * Will throw if not initialized.
     */
    virtual void EnqueueParticle(const I3Particle &particle, uint32_t identifier) = 0;

    /**
     * Adds a "barrier" to the particle queue. This will keep the
     * simulation running until all steps have been retrieved.
     * New particles can only be added to the queue after all
     * "old" steps have been retrieved.
     * 
     * Will throw if not initialized.
     */
    virtual void EnqueueBarrier() = 0;

    /**
     * Returns true an enqueued barrier is still active. And active
     * barrier means that no new particles can currently be added
     * to the queue. Steps have to be retrieved using GetConversionResult()
     * until this function returns false.
     * 
     * Will throw if not initialized.
     */
    virtual bool BarrierActive() const = 0;
    
    /**
     * Returns true if more steps are available for the current particle.
     * If the return value is false, the current simulation is finished
     * and a new particle may be set.
     * 
     * Will throw if not initialized.
     */
    virtual bool MoreStepsAvailable() const = 0;

    /**
     * Returns a bunch of steps stored in a vector<I3CLSimStep> or
     * a pair of <uint32_t, I3ParticleConstPtr>.
     * Returned particles were produced by the primary particle.
     * The uint32_t value is the primary's identifier as passed to
     * EnqueueParticle().
     *
     * Might block if no steps are available.
     * The steps may belong to various particles at the same time.
     * 
     * Will throw if not initialized.
     */
    virtual ConversionResult_t GetConversionResult(double timeout=NAN) = 0;
    
protected:
};

I3_POINTER_TYPEDEFS(I3CLSimParticleToStepConverter);

#endif //I3CLSIMPARTICLETOSTEPCONVERTER_H_INCLUDED
