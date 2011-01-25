#ifndef I3CLSIMPARTICLETOSTEPCONVERTERGEANT4_H_INCLUDED
#define I3CLSIMPARTICLETOSTEPCONVERTERGEANT4_H_INCLUDED

#include "clsim/I3CLSimParticleToStepConverter.h"
#include "dataclasses/physics/I3Particle.h"

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "clsim/I3CLSimQueue.h"

#include <map>
#include <string>

/**
 * @brief A particle-to-step converter using Geant4
 * for tracking.
 */
struct I3CLSimParticleToStepConverterGeant4 : public I3CLSimParticleToStepConverter
{
private:
    static boost::mutex thereCanBeOnlyOneGeant4_mutex;
    static bool thereCanBeOnlyOneGeant4;
    
public:
    typedef std::pair<I3CLSimParticleToStepConverter::ConversionResult_t, bool> FromGeant4Pair_t;
    
    static const uint32_t default_maxQueueItems;
    
    I3CLSimParticleToStepConverterGeant4(uint32_t maxQueueItems=default_maxQueueItems);
    virtual ~I3CLSimParticleToStepConverterGeant4();

    // this class:

    /**
     * Electrons and positrons with a higher enery than this
     * will be stored as secondaries and will not be used
     * for Cherenkov step generation.
     */
    void SetElectronPositronMinEnergyForSecondary(double val);
    
    /**
     * Electrons and positrons with a lower enery than this
     * will be stored as secondaries and will not be used
     * for Cherenkov step generation.
     */
    void SetElectronPositronMaxEnergyForSecondary(double val);
    
    double GetElectronPositronMinEnergyForSecondary() const;
    double GetElectronPositronMaxEnergyForSecondary() const;
    
    
    // inherited:
    
    /**
     * Sets the granularity of the bunch size for the
     * return vectors. The number of entries in a vector
     * will always be a multiple of this number.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetBunchSizeGranularity(uint64_t num);

    /**
     * Sets the maximum bunch size.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMaxBunchSize(uint64_t num);

    /**
     * Sets the medium properties.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties);
    
    /**
     * Initializes the simulation.
     * Will throw if already initialized.
     */
    virtual void Initialize();

    /**
     * Returns true if initialized.
     * Never throws.
     */
    virtual bool IsInitialized() const;
    
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
    virtual void EnqueueParticle(const I3Particle &particle, uint32_t identifier);
    
    /**
     * Adds a "barrier" to the particle queue. This will keep the
     * simulation running until all steps have been retrieved.
     * New particles can only be added to the queue after all
     * "old" steps have been retrieved.
     * 
     * Will throw if not initialized.
     */
    virtual void EnqueueBarrier();
    
    /**
     * Returns true an enqueued barrier is still active. And active
     * barrier means that no new particles can currently be added
     * to the queue. Steps have to be retrieved using GetConversionResult()
     * until this function returns false.
     * 
     * Will throw if not initialized.
     */
    virtual bool BarrierActive() const;
    
    /**
     * Returns true if more steps are available for the current particle.
     * If the return value is false, the current simulation is finished
     * and a new particle may be set.
     * 
     * Will throw if not initialized.
     */
    virtual bool MoreStepsAvailable() const;

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
    virtual I3CLSimParticleToStepConverter::ConversionResult_t GetConversionResult();
    
private:
    void LogGeant4Messages(bool allAsWarn=false) const;

    typedef std::pair<uint32_t, I3ParticleConstPtr> ToGeant4Pair_t;

    void Geant4Thread();
    void Geant4Thread_impl(boost::this_thread::disable_interruption &di);
    boost::shared_ptr<boost::thread> geant4ThreadObj_;
    boost::condition_variable_any geant4Started_cond_;
    boost::mutex geant4Started_mutex_;
    bool geant4Started_;

    mutable boost::mutex barrier_is_enqueued_mutex_;
    bool barrier_is_enqueued_;

    boost::shared_ptr<I3CLSimQueue<ToGeant4Pair_t> > queueToGeant4_;
    boost::shared_ptr<I3CLSimQueue<FromGeant4Pair_t> > queueFromGeant4_;
    mutable boost::shared_ptr<I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > > > queueFromGeant4Messages_;
    
    double electronPositronMinEnergyForSecondary_;
    double electronPositronMaxEnergyForSecondary_;

    bool initialized_;
    uint64_t bunchSizeGranularity_;
    uint64_t maxBunchSize_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
};

I3_POINTER_TYPEDEFS(I3CLSimParticleToStepConverterGeant4);

#endif //I3CLSIMPARTICLETOSTEPCONVERTERGEANT4_H_INCLUDED
