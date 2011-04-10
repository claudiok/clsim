#include "clsim/I3CLSimParticleToStepConverterGeant4.h"

#include "clsim/I3CLSimQueue.h"

#include <boost/thread/locks.hpp>

#include <limits>

// geant4 stuff
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4String.hh"

#include "G4PhysListFactory.hh"
#include "G4ParticleGun.hh"

#include "TrkEMPhysicsUHE.hh"
#include "TrkOpticalPhysics.hh"
#include "TrkEnergyCut.hh"

#include "TrkDetectorConstruction.hh"
#include "TrkPrimaryGeneratorAction.hh"
#include "TrkEventAction.hh"
#include "TrkStackingAction.hh"
//#include "TrkSteppingAction.hh"
#include "TrkUISessionToQueue.hh"

#include "I3CLSimI3ParticleGeantConverter.hh"
#include "clsim/I3CLSimStepStore.h"

#include "Randomize.hh"

// other headers
#include <stdlib.h>


// static definitions

boost::mutex I3CLSimParticleToStepConverterGeant4::thereCanBeOnlyOneGeant4_mutex;
bool I3CLSimParticleToStepConverterGeant4::thereCanBeOnlyOneGeant4=false;

const uint32_t I3CLSimParticleToStepConverterGeant4::default_maxQueueItems=5;
const std::string I3CLSimParticleToStepConverterGeant4::default_physicsListName="QGSP_BERT_EMV";
const double I3CLSimParticleToStepConverterGeant4::default_maxBetaChangePerStep=10.*perCent;
const double I3CLSimParticleToStepConverterGeant4::default_maxNumPhotonsPerStep=200.;


I3CLSimParticleToStepConverterGeant4::I3CLSimParticleToStepConverterGeant4(uint32_t randomSeed,
                                                                           std::string physicsListName,
                                                                           double maxBetaChangePerStep,
                                                                           double maxNumPhotonsPerStep,
                                                                           uint32_t maxQueueItems)
:
queueToGeant4_(new I3CLSimQueue<ToGeant4Pair_t>(0)),
queueFromGeant4_(new I3CLSimQueue<FromGeant4Pair_t>(maxQueueItems)),
queueFromGeant4Messages_(new I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > >(0)), // no maximum size
randomSeed_(randomSeed),
physicsListName_(physicsListName),
maxBetaChangePerStep_(maxBetaChangePerStep),
maxNumPhotonsPerStep_(maxNumPhotonsPerStep),
initialized_(false),
bunchSizeGranularity_(512),
maxBunchSize_(512000)
{
    // Geant4 keeps LOTS of global state and is inherently non-thread-safe.
    // To prevent users from using more than one instance of Geant4 within the
    // same process, we throw an exception if there is already another instance
    // of this class.
    
    {
        boost::unique_lock<boost::mutex> guard(thereCanBeOnlyOneGeant4_mutex);
        
        if (thereCanBeOnlyOneGeant4) 
            throw I3CLSimParticleToStepConverter_exception("There can be only one! ...instance of I3CLSimParticleToStepConverterGeant4.");
        
        thereCanBeOnlyOneGeant4=true;
    }
    
    if ((randomSeed_<0) || (randomSeed_>900000000))
        throw I3CLSimParticleToStepConverter_exception("Invalid random seed (has to be >=0 and <=900000000).");
    
    if ((maxBetaChangePerStep_<=0.) || (maxBetaChangePerStep_>1.))
        throw I3CLSimParticleToStepConverter_exception("Invalid maxBetaChangePerStep.");

    if ((maxNumPhotonsPerStep_<=0.))
        throw I3CLSimParticleToStepConverter_exception("Invalid maxNumPhotonsPerStep.");

    // set up stupid environment variables
    char *dummy = getenv("I3_PORTS");
    const std::string I3_PORTS = dummy?dummy:"";    
    const std::string DATA_BASEDIR = "share/geant4/data";
    
    if (!dummy) {
        G4cout << "The $I3_PORTS variable is not set! You have to provide the Geant4 data directories using their respective environment variables!" << G4endl;
    } else {
        // this does overwrite already existing variables!
        setenv("G4LEVELGAMMADATA", (I3_PORTS + "/" + DATA_BASEDIR + "/PhotonEvaporation2.0").c_str(), 1);
        setenv("G4RADIOACTIVEDATA", (I3_PORTS + "/" + DATA_BASEDIR + "/RadioactiveDecay3.2").c_str(), 1);
        setenv("G4LEDATA", (I3_PORTS + "/" + DATA_BASEDIR + "/G4EMLOW6.9").c_str(), 1);
        setenv("G4NEUTRONHPDATA", (I3_PORTS + "/" + DATA_BASEDIR + "/G4NDL3.13").c_str(), 1);
        setenv("G4ABLADATA", (I3_PORTS + "/" + DATA_BASEDIR + "/G4ABLA3.0").c_str(), 1);
    }
    
    /*
    // check physics list
    {
        G4PhysListFactory *factory = new G4PhysListFactory();
        G4VModularPhysicsList *physics = factory->GetReferencePhysList(physicsListName_.c_str());
        delete factory;
        if (!physics)
            throw I3CLSimParticleToStepConverter_exception("Invalid Physics List name.");
        delete physics;
    } 
     */
}

I3CLSimParticleToStepConverterGeant4::~I3CLSimParticleToStepConverterGeant4()
{
    LogGeant4Messages();

    if (geant4ThreadObj_)
    {
        if (geant4ThreadObj_->joinable())
        {
            log_debug("Stopping the Geant4 thread..");

            geant4ThreadObj_->interrupt();
            
            geant4ThreadObj_->join(); // wait for it indefinitely

            LogGeant4Messages();

            log_debug("Geant4 thread stopped.");
        }
        
        geant4ThreadObj_.reset();
    }
    
    {
        boost::unique_lock<boost::mutex> guard(thereCanBeOnlyOneGeant4_mutex);
        thereCanBeOnlyOneGeant4=false;
    }

    LogGeant4Messages();
}

void I3CLSimParticleToStepConverterGeant4::Initialize()
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 already initialized!");
    
    if (!mediumProperties_)
        throw I3CLSimParticleToStepConverter_exception("MediumProperties not set!");
    
    if (bunchSizeGranularity_ > maxBunchSize_)
        throw I3CLSimParticleToStepConverter_exception("BunchSizeGranularity must not be greater than MaxBunchSize!");
    
    if (maxBunchSize_%bunchSizeGranularity_ != 0)
        throw I3CLSimParticleToStepConverter_exception("MaxBunchSize is not a multiple of BunchSizeGranularity!");
    
    // making a copy of the medium properties
    {
        I3CLSimMediumPropertiesConstPtr copiedMediumProperties(new I3CLSimMediumProperties(*mediumProperties_));
        mediumProperties_ = copiedMediumProperties;
    }
    
    log_info("Starting the Geant4 thread..");
    geant4Started_=false;
    barrier_is_enqueued_=false;

    geant4ThreadObj_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&I3CLSimParticleToStepConverterGeant4::Geant4Thread, this)));

    // wait for startup
    {
        boost::unique_lock<boost::mutex> guard(geant4Started_mutex_);
        for (;;)
        {
            if (geant4Started_) break;
            geant4Started_cond_.wait(guard);
        }
    }        
    
    log_info("Geant4 thread started.");

    LogGeant4Messages();

    initialized_=true;
}

namespace {
    static double GetMaxRIndex(I3CLSimMediumPropertiesConstPtr mProp)
    {
        if (!mProp) return NAN;

        const double minWlen=mProp->GetMinWavelength();
        const double maxWlen=mProp->GetMaxWavelength();
        const unsigned int points=50;
        
        double maxVal=NAN;
        
        for (uint32_t i=0;i<mProp->GetLayersNum();++i)
        {
            I3CLSimWlenDependentValueConstPtr refIndFunc = mProp->GetPhaseRefractiveIndex(i);
            if (!refIndFunc) continue;
            
            for (unsigned int point=0;point<points;++point)
            {
                const double wlen = minWlen + (maxWlen-minWlen)*static_cast<double>(point)/static_cast<double>(points-1);
                const double rindex = refIndFunc->GetValue(wlen);

                //G4cout << "wlen=" << wlen/I3Units::nm << ", rindex=" << rindex << G4endl;
                
                if (isnan(maxVal) || (rindex > maxVal)) maxVal=rindex;
            }
            
        }
        
        return maxVal;
    }
}

#include "G4VStateDependent.hh"

class UserHookForAbortState : public G4VStateDependent
{
public:
    UserHookForAbortState() {;}
    ~UserHookForAbortState() {;}
    
    virtual G4bool Notify(G4ApplicationState requiredState)
    {
        if (requiredState!=G4State_Abort) return true;

        // this seems to be caused by an G4Exception call.
        // throw a real exception here, so we are able to catch it.
        throw std::runtime_error("Geant4 is not amused.");

        return true;
    }
};



void I3CLSimParticleToStepConverterGeant4::Geant4Thread()
{
    // do not interrupt this thread by default
    boost::this_thread::disable_interruption di;

    try {
        Geant4Thread_impl(di);
    } catch(...) { // any exceptions?
        log_warn("Geant4 thread died unexpectedly..");
        LogGeant4Messages(true); // this is maybe the last chance to log..
        
        throw;
    }
}

void I3CLSimParticleToStepConverterGeant4::Geant4Thread_impl(boost::this_thread::disable_interruption &di)
{
    // set up geant4
    
    // this thing stores all the steps generated by Geant4, sorted by the number
    // of Cherenkov photons they generate
    I3CLSimStepStorePtr stepStore(new I3CLSimStepStore( (isnan(maxNumPhotonsPerStep_)||(maxNumPhotonsPerStep_<0.))?0:(static_cast<uint32_t>(maxNumPhotonsPerStep_*1.5)) ));

    // keep this around to "catch" G4Eceptions and throw real exceptions
    UserHookForAbortState *theUserHookForAbortState = new UserHookForAbortState();
    G4StateManager *theStateManager = G4StateManager::GetStateManager();
    theStateManager->RegisterDependent(theUserHookForAbortState);
    
    // get the pointer to the UI manager and set our custom G4cout/G4cerr destination
    G4UImanager* UI = G4UImanager::GetUIpointer();
    TrkUISessionToQueue *theUISessionToQueue = new TrkUISessionToQueue(queueFromGeant4Messages_);
    UI->SetCoutDestination(theUISessionToQueue);

    // initialize the run manager
    G4RunManager* runManager = new G4RunManager;
    
    // set up the "detector" (a lot of water)
    runManager->SetUserInitialization(new TrkDetectorConstruction(mediumProperties_));
    
    // set up the physics list (something+Optical Physics)
    G4PhysListFactory *factory = new G4PhysListFactory();
    G4VModularPhysicsList *physics = factory->GetReferencePhysList(physicsListName_.c_str());
    delete factory;
    
    physics->RegisterPhysics(new TrkOpticalPhysics("Optical",
                                                   maxBetaChangePerStep_,
                                                   maxNumPhotonsPerStep_));
    //physics->RegisterPhysics(new TrkEMPhysicsUHE("EMUHE"));
    //physics->RegisterPhysics(new TrkEnergyCut("TrkEnergyCut"));
    
    physics->SetDefaultCutValue(0.25*mm);
    runManager->SetUserInitialization(physics);
    
    // instantiate some Geant4 helper classes
    TrkPrimaryGeneratorAction *thePrimaryGenerator = new TrkPrimaryGeneratorAction();
    
    runManager->SetUserAction(thePrimaryGenerator);     // runManager now owns this pointer
    runManager->SetUserAction(new TrkStackingAction());   // runManager now owns this pointer

    const double maxRefractiveIndex = GetMaxRIndex(mediumProperties_);
    G4cout << "overall maximum refractive index is " << maxRefractiveIndex << G4endl;
    if (isnan(maxRefractiveIndex)) log_fatal("No maximum refractive index could be found");
    
    TrkEventAction *theEventAction = new TrkEventAction(maxBunchSize_,
                                                        stepStore,
                                                        this->GetParticleParameterizationSeries(),
                                                        queueFromGeant4_,
                                                        di,
                                                        maxRefractiveIndex);
    runManager->SetUserAction(theEventAction);      // runManager now owns this pointer
    
    //UI->ApplyCommand("/physics_engine/tailor/SyncRadiation on");
    //UI->ApplyCommand("/physics_engine/tailor/GammaNuclear on");
    UI->ApplyCommand("/physics_engine/tailor/MuonNuclear on");
    runManager->Initialize();

    //UI->ApplyCommand("/run/particle/dumpCutValues");
    
    CLHEP::HepRandom::setTheSeed(randomSeed_); // value [0,900000000]
    
    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(geant4Started_mutex_);
        geant4Started_=true;
    }
    geant4Started_cond_.notify_all();

    // when flushing the last steps for an event, it may be necessary
    // to include dummy steps to fill the number of steps to be a multiple
    // of bunchSizeGranularity_. This is a step with weight==0 and numPhotons==0,
    // so it should not contribute to the final results!
    I3CLSimStep NoOpStepTemplate;
    NoOpStepTemplate.SetPos(I3Position(0.,0.,0.));
    NoOpStepTemplate.SetDir(I3Direction(0.,0.,-1.));
    NoOpStepTemplate.SetTime(0.);
    NoOpStepTemplate.SetLength(0.);
    NoOpStepTemplate.SetNumPhotons(0);
    NoOpStepTemplate.SetWeight(0.);
    NoOpStepTemplate.SetBeta(1.);
    
    // make a copy of the list of available parameterizations
    const I3CLSimParticleParameterizationSeries parameterizations = this->GetParticleParameterizationSeries();

    // start the main loop
    for (;;)
    {
        I3ParticleConstPtr particle;
        uint32_t particleIdentifier;

        {
            boost::this_thread::restore_interruption ri(di);
            try {
                ToGeant4Pair_t val = queueToGeant4_->Get();
                particleIdentifier = val.first;
                particle = val.second;
            }
            catch(boost::thread_interrupted &i)
            {
                G4cout << "G4 thread was interrupted. closing." << G4endl;
                break;
            }
        }
        
        bool interruptionOccured=false;
        while (stepStore->size() >= maxBunchSize_)
        {
            I3CLSimStepSeriesPtr steps(new I3CLSimStepSeries());
            stepStore->pop_bunch_to_vector(maxBunchSize_, *steps);
            
            {
                boost::this_thread::restore_interruption ri(di);
                try {
                    queueFromGeant4_->Put(std::make_pair(steps, false));
                } catch(boost::thread_interrupted &i) {
                    G4cout << "G4 thread was interrupted. shutting down Geant4!" << G4endl;
                    
                    interruptionOccured = true;
                    break;
                }
            }
            
            //G4cout << " -> flushed " << maxBunchSize_ << " steps, " << stepStore->size() << " steps left" << G4endl;
        }            
        
        if (interruptionOccured) break;
        
        if (!particle) {
            //G4cout << "G4 thread got NULL! flushing " << stepStore->size() << " steps." << G4endl;

            if (stepStore->empty()) {
                // nothing to send. send an empty step vector along with
                // the command to disable the barrier
                
                I3CLSimStepSeriesPtr steps(new I3CLSimStepSeries());

                {
                    boost::this_thread::restore_interruption ri(di);
                    try {
                        queueFromGeant4_->Put(std::make_pair(steps, true /* this is the LAST reply before the barrier is reached! */));
                    } catch(boost::thread_interrupted &i) {
                        G4cout << "G4 thread was interrupted. closing." << G4endl;
                        break;
                    }

                    
                    //boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
                    //barrier_is_enqueued_=false;
                }
                
                continue; // nothing to send
            }
            
            
            // send fully-sized bunches first
            
            bool interruptionOccured=false;
            while (stepStore->size() >= maxBunchSize_)
            {
                I3CLSimStepSeriesPtr steps(new I3CLSimStepSeries());
                stepStore->pop_bunch_to_vector(maxBunchSize_, *steps);
                
                {
                    boost::this_thread::restore_interruption ri(di);
                    try {
                        queueFromGeant4_->Put(std::make_pair(steps, false));
                    } catch(boost::thread_interrupted &i) {
                        G4cout << "G4 thread was interrupted. shutting down Geant4!" << G4endl;
                        
                        interruptionOccured = true;
                        break;
                    }
                }
                
                //G4cout << " -> flushed " << maxBunchSize_ << " steps, " << stepStore->size() << " steps left" << G4endl;
            }            
            
            if (interruptionOccured) break;
            
            
            // flush the rest (size < full bunch size)
            
            I3CLSimStepSeriesPtr steps(new I3CLSimStepSeries());
            const std::size_t numStepsWithDummyFill = bunchSizeGranularity_>1?(((stepStore->size()/bunchSizeGranularity_)+1)*bunchSizeGranularity_):stepStore->size();

            //G4cout << " -> " << stepStore->size() << " steps left, padding to " << numStepsWithDummyFill << G4endl;
            
            stepStore->pop_bunch_to_vector(numStepsWithDummyFill, *steps, NoOpStepTemplate);
            
            if (!stepStore->empty())
                log_fatal("Internal logic error. step store should be empty.");
            
            {
                boost::this_thread::restore_interruption ri(di);
                try {
                    queueFromGeant4_->Put(std::make_pair(steps, true /* this is the LAST reply before the barrier is reached! */));
                } catch(boost::thread_interrupted &i) {
                    G4cout << "G4 thread was interrupted. closing." << G4endl;
                    break;
                }
            }
            
            //G4cout << " -> flush complete." << G4endl;

            // nothing to send to Geant4, so start from the beginning
            continue;
        }
        
        //G4cout << "G4 thread got particle id " << particleIdentifier << ", type: " << particle->GetTypeString() << G4endl;
        
        // check if there is a parameterization for this particle, so we
        // may not even have to send it to Geant4
        {
            bool parameterizationIsAvailable=false;
            
            for (I3CLSimParticleParameterizationSeries::const_iterator it=parameterizations.begin();
                 it!=parameterizations.end(); ++it)
            {
                const I3CLSimParticleParameterization &parameterization = *it;
                
                if (parameterization.IsValidForParticle(*particle))
                {
                    parameterizationIsAvailable=true;
                    
                    G4cout << "Geant4: sending a " << particle->GetTypeString() << " with id " << particleIdentifier << " and E=" << particle->GetEnergy()/I3Units::GeV << "GeV to a parameterization handler." << G4endl;

                    // call the converter
                    if (!parameterization.converter) log_fatal("Internal error: parameteriation has NULL converter");
                    if (!parameterization.converter->IsInitialized()) log_fatal("Internal error: parameterization converter is not initialized.");
                    if (parameterization.converter->BarrierActive()) log_fatal("Logic error: parameterization converter has active barrier.");
                    
                    parameterization.converter->EnqueueParticle(*particle, particleIdentifier);
                    parameterization.converter->EnqueueBarrier();
                    
                    while (parameterization.converter->BarrierActive())
                    {
                        I3CLSimStepSeriesConstPtr res =
                        parameterization.converter->GetConversionResult();
                        
                        if (!res) {
                            log_warn("NULL result from parameterization GetConversionResult(). ignoring.");
                            continue;
                        }
                        if (res->size()==0) continue; // ignore empty vectors
                        
                        BOOST_FOREACH(const I3CLSimStep &step, *res)
                        {
                            stepStore->insert_copy(step.GetNumPhotons(), step);
                        }
                    }
                    
                    
                    
                    break;
                }
            }
            
            if (parameterizationIsAvailable) continue;

        }
        
        // no parameterization was found, use default Geant4
        
        // configure the Geant4 particle gun
        {
            G4ParticleGun *particleGun = thePrimaryGenerator->GetParticleGun();
            if (!particleGun) log_fatal("Internal error: G4ParticleGun instance is NULL!");
            
            const bool ret = I3CLSimI3ParticleGeantConverter::SetParticleGun(particleGun, *particle);
            
            if (!ret) {
                G4cerr << "Could not configure Geant4 to shoot a " << particle->GetTypeString() << "! Ignoring." << G4endl;

                continue;
            }
        }
        
        // set the current particle ID
        theEventAction->SetExternalParticleID(particleIdentifier);

        log_debug("Geant4: no parameterization for %s with E=%fGeV", particle->GetTypeString().c_str(), particle->GetEnergy()/I3Units::GeV);
        
        G4cout << "Geant4: shooting a " << particle->GetTypeString() << " with id " << particleIdentifier << " and E=" << particle->GetEnergy()/I3Units::GeV << "GeV." << G4endl;

        // we don't need the particle anymore.
        particle.reset();
        
        // turn on the Geant4 beam!
        runManager->BeamOn(1);
        
        // check if AbortRun was requested beacause of a thread interruption.
        if (theEventAction->AbortWasRequested()) break;
        
    }

    G4cout << "G4 thread terminating..." << G4endl;

    UI->SetCoutDestination(NULL);
    delete theUISessionToQueue;
    
    // job termination
    delete runManager;

    log_debug("G4 thread terminated.");
}

bool I3CLSimParticleToStepConverterGeant4::IsInitialized() const
{
    LogGeant4Messages();

    return initialized_;
}

void I3CLSimParticleToStepConverterGeant4::SetBunchSizeGranularity(uint64_t num)
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 already initialized!");
    
    if (num<=0)
        throw I3CLSimParticleToStepConverter_exception("BunchSizeGranularity of 0 is invalid!");
    
    bunchSizeGranularity_=num;
}

void I3CLSimParticleToStepConverterGeant4::SetMaxBunchSize(uint64_t num)
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 already initialized!");

    if (num<=0)
        throw I3CLSimParticleToStepConverter_exception("MaxBunchSize of 0 is invalid!");

    maxBunchSize_=num;
}

void I3CLSimParticleToStepConverterGeant4::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 already initialized!");

    mediumProperties_=mediumProperties;
}

void I3CLSimParticleToStepConverterGeant4::EnqueueParticle(const I3Particle &particle, uint32_t identifier)
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
        if (barrier_is_enqueued_)
            throw I3CLSimParticleToStepConverter_exception("A barrier is enqueued! You must receive all steps before enqueuing a new particle.");
    }
    
    I3ParticleConstPtr particleCopy(new I3Particle(particle));
    queueToGeant4_->Put(make_pair(identifier, particleCopy));
    
    LogGeant4Messages();
}

void I3CLSimParticleToStepConverterGeant4::EnqueueBarrier()
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
        if (barrier_is_enqueued_)
            throw I3CLSimParticleToStepConverter_exception("A barrier is already enqueued!");
        
        barrier_is_enqueued_=true;

        // we use a NULL pointer as the barrier
        queueToGeant4_->Put(make_pair(0, I3ParticleConstPtr()));
    }
    
    LogGeant4Messages();
}

bool I3CLSimParticleToStepConverterGeant4::BarrierActive() const
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);

        return barrier_is_enqueued_;
    }
}

bool I3CLSimParticleToStepConverterGeant4::MoreStepsAvailable() const
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 is not initialized!");

    return (!queueFromGeant4_->empty());
}

I3CLSimStepSeriesConstPtr I3CLSimParticleToStepConverterGeant4::GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterGeant4 is not initialized!");

    barrierWasReset=false;
    
    FromGeant4Pair_t ret;
    if (!isnan(timeout))
        ret = queueFromGeant4_->Get(timeout/I3Units::second, FromGeant4Pair_t(I3CLSimStepSeriesConstPtr(), false));
    else
        ret = queueFromGeant4_->Get(); // no timeout
    
    if (ret.second)
    {
        {
            boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
            if (!barrier_is_enqueued_)
                log_error("Internal logic error. Barrier is not set as enqueued, yet a finalization message was received.");
            barrierWasReset=true;
            barrier_is_enqueued_=false;
        }
    }
    
    LogGeant4Messages();
    
    return ret.first;
}

void I3CLSimParticleToStepConverterGeant4::LogGeant4Messages(bool allAsWarn) const
{
    if (!queueFromGeant4Messages_) return;
    
    while (!queueFromGeant4Messages_->empty())
    {
        boost::shared_ptr<std::pair<const std::string, bool> > str = queueFromGeant4Messages_->Get();
        
        if (!str) 
        {
            log_warn("Geant4 said: (null)");
        }
        else
        {
            std::string out;
            if (str->first.size()==0) {
                out = str->first;
            } else if (str->first[str->first.size()-1] == '\n') {
                out = str->first.substr(0, str->first.size()-1);
            } else {
                out = str->first;
            }
            
            if (str->second) {
                log_warn("Geant4 warns: %s", out.c_str());
            } else {
                if (allAsWarn)
                    log_warn("Geant4 says:  %s", out.c_str());
                else
                    log_info("Geant4 says:  %s", out.c_str());
            }
        }
    }
}

