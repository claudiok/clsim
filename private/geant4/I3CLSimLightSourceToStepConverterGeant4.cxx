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
 * @file I3CLSimLightSourceToStepConverterGeant4.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

#include "clsim/I3CLSimQueue.h"

#include <boost/thread/locks.hpp>
#include <boost/foreach.hpp>

#include <limits>
#include <deque>
#include <boost/tuple/tuple.hpp>

#ifdef HAS_GEANT4
// geant4 stuff
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4String.hh"

#include "G4PhysListFactory.hh"
#include "G4ParticleGun.hh"

#include "G4Version.hh"

#include "TrkEMPhysicsUHE.hh"
#include "TrkOpticalPhysics.hh"

#include "TrkDetectorConstruction.hh"
#include "TrkPrimaryGeneratorAction.hh"
#include "TrkEventAction.hh"
#include "TrkStackingAction.hh"
//#include "TrkSteppingAction.hh"
#include "TrkUISessionToQueue.hh"

#include "I3CLSimI3ParticleGeantConverter.hh"
#endif

#include "clsim/I3CLSimStepStore.h"

#ifdef HAS_GEANT4
#include "Randomize.hh"
#endif

// other headers
#include <stdlib.h>


// static definitions

boost::mutex I3CLSimLightSourceToStepConverterGeant4::thereCanBeOnlyOneGeant4_mutex;
bool I3CLSimLightSourceToStepConverterGeant4::thereCanBeOnlyOneGeant4=false;

const uint32_t I3CLSimLightSourceToStepConverterGeant4::default_maxQueueItems=5;
const std::string I3CLSimLightSourceToStepConverterGeant4::default_physicsListName="QGSP_BERT_EMV";
const double I3CLSimLightSourceToStepConverterGeant4::default_maxBetaChangePerStep=10.*I3Units::perCent;
const uint32_t I3CLSimLightSourceToStepConverterGeant4::default_maxNumPhotonsPerStep=200;
#ifdef HAS_GEANT4
const bool I3CLSimLightSourceToStepConverterGeant4::canUseGeant4=true;
#else
const bool I3CLSimLightSourceToStepConverterGeant4::canUseGeant4=false;
#endif


I3CLSimLightSourceToStepConverterGeant4::I3CLSimLightSourceToStepConverterGeant4(std::string physicsListName,
                                                                           double maxBetaChangePerStep,
                                                                           uint32_t maxNumPhotonsPerStep,
                                                                           uint32_t maxQueueItems)
:
queueToGeant4_(new I3CLSimQueue<ToGeant4Pair_t>(0)),
queueFromGeant4_(new I3CLSimQueue<FromGeant4Pair_t>(maxQueueItems)),
queueFromGeant4Messages_(new I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > >(0)), // no maximum size
physicsListName_(physicsListName),
maxBetaChangePerStep_(maxBetaChangePerStep),
maxNumPhotonsPerStep_(maxNumPhotonsPerStep),
initialized_(false),
bunchSizeGranularity_(512),
maxBunchSize_(512000)
{
    // Geant4 keeps LOTS of global state and is inherently non-thread-safe.
    // To prevent users from using more than one instance of Geant4 within the
    // same process, we throw an exception in case another instance
    // of this class already exists.
    
    {
        boost::unique_lock<boost::mutex> guard(thereCanBeOnlyOneGeant4_mutex);
        
        if (thereCanBeOnlyOneGeant4) 
            throw I3CLSimLightSourceToStepConverter_exception("There can be only one! ...instance of I3CLSimLightSourceToStepConverterGeant4.");
        
        thereCanBeOnlyOneGeant4=true;
    }
    
    if ((maxBetaChangePerStep_<=0.) || (maxBetaChangePerStep_>1.))
        throw I3CLSimLightSourceToStepConverter_exception("Invalid maxBetaChangePerStep.");

    if ((maxNumPhotonsPerStep_<=0.))
        throw I3CLSimLightSourceToStepConverter_exception("Invalid maxNumPhotonsPerStep.");

    // check for the braindead Geant4 environment variables
    if ((!getenv("G4LEVELGAMMADATA")) ||
        (!getenv("G4RADIOACTIVEDATA")) ||
        (!getenv("G4LEDATA")) ||
        (!getenv("G4NEUTRONHPDATA")) ||
        (!getenv("G4ABLADATA")))
    {
        log_info("Geant4 requires the following environment variables to be set: \"G4LEVELGAMMADATA\", \"G4RADIOACTIVEDATA\", \"G4LEDATA\", \"G4NEUTRONHPDATA\" and \"G4ABLADATA\"");
    }

    // geant 4.9.5 needs some more
#if G4VERSION_NUMBER >= 950
    if ((!getenv("G4NEUTRONXSDATA")) ||
        (!getenv("G4PIIDATA")) ||
        (!getenv("G4REALSURFACEDATA")))
    {
        log_info("Geant4.9.5 requires the following environment variables to be set: \"G4NEUTRONXSDATA\", \"G4PIIDATA\", \"G4REALSURFACEDATA\"");
    }
#endif
    
}

I3CLSimLightSourceToStepConverterGeant4::~I3CLSimLightSourceToStepConverterGeant4()
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

void I3CLSimLightSourceToStepConverterGeant4::Initialize()
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 already initialized!");

    if (!randomService_)
        throw I3CLSimLightSourceToStepConverter_exception("RandomService not set!");

    if (!wlenBias_)
        throw I3CLSimLightSourceToStepConverter_exception("WlenBias not set!");

    if (!mediumProperties_)
        throw I3CLSimLightSourceToStepConverter_exception("MediumProperties not set!");
    
    if (bunchSizeGranularity_ > maxBunchSize_)
        throw I3CLSimLightSourceToStepConverter_exception("BunchSizeGranularity must not be greater than MaxBunchSize!");
    
    if (maxBunchSize_%bunchSizeGranularity_ != 0)
        throw I3CLSimLightSourceToStepConverter_exception("MaxBunchSize is not a multiple of BunchSizeGranularity!");
    
    // make sure none of the parameterizations are initialized
    const I3CLSimLightSourceParameterizationSeries &parameterizations = this->GetLightSourceParameterizationSeries();
    for (I3CLSimLightSourceParameterizationSeries::const_iterator it=parameterizations.begin();
         it!=parameterizations.end(); ++it)
    {
        const I3CLSimLightSourceParameterization &parameterization = *it;
        if (!parameterization.converter) log_fatal("Internal error: parameteriation has NULL converter");
        
        // all parameterizations could have the same converter,
        // but nevertheless, none of them must be initialized yet.
        if (parameterization.converter->IsInitialized())
            log_fatal("A parameterization converter is already initialized. Do not call their Initialize() method yourself!");
    }

    // now initialize them and set the medium properties and bias factors
    for (I3CLSimLightSourceParameterizationSeries::const_iterator it=parameterizations.begin();
         it!=parameterizations.end(); ++it)
    {
        const I3CLSimLightSourceParameterization &parameterization = *it;
        if (parameterization.converter->IsInitialized()) continue; // skip initialized converters
        
        parameterization.converter->SetRandomService(randomService_);
        parameterization.converter->SetMediumProperties(mediumProperties_);
        parameterization.converter->SetWlenBias(wlenBias_);
        parameterization.converter->SetBunchSizeGranularity(1); // we do not send the bunches directly, the steps are integrated in the step store first, so granularity does not matter
        parameterization.converter->SetMaxBunchSize(maxBunchSize_); // use the same bunch size for the parameterizations
        parameterization.converter->Initialize();
    }
    
    
    // making a copy of the medium properties
    {
        I3CLSimMediumPropertiesConstPtr copiedMediumProperties(new I3CLSimMediumProperties(*mediumProperties_));
        mediumProperties_ = copiedMediumProperties;
    }

    
    log_debug("Starting the Geant4 thread..");
    geant4Started_=false;
    barrier_is_enqueued_=false;

    geant4ThreadObj_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&I3CLSimLightSourceToStepConverterGeant4::Geant4Thread, this)));

    // wait for startup
    {
        boost::unique_lock<boost::mutex> guard(geant4Started_mutex_);
        for (;;)
        {
            if (geant4Started_) break;
            geant4Started_cond_.wait(guard);
        }
    }        
    
    log_debug("Geant4 thread started.");

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
            I3CLSimFunctionConstPtr refIndFunc = mProp->GetPhaseRefractiveIndex(i);
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

#ifdef HAS_GEANT4

#include "G4VStateDependent.hh"

class UserHookForAbortState : public G4VStateDependent
{
public:
    UserHookForAbortState() {;}
    virtual ~UserHookForAbortState() {;}
    
    virtual G4bool Notify(G4ApplicationState requiredState)
    {
        if (requiredState!=G4State_Abort) return true;

        // this seems to be caused by an G4Exception call.
        // throw a real exception here, so we are able to catch it.
        throw std::runtime_error("Geant4 is not amused.");

        return true;
    }
};

#endif


void I3CLSimLightSourceToStepConverterGeant4::Geant4Thread()
{
    // do not interrupt this thread by default
    boost::this_thread::disable_interruption di;

    try {
        Geant4Thread_impl(di);
    } catch(...) { // any exceptions?
        log_warn("Geant4 thread died unexpectedly..");
        LogGeant4Messages(true); // this is maybe the last chance to log..
        
        //throw; // don't bother cleaning up, we can't continue and the process is going to die anyway
    }
}

void I3CLSimLightSourceToStepConverterGeant4::Geant4Thread_impl(boost::this_thread::disable_interruption &di)
{
    // set up geant4
    
    // this thing stores all the steps generated by Geant4, sorted by the number
    // of Cherenkov photons they generate
    I3CLSimStepStorePtr stepStore(new I3CLSimStepStore( (isnan(maxNumPhotonsPerStep_)||(maxNumPhotonsPerStep_<0.))?0:(static_cast<uint32_t>(maxNumPhotonsPerStep_*1.5)) ));

    // this stores all particles that will be ent to parametrizations
    shared_ptr<std::deque<boost::tuple<I3CLSimLightSourceConstPtr, uint32_t, const I3CLSimLightSourceParameterization> > > sendToParameterizationQueue
    (new std::deque<boost::tuple<I3CLSimLightSourceConstPtr, uint32_t, const I3CLSimLightSourceParameterization> >());
    
#ifdef HAS_GEANT4
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
                                                   maxNumPhotonsPerStep_,
                                                   wlenBias_));
    //physics->RegisterPhysics(new TrkEMPhysicsUHE("EMUHE"));
    
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
                                                        sendToParameterizationQueue,
                                                        this->GetLightSourceParameterizationSeries(),
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
#endif
    
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
    const I3CLSimLightSourceParameterizationSeries parameterizations = this->GetLightSourceParameterizationSeries();

    // start the main loop
    for (;;)
    {
        //I3ParticleConstPtr particle;
        I3CLSimLightSourceConstPtr lightSource;
        uint32_t lightSourceIdentifier;

        {
            boost::this_thread::restore_interruption ri(di);
            try {
                ToGeant4Pair_t val = queueToGeant4_->Get();
                lightSourceIdentifier = val.first;
                lightSource = val.second;
            }
            catch(boost::thread_interrupted &i)
            {
                log_debug("G4 thread was interrupted. closing.");
                break;
            }
        }
        
        
        {
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
                        log_debug("G4 thread was interrupted. shutting down Geant4!");
                        
                        interruptionOccured = true;
                        break;
                    }
                }
            }            
            if (interruptionOccured) break;
        }
        
        
        if (!lightSource) {
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
                        log_debug("G4 thread was interrupted. closing.");
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
                        log_debug("G4 thread was interrupted. shutting down Geant4!");
                        
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
                    log_debug("G4 thread was interrupted. closing.");
                    break;
                }
            }
            
            //G4cout << " -> flush complete." << G4endl;

            // nothing to send to Geant4, so start from the beginning
            continue;
        }

        if (lightSource->GetType() == I3CLSimLightSource::Unknown)
        {
            log_warn("Ignoring a light source with type \"Unknown\".");
            continue;
        }
        
        if ((lightSource->GetType() == I3CLSimLightSource::Particle) && (lightSource->GetParticle().GetType() == I3Particle::unknown))
        {
            log_warn("Ignoring a particle with type \"unknown\".");
            continue;
        }

        // check if there is a parameterization for this particle, so we
        // may not even have to send it to Geant4
        bool parameterizationIsAvailable=false;

        // empty the queue        
        sendToParameterizationQueue->clear();

        for (I3CLSimLightSourceParameterizationSeries::const_iterator it=parameterizations.begin();
             it!=parameterizations.end(); ++it)
        {
            const I3CLSimLightSourceParameterization &parameterization = *it;
            
            if (parameterization.IsValidForLightSource(*lightSource))
            {
                sendToParameterizationQueue->push_back(boost::make_tuple(lightSource, lightSourceIdentifier, parameterization));
                parameterizationIsAvailable=true;
                break;
            }
        }
        
        
        if (!parameterizationIsAvailable) 
        {
#ifdef HAS_GEANT4
            // no parameterization was found, use default Geant4
            if (lightSource->GetType() != I3CLSimLightSource::Particle)
            {
                G4cerr << "Geant4 can only handle light sources of type \"Particle\". Please add parameterizations for light source of other types (e.g. flashers). Ignoring non-particle light source." << G4endl;

                continue;
            }
            
            const I3Particle &particle = lightSource->GetParticle();

            // configure the Geant4 particle gun
            {
                G4ParticleGun *particleGun = thePrimaryGenerator->GetParticleGun();
                if (!particleGun) log_fatal("Internal error: G4ParticleGun instance is NULL!");
                
                const bool ret = I3CLSimI3ParticleGeantConverter::SetParticleGun(particleGun, particle);
                
                if (!ret) {
                    G4cerr << "Could not configure Geant4 to shoot a " << particle.GetTypeString() << "! Ignoring." << G4endl;

                    continue;
                }
            }
            
            // set the current particle ID
            theEventAction->SetExternalParticleID(lightSourceIdentifier);

            log_trace("Geant4: no parameterization for %s with E=%fGeV", particle.GetTypeString().c_str(), particle.GetEnergy()/I3Units::GeV);
            
            G4cout << "Geant4: shooting a " << particle.GetTypeString() << " with id " << lightSourceIdentifier << " and E=" << particle.GetEnergy()/I3Units::GeV << "GeV." << G4endl;

            // we don't need the light source anymore.
            lightSource.reset();
            
            // turn on the Geant4 beam!
            // (this fills the stepStore with steps and our particle list with
            // output particles for the available parameterizations..)
            runManager->BeamOn(1);
            
            // check if AbortRun was requested beacause of a thread interruption.
            if (theEventAction->AbortWasRequested()) break;
#else
            log_fatal("Geant4 is not available, but there is no parameterization for the current particle.");
#endif
        }
        else
        {
            // we don't need the light source anymore.
            lightSource.reset();
        }
        

        // loop over all particles that should be sent to a parameterization and send them.
        // They were added by Geant4 (our by this code directly) to sendToParameterizationQueue
        bool interruptionOccured=false;
        typedef boost::tuple<I3CLSimLightSourceConstPtr, uint32_t, const I3CLSimLightSourceParameterization> particleAndIndexAndParamTuple_t;
        BOOST_FOREACH(const particleAndIndexAndParamTuple_t &particleAndIndexPair, *sendToParameterizationQueue)
        {
            lightSource = particleAndIndexPair.get<0>();              // re-use the lightSource and lightSourceIdentifier variables (the original
            lightSourceIdentifier = particleAndIndexPair.get<1>();    // ones are not needed anymore)
            const I3CLSimLightSourceParameterization &parameterization = particleAndIndexPair.get<2>();
            
            if (!parameterization.IsValidForLightSource(*lightSource))
                log_fatal("internal error. Parameterization in queue is not valid for the light source that came with it..");

            // call the converter
            if (!parameterization.converter) log_fatal("Internal error: parameteriation has NULL converter");
            if (!parameterization.converter->IsInitialized()) log_fatal("Internal error: parameterization converter is not initialized.");
            if (parameterization.converter->BarrierActive()) log_fatal("Logic error: parameterization converter has active barrier.");
                    
            parameterization.converter->EnqueueLightSource(*lightSource, lightSourceIdentifier);
            parameterization.converter->EnqueueBarrier();
            
            // get steps from the parameterization until the barrier is reached
            for (;;)
            {
                I3CLSimStepSeriesConstPtr res;
                bool barrierHasBeenReached=false;
                
                {
                    boost::this_thread::restore_interruption ri(di);
                    try {
                        // this blocks if there are no steps yet and the
                        // parameterization code is still working.
                        res = parameterization.converter->GetConversionResultWithBarrierInfo(barrierHasBeenReached);
                    } catch(boost::thread_interrupted &i) {
                        log_debug("G4 thread was interrupted. shutting down Geant4!");
                        interruptionOccured = true;
                        break;
                    }
                }

                
                if (!res) {
                    log_debug("NULL result from parameterization GetConversionResult(). ignoring.");
                } else {
                    // add steps from the parameterization to the step store
                    BOOST_FOREACH(const I3CLSimStep &step, *res)
                    {
                        stepStore->insert_copy(step.GetNumPhotons(), step);
                    }
                }
                
                // we don't need the results anymore
                res.reset();
                
                // push steps out if there are enough of them
                while (stepStore->size() >= maxBunchSize_)
                {
                    I3CLSimStepSeriesPtr steps(new I3CLSimStepSeries());
                    stepStore->pop_bunch_to_vector(maxBunchSize_, *steps);
                    
                    {
                        boost::this_thread::restore_interruption ri(di);
                        try {
                            // this blocks if the queue from Geant4 to
                            // OpenCL is full.
                            queueFromGeant4_->Put(std::make_pair(steps, false));
                        } catch(boost::thread_interrupted &i) {
                            log_debug("G4 thread was interrupted. shutting down Geant4!");
                            interruptionOccured = true;
                            break;
                        }
                    }
                }            
                if (interruptionOccured) break;

                if (barrierHasBeenReached) break; // get out of the loop if the barrier has been reached
            }

            if (interruptionOccured) break;
        }

        if (interruptionOccured) break;

        // empty the queue and clean up  
        sendToParameterizationQueue->clear();
        lightSource.reset();
        
    }

#ifdef HAS_GEANT4
    G4cout << "G4 thread terminating..." << G4endl;

    UI->SetCoutDestination(NULL);
    delete theUISessionToQueue;
    
    // job termination
    delete runManager;
#endif

    log_debug("G4 thread terminated.");
}

bool I3CLSimLightSourceToStepConverterGeant4::IsInitialized() const
{
    LogGeant4Messages();

    return initialized_;
}

void I3CLSimLightSourceToStepConverterGeant4::SetBunchSizeGranularity(uint64_t num)
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 already initialized!");
    
    if (num<=0)
        throw I3CLSimLightSourceToStepConverter_exception("BunchSizeGranularity of 0 is invalid!");
    
    bunchSizeGranularity_=num;
}

void I3CLSimLightSourceToStepConverterGeant4::SetMaxBunchSize(uint64_t num)
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 already initialized!");

    if (num<=0)
        throw I3CLSimLightSourceToStepConverter_exception("MaxBunchSize of 0 is invalid!");

    maxBunchSize_=num;
}

void I3CLSimLightSourceToStepConverterGeant4::SetRandomService(I3RandomServicePtr random)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 already initialized!");
    
    randomService_=random;
    
    // TODO: eventually Geant4 should use the IceTray rng!!
    randomSeed_ = randomService_->Integer(900000000);
}

void I3CLSimLightSourceToStepConverterGeant4::SetWlenBias(I3CLSimFunctionConstPtr wlenBias)
{
    LogGeant4Messages();
    
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 already initialized!");
    
    wlenBias_=wlenBias;
}

void I3CLSimLightSourceToStepConverterGeant4::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 already initialized!");

    mediumProperties_=mediumProperties;
}

void I3CLSimLightSourceToStepConverterGeant4::EnqueueLightSource(const I3CLSimLightSource &lightSource, uint32_t identifier)
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
        if (barrier_is_enqueued_)
            throw I3CLSimLightSourceToStepConverter_exception("A barrier is enqueued! You must receive all steps before enqueuing a new particle.");
    }
    
    I3CLSimLightSourceConstPtr lightSourceCopy(new I3CLSimLightSource(lightSource));
    queueToGeant4_->Put(std::make_pair(identifier, lightSourceCopy));
    
    LogGeant4Messages();
}

void I3CLSimLightSourceToStepConverterGeant4::EnqueueBarrier()
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
        if (barrier_is_enqueued_)
            throw I3CLSimLightSourceToStepConverter_exception("A barrier is already enqueued!");
        
        barrier_is_enqueued_=true;

        // we use a NULL pointer as the barrier
        queueToGeant4_->Put(std::make_pair(0, I3CLSimLightSourceConstPtr()));
    }
    
    LogGeant4Messages();
}

bool I3CLSimLightSourceToStepConverterGeant4::BarrierActive() const
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);

        return barrier_is_enqueued_;
    }
}

bool I3CLSimLightSourceToStepConverterGeant4::MoreStepsAvailable() const
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 is not initialized!");

    return (!queueFromGeant4_->empty());
}

I3CLSimStepSeriesConstPtr I3CLSimLightSourceToStepConverterGeant4::GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
{
    LogGeant4Messages();

    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterGeant4 is not initialized!");

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

void I3CLSimLightSourceToStepConverterGeant4::LogGeant4Messages(bool allAsWarn) const
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
                    log_debug("Geant4 says:  %s", out.c_str());
            }
        }
    }
}

