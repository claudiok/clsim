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
 * @file I3CLSimLightSourcePropagatorGeant4.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "clsim/I3CLSimLightSourceToStepConverter.h"
#include "geant4/I3CLSimLightSourcePropagatorGeant4.h"

#include <limits>
#include <deque>
#include <boost/tuple/tuple.hpp>
#include <cmath>

// geant4 stuff
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4String.hh"

#include "G4PhysListFactory.hh"
#include "G4ParticleGun.hh"

#include "G4Version.hh"

#include "TrkOpticalPhysics.hh"

#include "TrkDetectorConstruction.hh"
#include "TrkPrimaryGeneratorAction.hh"
#include "TrkEventAction.hh"
#include "TrkStackingAction.hh"
//#include "TrkSteppingAction.hh"
#include "TrkUISessionToQueue.hh"

#include "I3CLSimI3ParticleGeantConverter.hh"
#include "Randomize.hh"

#include "dataclasses/physics/I3Particle.h"
#include "clsim/I3CLSimLightSource.h"
#include "clsim/function/I3CLSimFunction.h"
#include "clsim/I3CLSimMediumProperties.h"

#if G4VERSION_NUMBER >= 1000
using CLHEP::mm;
#endif

// other headers
#include <stdlib.h>


// static definitions

std::atomic<bool> I3CLSimLightSourcePropagatorGeant4::thereCanBeOnlyOneGeant4(false);

const std::string I3CLSimLightSourcePropagatorGeant4::default_physicsListName="QGSP_BERT_EMV";
const double I3CLSimLightSourcePropagatorGeant4::default_maxBetaChangePerStep=10.*I3Units::perCent;
const uint32_t I3CLSimLightSourcePropagatorGeant4::default_maxNumPhotonsPerStep=200;

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

I3CLSimLightSourcePropagatorGeant4::I3CLSimLightSourcePropagatorGeant4(std::string physicsListName,
                                                                           double maxBetaChangePerStep,
                                                                           uint32_t maxNumPhotonsPerStep
                                                                           )
:
queueFromGeant4Messages_(new I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > >(UINT_MAX)), // no maximum size
physicsListName_(physicsListName),
maxBetaChangePerStep_(maxBetaChangePerStep),
maxNumPhotonsPerStep_(maxNumPhotonsPerStep),
initialized_(false)
{
    // Geant4 keeps LOTS of global state and is inherently non-thread-safe.
    // To prevent users from using more than one instance of Geant4 within the
    // same process, we throw an exception in case another instance
    // of this class already exists.
    
    {
        if (thereCanBeOnlyOneGeant4) 
            throw I3CLSimLightSourceToStepConverter_exception("There can be only one! ...instance of I3CLSimLightSourcePropagatorGeant4.");
        
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
    
    // get the pointer to the UI manager and set our custom G4cout/G4cerr destination
    G4UImanager* UI = G4UImanager::GetUIpointer();
    TrkUISessionToQueue *theUISessionToQueue = new TrkUISessionToQueue(queueFromGeant4Messages_);
    UI->SetCoutDestination(theUISessionToQueue);
    
    // keep this around to "catch" G4Eceptions and throw real exceptions
    UserHookForAbortState *theUserHookForAbortState = new UserHookForAbortState();
    G4StateManager *theStateManager = G4StateManager::GetStateManager();
    theStateManager->RegisterDependent(theUserHookForAbortState);
}

I3CLSimLightSourcePropagatorGeant4::~I3CLSimLightSourcePropagatorGeant4()
{
    LogGeant4Messages();
    
    {
        thereCanBeOnlyOneGeant4=false;
    }

    LogGeant4Messages();
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
                
                if (std::isnan(maxVal) || (rindex > maxVal)) maxVal=rindex;
            }
            
        }
        
        return maxVal;
    }
}

void I3CLSimLightSourcePropagatorGeant4::Initialize()
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourcePropagatorGeant4 already initialized!");

    if (!randomService_)
        throw I3CLSimLightSourceToStepConverter_exception("RandomService not set!");

    if (!wlenBias_)
        throw I3CLSimLightSourceToStepConverter_exception("WlenBias not set!");

    if (!mediumProperties_)
        throw I3CLSimLightSourceToStepConverter_exception("MediumProperties not set!");
    
    // making a copy of the medium properties
    {
        I3CLSimMediumPropertiesConstPtr copiedMediumProperties(new I3CLSimMediumProperties(*mediumProperties_));
        mediumProperties_ = copiedMediumProperties;
    }

    LogGeant4Messages();

    // initialize the run manager
    runManager_.reset(new G4RunManager);
    
    // set up the "detector" (a lot of water)
    runManager_->SetUserInitialization(new TrkDetectorConstruction(mediumProperties_));
    
    // set up the physics list (something+Optical Physics)
    std::unique_ptr<G4PhysListFactory> factory(new G4PhysListFactory);
    G4VModularPhysicsList *physics = factory->GetReferencePhysList(physicsListName_.c_str());
    
    physics->RegisterPhysics(new TrkOpticalPhysics("Optical",
                                                   maxBetaChangePerStep_,
                                                   maxNumPhotonsPerStep_,
                                                   wlenBias_));
    
    physics->SetDefaultCutValue(0.25*mm);
    runManager_->SetUserInitialization(physics);
    
    // instantiate some Geant4 helper classes
    TrkPrimaryGeneratorAction *thePrimaryGenerator = new TrkPrimaryGeneratorAction();
    
    runManager_->SetUserAction(thePrimaryGenerator);     // runManager now owns this pointer
    runManager_->SetUserAction(new TrkStackingAction());   // runManager now owns this pointer

    const double maxRefractiveIndex = GetMaxRIndex(mediumProperties_);
    G4cout << "overall maximum refractive index is " << maxRefractiveIndex << G4endl;
    if (std::isnan(maxRefractiveIndex)) log_fatal("No maximum refractive index could be found");
    
    TrkEventAction *theEventAction = new TrkEventAction(maxRefractiveIndex);
    runManager_->SetUserAction(theEventAction);      // runManager now owns this pointer
    
    G4UImanager* UI = G4UImanager::GetUIpointer();
    //UI->ApplyCommand("/physics_engine/tailor/SyncRadiation on");
    //UI->ApplyCommand("/physics_engine/tailor/GammaNuclear on");
    UI->ApplyCommand("/physics_engine/tailor/MuonNuclear on");
    runManager_->Initialize();

    //UI->ApplyCommand("/run/particle/dumpCutValues");
    
    CLHEP::HepRandom::setTheSeed(randomSeed_); // value [0,900000000]

    initialized_=true;
}

void I3CLSimLightSourcePropagatorGeant4::Convert(I3CLSimLightSourceConstPtr &lightSource, uint32_t lightSourceIdentifier,
    secondary_callback emitSecondary, step_callback emitStep)
{
    const I3Particle &particle = lightSource->GetParticle();

    // configure the Geant4 particle gun
    {
        TrkPrimaryGeneratorAction *thePrimaryGenerator =
            const_cast<TrkPrimaryGeneratorAction*>(dynamic_cast<const TrkPrimaryGeneratorAction*>(
            runManager_->GetUserPrimaryGeneratorAction()));
        G4ParticleGun *particleGun = thePrimaryGenerator->GetParticleGun();
        if (!particleGun) log_fatal("Internal error: G4ParticleGun instance is NULL!");
        
        const bool ret = I3CLSimI3ParticleGeantConverter::SetParticleGun(particleGun, particle);
        
        if (!ret) {
            G4cerr << "Could not configure Geant4 to shoot a " << particle.GetTypeString() << "! Ignoring." << G4endl;

            return;
        }
    }
    
    // set the current particle ID
    TrkEventAction *theEventAction = const_cast<TrkEventAction*>(dynamic_cast<const TrkEventAction*>(runManager_->GetUserEventAction()));
    theEventAction->SetExternalParticleID(lightSourceIdentifier);
    theEventAction->SetSecondaryCallback(emitSecondary);
    theEventAction->SetStepCallback(emitStep);
    
    G4cout << "Geant4: shooting a " << particle.GetTypeString() << " with id " << lightSourceIdentifier << " and E=" << particle.GetEnergy()/I3Units::GeV << "GeV." << G4endl;
    
    // turn on the Geant4 beam!
    // (this fills the stepStore with steps and our particle list with
    // output particles for the available parameterizations..)
    runManager_->BeamOn(1);
}

bool I3CLSimLightSourcePropagatorGeant4::IsInitialized() const
{
    LogGeant4Messages();

    return initialized_;
}

void I3CLSimLightSourcePropagatorGeant4::SetRandomService(I3RandomServicePtr random)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourcePropagatorGeant4 already initialized!");
    
    randomService_=random;
    
    // TODO: eventually Geant4 should use the IceTray rng!!
    randomSeed_ = randomService_->Integer(900000000);
}

void I3CLSimLightSourcePropagatorGeant4::SetWlenBias(I3CLSimFunctionConstPtr wlenBias)
{
    LogGeant4Messages();
    
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourcePropagatorGeant4 already initialized!");
    
    wlenBias_=wlenBias;
}

void I3CLSimLightSourcePropagatorGeant4::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    LogGeant4Messages();

    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourcePropagatorGeant4 already initialized!");

    mediumProperties_=mediumProperties;
}

void I3CLSimLightSourcePropagatorGeant4::LogGeant4Messages(bool allAsWarn) const
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

