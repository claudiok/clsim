/*
 *  TrkOpticalPhysics.cc
 *  ShowerSim
 *
 *  Created by Claudio Kopper on 27.09.2006.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "TrkOpticalPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>  

#include "G4ProcessManager.hh"

#include "TrkCerenkov.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4OpticalPhoton.hh"

TrkOpticalPhysics::TrkOpticalPhysics(const G4String& name, 
                                     double maxBetaChangePerStep,
                                     uint32_t maxNumPhotonsPerStep,
                                     I3CLSimFunctionConstPtr wlenBias
                                     )
:  G4VPhysicsConstructor(name),
maxBetaChangePerStep_(maxBetaChangePerStep),
maxNumPhotonsPerStep_(maxNumPhotonsPerStep),
wlenBias_(wlenBias)
{
    G4cout << "<<<< Optical Processes (TrkCerenkov)" << G4endl;
}

TrkOpticalPhysics::~TrkOpticalPhysics()
{
}


void TrkOpticalPhysics::ConstructParticle()
{
}

void TrkOpticalPhysics::ConstructProcess()
{
    // Construct Processes
    
    theCerenkovProcess=new TrkCerenkov();

    theCerenkovProcess->SetMaxBetaChangePerStep(maxBetaChangePerStep_);
    theCerenkovProcess->SetMaxNumPhotonsPerStep(maxNumPhotonsPerStep_);
    theCerenkovProcess->SetWlenBiasFunction(wlenBias_);
    
    // Add the processes to their respective particles
    G4ProcessManager * pManager = NULL;
    
    
    theParticleIterator->reset();
    while ((*theParticleIterator)())
    {
        G4ParticleDefinition* particle = theParticleIterator->value();
        pManager = particle->GetProcessManager();
        if (theCerenkovProcess->IsApplicable(*particle))
        {
            pManager->AddProcess(theCerenkovProcess);
            pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
        }
    }
    
}
