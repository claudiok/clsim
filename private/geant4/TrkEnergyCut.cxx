/*
 *  TrkOpticalPhysics.cc
 *  ShowerSim
 *
 *  Created by Claudio Kopper on 27.09.2006.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "TrkEnergyCut.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>  

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4UserSpecialCuts.hh"

TrkEnergyCut::TrkEnergyCut(const G4String& name)
  :  G4VPhysicsConstructor(name)
{
    G4cout << "<<<< Energy cut Process (kills particles below an energy threshold)" << G4endl;
}

TrkEnergyCut::~TrkEnergyCut()
{
}


void TrkEnergyCut::ConstructParticle()
{
    // no particles to construct
}

void TrkEnergyCut::ConstructProcess()
{
    // Construct Processes
    theSpecialCutProcess = new G4UserSpecialCuts();
    
    // Add the processes to their respective particles
    G4ProcessManager * pManager = 0;
    
    theParticleIterator->reset();
    while ((*theParticleIterator)())
    {
        G4ParticleDefinition* particle = theParticleIterator->value();
        pManager = particle->GetProcessManager();
        if (theSpecialCutProcess->IsApplicable(*particle))
        {
            pManager->AddProcess(theSpecialCutProcess,-1,-1,9000);
        }
    }

}
