/*
 *  TrkOpticalPhysics.cc
 *  ShowerSim
 *
 *  Created by Claudio Kopper on 27.09.2006.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "TrkEMPhysicsUHE.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>  

#include "G4ProcessManager.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4OpticalPhoton.hh"

#include "G4LossTableManager.hh"

TrkEMPhysicsUHE::TrkEMPhysicsUHE(const G4String& name)
  :  G4VPhysicsConstructor(name)
{
	G4cout << "<<<< UHE EM Physics and Loss Table Extension" << G4endl;
}

TrkEMPhysicsUHE::~TrkEMPhysicsUHE()
{
}


void TrkEMPhysicsUHE::ConstructParticle()
{
	// no additional particles
}

void TrkEMPhysicsUHE::ConstructProcess()
{
	// Construct Processes
	
	//extend binning of PhysicsTables
	G4LossTableManager::Instance()->SetMaxEnergy(10.*PeV);
	G4LossTableManager::Instance()->SetDEDXBinning(440);
	G4LossTableManager::Instance()->SetLambdaBinning(440);
	G4LossTableManager::Instance()->SetVerbose(0);
}
