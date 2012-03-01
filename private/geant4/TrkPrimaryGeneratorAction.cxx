#include "TrkPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

TrkPrimaryGeneratorAction::TrkPrimaryGeneratorAction()
{
    particleGun = new G4ParticleGun();
}

TrkPrimaryGeneratorAction::~TrkPrimaryGeneratorAction()
{
    delete particleGun;
}

void TrkPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    particleGun->GeneratePrimaryVertex(anEvent);
}


