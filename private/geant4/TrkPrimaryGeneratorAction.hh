#ifndef TrkPrimaryGeneratorAction_h
#define TrkPrimaryGeneratorAction_h 1

#include "G4GeneralParticleSource.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class TrkPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    TrkPrimaryGeneratorAction();
    ~TrkPrimaryGeneratorAction();

public:
    void GeneratePrimaries(G4Event* anEvent);

    inline G4ParticleGun *GetParticleGun() {return particleGun;}

private:
    G4ParticleGun* particleGun;
};

#endif


