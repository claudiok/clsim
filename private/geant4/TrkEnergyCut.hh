#ifndef TrkEnergyCut_hh
#define TrkEnergyCut_hh

#include "G4VPhysicsConstructor.hh"

class G4UserSpecialCuts;

class TrkEnergyCut : public G4VPhysicsConstructor
{
public:
    TrkEnergyCut(const G4String& name);
    virtual ~TrkEnergyCut();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

protected:
    G4UserSpecialCuts *theSpecialCutProcess;
};

#endif
