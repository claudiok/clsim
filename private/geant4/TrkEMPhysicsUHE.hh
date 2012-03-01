#ifndef TrkEMPhysicsUHE_hh
#define TrkEMPhysicsUHE_hh

#include "G4VPhysicsConstructor.hh"

class TrkEMPhysicsUHE : public G4VPhysicsConstructor
{
public:
    TrkEMPhysicsUHE(const G4String& name);
    virtual ~TrkEMPhysicsUHE();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

protected:
};

#endif
