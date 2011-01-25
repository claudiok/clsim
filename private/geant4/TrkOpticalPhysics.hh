#ifndef TrkOpticalPhysics_hh
#define TrkOpticalPhysics_hh

#include "G4VPhysicsConstructor.hh"

class TrkCerenkov;

class TrkOpticalPhysics : public G4VPhysicsConstructor
{
public:
    TrkOpticalPhysics(const G4String& name,
                      double maxBetaChangePerStep,
                      double maxNumPhotonsPerStep);
    virtual ~TrkOpticalPhysics();
    
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    
protected:
    TrkCerenkov* theCerenkovProcess;
    
private:
    double maxBetaChangePerStep_;
    double maxNumPhotonsPerStep_;
};

#endif
