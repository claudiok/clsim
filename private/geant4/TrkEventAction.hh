#ifndef TrkEventAction_h
#define TrkEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "clsim/I3CLSimStepStore.h"
#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimQueue.h"
#include "clsim/I3CLSimParticleToStepConverterGeant4.h"

#include <boost/thread.hpp>

class G4Event;

class TrkEventAction : public G4UserEventAction
{
public:
    TrkEventAction(uint64_t maxBunchSize,
                   I3CLSimStepStorePtr stepStore,
                   double electronPositronMinEnergyForSecondary,
                   double electronPositronMaxEnergyForSecondary,
                   boost::shared_ptr<I3CLSimQueue<I3CLSimParticleToStepConverterGeant4::FromGeant4Pair_t> > queueFromGeant4,
                   boost::this_thread::disable_interruption &threadDisabledInterruptionState,
                   double maxRefractiveIndex);
    ~TrkEventAction();
    
public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    
    inline bool AbortWasRequested() {return abortRequested_;}
    inline void SetExternalParticleID(uint32_t val) {currentExternalParticleID_=val;}
    
private:
    bool abortRequested_;
    uint64_t maxBunchSize_;
    I3CLSimStepStorePtr stepStore_;
    uint32_t currentExternalParticleID_;
    
    double electronPositronMinEnergyForSecondary_;
    double electronPositronMaxEnergyForSecondary_;
    
    boost::shared_ptr<I3CLSimQueue<I3CLSimParticleToStepConverterGeant4::FromGeant4Pair_t> > queueFromGeant4_;
    boost::this_thread::disable_interruption &threadDisabledInterruptionState_;
    
    double maxRefractiveIndex_;
};

#endif
