#include "TrkEventAction.hh"
#include "TrkUserEventInformation.hh"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"

#include "G4ParticleTable.hh"

TrkEventAction::TrkEventAction(uint64_t maxBunchSize,
                               I3CLSimStepStorePtr stepStore,
                               const I3CLSimParticleParameterizationSeries &parameterizationAvailable,
                               boost::shared_ptr<I3CLSimQueue<I3CLSimParticleToStepConverterGeant4::FromGeant4Pair_t> > queueFromGeant4,
                               boost::this_thread::disable_interruption &threadDisabledInterruptionState,
                               double maxRefractiveIndex)
:
abortRequested_(false),
maxBunchSize_(maxBunchSize),
stepStore_(stepStore),
parameterizationAvailable_(parameterizationAvailable),
queueFromGeant4_(queueFromGeant4),
threadDisabledInterruptionState_(threadDisabledInterruptionState),
maxRefractiveIndex_(maxRefractiveIndex)
{
}

TrkEventAction::~TrkEventAction()
{
}

void TrkEventAction::BeginOfEventAction(const G4Event* anEvent)
{
	// New event, add the user information object
	TrkUserEventInformation* eventInformation = 
    new TrkUserEventInformation(maxBunchSize_,
                                stepStore_,
                                parameterizationAvailable_,
                                queueFromGeant4_,
                                threadDisabledInterruptionState_,
                                currentExternalParticleID_,
                                maxRefractiveIndex_);

	G4EventManager::GetEventManager()->SetUserInformation(eventInformation);
	
    eventInformation->StartClock();
}

void TrkEventAction::EndOfEventAction(const G4Event* anEvent)
{
	//G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
	
	TrkUserEventInformation* eventInformation =
	(TrkUserEventInformation*)anEvent->GetUserInformation();

    abortRequested_ = eventInformation->abortRequested;
}
