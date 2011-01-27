#include "TrkStackingAction.hh"
#include "TrkUserEventInformation.hh"

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

#include "G4UnitsTable.hh"

TrkStackingAction::TrkStackingAction()
{
}

TrkStackingAction::~TrkStackingAction()
{
}

G4ClassificationOfNewTrack TrkStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
	TrkUserEventInformation* eventInformation =
	(TrkUserEventInformation*)G4EventManager::GetEventManager()
	->GetConstCurrentEvent()->GetUserInformation();

	{
        const double maxRefractiveIndex = eventInformation->maxRefractiveIndex;
        const double BetaInverse = c_light/aTrack->GetVelocity();
        
        // below the Cherekov threshold?
        if (BetaInverse > maxRefractiveIndex)
            return fKill;
    }
    
    /*
    // is it the primary particle?
    if (aTrack->GetParentID()!=0)
    {
        // it's NOT the primary! check if this new particle is inside the (extended) can
        
        const double canHeight = 1000.*m;
        const double canRadius = 750.*m;
        
        G4ThreeVector particlePosRelToCan = aTrack->GetPosition(); // - SSimDetectorConstruction::canPosition;
        
        if (fabs(particlePosRelToCan.z()) > canHeight/2.) {
            // we are above or below the can
            return fKill;
        } else {
            G4double posRadius = std::sqrt(particlePosRelToCan.x()*particlePosRelToCan.x() + particlePosRelToCan.y()*particlePosRelToCan.y());
            if (posRadius > canRadius) {
                return fKill;
            }
        }
        
        //return fUrgent;
    }
    */
    
	if((aTrack->GetDefinition()==G4Electron::ElectronDefinition()) ||
	   (aTrack->GetDefinition()==G4Positron::PositronDefinition()))
	{ 
		G4double E = aTrack->GetKineticEnergy();

		const double lowCutoff = eventInformation->electronPositronMinEnergyForSecondary;
		const double highCutoff = eventInformation->electronPositronMaxEnergyForSecondary;
		
		if ((E >= lowCutoff) && (E <= highCutoff)) 
        {
            const G4ThreeVector &trackPos = aTrack->GetPosition();
            const G4double trackTime = aTrack->GetGlobalTime();
            const G4double trackEnergy = aTrack->GetKineticEnergy();
            const G4ThreeVector &trackDir = aTrack->GetMomentumDirection();
            
            I3ParticlePtr particle(new I3Particle());
            
            if (aTrack->GetDefinition()==G4Electron::ElectronDefinition())
                particle->SetType(I3Particle::EMinus);
            if (aTrack->GetDefinition()==G4Positron::PositronDefinition())
                particle->SetType(I3Particle::EPlus);
            
            particle->SetPos(trackPos.x()*I3Units::m/m,trackPos.y()*I3Units::m/m,trackPos.z()*I3Units::m/m);
            particle->SetDir(trackDir.x(),trackDir.y(),trackDir.z());
            particle->SetTime(trackTime*I3Units::ns/ns);
            particle->SetEnergy(trackEnergy*I3Units::GeV/GeV);
         
            std::pair<uint32_t, I3ParticleConstPtr> particleToHost(eventInformation->currentExternalParticleID, particle);

            particle.reset();
            
            {
                boost::this_thread::restore_interruption ri(eventInformation->threadDisabledInterruptionState);
                try {
                    eventInformation->queueFromGeant4->Put(std::make_pair(particleToHost, false));
                } catch(boost::thread_interrupted &i) {
                    G4cout << "G4 thread was interrupted. shutting down Geant4!" << G4endl;
                    
                    eventInformation->abortRequested = true;
                    G4RunManager::GetRunManager()->AbortRun();
                }
            }
            
            //G4cout << "Geant4 just sent a particle of " << E/MeV << "MeV." << G4endl;
            
            return fKill;
		}
	}
    
	return fUrgent;
}

void TrkStackingAction::NewStage()
{
	//G4cout << "New stage!" << G4endl;
}

void TrkStackingAction::PrepareNewEvent()
{
	//G4cout << "Prepare new event!" << G4endl;
}

