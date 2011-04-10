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

#include "I3CLSimI3ParticleGeantConverter.hh"

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

    // see if there are eny parameterizations available for this particle
    const I3CLSimParticleParameterizationSeries &parameterizations = eventInformation->parameterizationAvailable;

    const I3Particle::ParticleType trackI3ParticleType =
    I3CLSimI3ParticleGeantConverter::ConvertPDGEncodingToI3ParticleType(aTrack->GetDefinition()->GetPDGEncoding());
    
    const G4double trackEnergy = aTrack->GetKineticEnergy();

    if (trackI3ParticleType==I3Particle::unknown)
    {
        // there are no parameterizations for particles unknown to IceTray
        return fUrgent;
    }
    
    for (I3CLSimParticleParameterizationSeries::const_iterator it=parameterizations.begin();
         it!=parameterizations.end(); ++it)
    {
        const I3CLSimParticleParameterization &parameterization = *it;
        
        if (parameterization.IsValid(trackI3ParticleType, trackEnergy*I3Units::GeV/GeV))
        {
            I3CLSimStepStorePtr stepStore = eventInformation->stepStore;

            I3ParticlePtr particle(new I3Particle());
            
            const G4ThreeVector &trackPos = aTrack->GetPosition();
            const G4double trackTime = aTrack->GetGlobalTime();
            const G4ThreeVector &trackDir = aTrack->GetMomentumDirection();

            particle->SetType(trackI3ParticleType);
            particle->SetPos(trackPos.x()*I3Units::m/m,trackPos.y()*I3Units::m/m,trackPos.z()*I3Units::m/m);
            particle->SetDir(trackDir.x(),trackDir.y(),trackDir.z());
            particle->SetTime(trackTime*I3Units::ns/ns);
            particle->SetEnergy(trackEnergy*I3Units::GeV/GeV);

            //G4cout << "Geant4: sending a " << particle->GetTypeString() << " with id " << eventInformation->currentExternalParticleID << " and E=" << particle->GetEnergy()/I3Units::GeV << "GeV to a parameterization handler." << G4endl;
            
            // call the converter
            if (!parameterization.converter) log_fatal("Internal error: parameteriation has NULL converter");
            if (!parameterization.converter->IsInitialized()) log_fatal("Internal error: parameterization converter is not initialized.");
            if (parameterization.converter->BarrierActive()) log_fatal("Logic error: parameterization converter has active barrier.");
            
            parameterization.converter->EnqueueParticle(*particle, eventInformation->currentExternalParticleID);
            parameterization.converter->EnqueueBarrier();
            
            while (parameterization.converter->BarrierActive())
            {
                I3CLSimStepSeriesConstPtr res =
                parameterization.converter->GetConversionResult();
                
                if (!res) {
                    log_warn("NULL result from parameterization GetConversionResult(). ignoring.");
                    continue;
                }
                if (res->size()==0) continue; // ignore empty vectors
                
                BOOST_FOREACH(const I3CLSimStep &step, *res)
                {
                    stepStore->insert_copy(step.GetNumPhotons(), step);
                }
            }
            
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

