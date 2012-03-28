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

    // see if there are eny parameterizations available for this particle
    const I3CLSimLightSourceParameterizationSeries &parameterizations = eventInformation->parameterizationAvailable;

#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
    const I3Particle::ParticleType trackI3ParticleType =
    I3CLSimI3ParticleGeantConverter::ConvertPDGEncodingToI3ParticleType(aTrack->GetDefinition()->GetPDGEncoding());
#endif
    
    const G4double trackEnergy = aTrack->GetKineticEnergy();

#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
    if (trackI3ParticleType==I3Particle::unknown)
    {
        // there are no parameterizations for particles unknown to IceTray
        return fUrgent;
    }
#endif
    
    for (I3CLSimLightSourceParameterizationSeries::const_iterator it=parameterizations.begin();
         it!=parameterizations.end(); ++it)
    {
        const I3CLSimLightSourceParameterization &parameterization = *it;
        
#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        if (parameterization.IsValid(trackI3ParticleType, trackEnergy*I3Units::GeV/GeV))
#else
        if (parameterization.IsValidForPdgEncoding(aTrack->GetDefinition()->GetPDGEncoding(), trackEnergy*I3Units::GeV/GeV))
#endif
        {
            shared_ptr<std::deque<boost::tuple<I3CLSimLightSourceConstPtr, uint32_t, const I3CLSimLightSourceParameterization> > > sendToParameterizationQueue = eventInformation->sendToParameterizationQueue;

            if (!sendToParameterizationQueue) 
                log_fatal("internal error: sendToParameterizationQueue==NULL");
            
            I3Particle particle;
            
            const G4ThreeVector &trackPos = aTrack->GetPosition();
            const G4double trackTime = aTrack->GetGlobalTime();
            const G4ThreeVector &trackDir = aTrack->GetMomentumDirection();

#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
            particle.SetPdgEncoding(aTrack->GetDefinition()->GetPDGEncoding());
#else
            particle.SetType(trackI3ParticleType);
#endif
            particle.SetPos(trackPos.x()*I3Units::m/m,trackPos.y()*I3Units::m/m,trackPos.z()*I3Units::m/m);
            particle.SetDir(trackDir.x(),trackDir.y(),trackDir.z());
            particle.SetTime(trackTime*I3Units::ns/ns);
            particle.SetEnergy(trackEnergy*I3Units::GeV/GeV);

            I3CLSimLightSourcePtr lightSource(new I3CLSimLightSource(particle));
            sendToParameterizationQueue->push_back(boost::make_tuple(lightSource, eventInformation->currentExternalParticleID, parameterization));

            return fKill;
        }
    }
    
        
    return fUrgent;
}

void TrkStackingAction::NewStage()
{
}

void TrkStackingAction::PrepareNewEvent()
{
}

