/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file TrkStackingAction.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

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
#include "G4Version.hh"

#if G4VERSION_NUMBER >= 1000
using CLHEP::m;
using CLHEP::ns;
using CLHEP::GeV;
#endif

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

    // changed Apr 19, 2012: do not automatically assume a particle
    // below the Cherenkov threshold should be killed. It might decay
    // later while we can still se it.
    //
    // {
    //     const double maxRefractiveIndex = eventInformation->maxRefractiveIndex;
    //     const double BetaInverse = c_light/aTrack->GetVelocity();
    //
    //     // below the Cherekov threshold?
    //     if (BetaInverse > maxRefractiveIndex)
    //      return fKill;
    // }

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
            boost::shared_ptr<std::deque<boost::tuple<I3CLSimLightSourceConstPtr, uint32_t, const I3CLSimLightSourceParameterization> > > sendToParameterizationQueue = eventInformation->sendToParameterizationQueue;

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

