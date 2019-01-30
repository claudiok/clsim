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

    const G4double trackEnergy = aTrack->GetKineticEnergy();
    
    // FIXME: find a way to query for parameterizations / other propagators
    // without copying G4 data into an I3Particle
    
    I3Particle particle;
    
    const G4ThreeVector &trackPos = aTrack->GetPosition();
    const G4double trackTime = aTrack->GetGlobalTime();
    const G4ThreeVector &trackDir = aTrack->GetMomentumDirection();

    particle.SetPdgEncoding(aTrack->GetDefinition()->GetPDGEncoding());
    particle.SetPos(trackPos.x()*I3Units::m/m,trackPos.y()*I3Units::m/m,trackPos.z()*I3Units::m/m);
    particle.SetDir(trackDir.x(),trackDir.y(),trackDir.z());
    particle.SetTime(trackTime*I3Units::ns/ns);
    particle.SetEnergy(trackEnergy*I3Units::GeV/GeV);

    I3CLSimLightSourceConstPtr lightSource = boost::make_shared<I3CLSimLightSource>(particle);
    
    // Forget this particle if it can be handled elsewhere
    if (eventInformation->emitSecondary(lightSource, eventInformation->currentExternalParticleID))
        return fKill;
    else
        return fUrgent;
}

void TrkStackingAction::NewStage()
{
}

void TrkStackingAction::PrepareNewEvent()
{
}

