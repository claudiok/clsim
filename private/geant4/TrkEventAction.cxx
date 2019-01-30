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
 * @file TrkEventAction.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

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

TrkEventAction::TrkEventAction(double maxRefractiveIndex)
:
abortRequested_(false),
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
    new TrkUserEventInformation(emitSecondary_,
                                emitStep_,
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
