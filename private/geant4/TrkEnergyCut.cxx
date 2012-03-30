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
 * @file TrkEnergyCut.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "TrkEnergyCut.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>  

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4UserSpecialCuts.hh"

TrkEnergyCut::TrkEnergyCut(const G4String& name)
  :  G4VPhysicsConstructor(name)
{
    G4cout << "<<<< Energy cut Process (kills particles below an energy threshold)" << G4endl;
}

TrkEnergyCut::~TrkEnergyCut()
{
}


void TrkEnergyCut::ConstructParticle()
{
    // no particles to construct
}

void TrkEnergyCut::ConstructProcess()
{
    // Construct Processes
    theSpecialCutProcess = new G4UserSpecialCuts();
    
    // Add the processes to their respective particles
    G4ProcessManager * pManager = 0;
    
    theParticleIterator->reset();
    while ((*theParticleIterator)())
    {
        G4ParticleDefinition* particle = theParticleIterator->value();
        pManager = particle->GetProcessManager();
        if (theSpecialCutProcess->IsApplicable(*particle))
        {
            pManager->AddProcess(theSpecialCutProcess,-1,-1,9000);
        }
    }

}
