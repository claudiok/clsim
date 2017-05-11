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
 * @file TrkOpticalPhysics.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "TrkOpticalPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>  

#include "G4ProcessManager.hh"

#include "TrkCerenkov.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4OpticalPhoton.hh"
#include "G4Version.hh"

TrkOpticalPhysics::TrkOpticalPhysics(const G4String& name, 
                                     double maxBetaChangePerStep,
                                     uint32_t maxNumPhotonsPerStep,
                                     I3CLSimFunctionConstPtr wlenBias
                                     )
:  G4VPhysicsConstructor(name),
maxBetaChangePerStep_(maxBetaChangePerStep),
maxNumPhotonsPerStep_(maxNumPhotonsPerStep),
wlenBias_(wlenBias)
{
    G4cout << "<<<< Optical Processes (TrkCerenkov)" << G4endl;
}

TrkOpticalPhysics::~TrkOpticalPhysics()
{
}


void TrkOpticalPhysics::ConstructParticle()
{
}

void TrkOpticalPhysics::ConstructProcess()
{
    // Construct Processes
    
    theCerenkovProcess=new TrkCerenkov();

    theCerenkovProcess->SetMaxBetaChangePerStep(maxBetaChangePerStep_);
    theCerenkovProcess->SetMaxNumPhotonsPerStep(maxNumPhotonsPerStep_);
    theCerenkovProcess->SetWlenBiasFunction(wlenBias_);
    
    // Add the processes to their respective particles
    G4ProcessManager * pManager = NULL;
    
    
#if G4VERSION_NUMBER >= 1000
    G4ParticleTable::G4PTblDicIterator* theParticleIterator;
    theParticleIterator = theParticleTable->GetIterator();
#endif
    theParticleIterator->reset();
    while ((*theParticleIterator)())
    {
        G4ParticleDefinition* particle = theParticleIterator->value();
        pManager = particle->GetProcessManager();
        if (theCerenkovProcess->IsApplicable(*particle))
        {
            pManager->AddProcess(theCerenkovProcess);
            pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
        }
    }
    
}
