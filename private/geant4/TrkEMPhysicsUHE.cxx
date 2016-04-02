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
 * @file TrkEMPhysicsUHE.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "TrkEMPhysicsUHE.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>  

#include "G4ProcessManager.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4OpticalPhoton.hh"

#include "G4LossTableManager.hh"

TrkEMPhysicsUHE::TrkEMPhysicsUHE(const G4String& name)
  :  G4VPhysicsConstructor(name)
{
    G4cout << "<<<< UHE EM Physics and Loss Table Extension" << G4endl;
}

TrkEMPhysicsUHE::~TrkEMPhysicsUHE()
{
}


void TrkEMPhysicsUHE::ConstructParticle()
{
    // no additional particles
}

void TrkEMPhysicsUHE::ConstructProcess()
{
    // Construct Processes
    
    //extend binning of PhysicsTables
    G4LossTableManager::Instance()->SetMaxEnergy(10.*PeV);
    G4LossTableManager::Instance()->SetDEDXBinning(440);
    G4LossTableManager::Instance()->SetLambdaBinning(440);
    G4LossTableManager::Instance()->SetVerbose(0);
}
