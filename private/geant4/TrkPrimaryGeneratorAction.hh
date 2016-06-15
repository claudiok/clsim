#ifndef TRKPRIMARYGENERATORACTION_H_INCLUDED
#define TRKPRIMARYGENERATORACTION_H_INCLUDED
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
 * @file TrkPrimaryGeneratorAction.hh
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef TrkPrimaryGeneratorAction_h
#define TrkPrimaryGeneratorAction_h 1

#include "G4GeneralParticleSource.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class TrkPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    TrkPrimaryGeneratorAction();
    virtual ~TrkPrimaryGeneratorAction();

public:
    void GeneratePrimaries(G4Event* anEvent);

    inline G4ParticleGun *GetParticleGun() {return particleGun;}

private:
    G4ParticleGun* particleGun;
};

#endif



#endif  // TRKPRIMARYGENERATORACTION_H_INCLUDED
