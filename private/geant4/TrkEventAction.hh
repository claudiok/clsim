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
 * @file TrkEventAction.hh
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef TrkEventAction_h
#define TrkEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
// #include "G4ThreeVector.hh"

#include "geant4/I3CLSimLightSourcePropagatorGeant4.h"

class G4Event;

class TrkEventAction : public G4UserEventAction
{
public:
    TrkEventAction(double maxRefractiveIndex);
    virtual ~TrkEventAction();
    
public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    
    inline bool AbortWasRequested() {return abortRequested_;}
    inline void SetExternalParticleID(uint32_t val) {currentExternalParticleID_=val;}
    inline void SetSecondaryCallback(const I3CLSimLightSourcePropagator::secondary_callback &callback) {emitSecondary_=callback;}
    inline void SetStepCallback(const I3CLSimLightSourcePropagator::step_callback &callback) {emitStep_=callback;}
    
private:
    bool abortRequested_;
    
    I3CLSimLightSourcePropagator::secondary_callback emitSecondary_;
    I3CLSimLightSourcePropagator::step_callback emitStep_;
    uint32_t currentExternalParticleID_;
    
    double maxRefractiveIndex_;
};

#endif
