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
#include "G4ThreeVector.hh"

#include "clsim/I3CLSimStepStore.h"
#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimQueue.h"
#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

#include "clsim/I3CLSimLightSource.h"
#include "clsim/I3CLSimLightSourceParameterization.h"

#include <deque>
#include <boost/tuple/tuple.hpp>

#include <boost/thread.hpp>

class G4Event;

class TrkEventAction : public G4UserEventAction
{
public:
    TrkEventAction(uint64_t maxBunchSize,
                   I3CLSimStepStorePtr stepStore,
                   boost::shared_ptr<std::deque<boost::tuple<I3CLSimLightSourceConstPtr, uint32_t, const I3CLSimLightSourceParameterization> > > sendToParameterizationQueue,
                   const I3CLSimLightSourceParameterizationSeries &parameterizationAvailable,
                   boost::shared_ptr<I3CLSimQueue<I3CLSimLightSourceToStepConverterGeant4::FromGeant4Pair_t> > queueFromGeant4,
                   boost::this_thread::disable_interruption &threadDisabledInterruptionState,
                   double maxRefractiveIndex);
    virtual ~TrkEventAction();
    
public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    
    inline bool AbortWasRequested() {return abortRequested_;}
    inline void SetExternalParticleID(uint32_t val) {currentExternalParticleID_=val;}
    
private:
    bool abortRequested_;
    uint64_t maxBunchSize_;
    I3CLSimStepStorePtr stepStore_;
    boost::shared_ptr<std::deque<boost::tuple<I3CLSimLightSourceConstPtr, uint32_t, const I3CLSimLightSourceParameterization> > > sendToParameterizationQueue_;
    uint32_t currentExternalParticleID_;
    
    I3CLSimLightSourceParameterizationSeries parameterizationAvailable_;
    
    boost::shared_ptr<I3CLSimQueue<I3CLSimLightSourceToStepConverterGeant4::FromGeant4Pair_t> > queueFromGeant4_;
    boost::this_thread::disable_interruption &threadDisabledInterruptionState_;
    
    double maxRefractiveIndex_;
};

#endif
