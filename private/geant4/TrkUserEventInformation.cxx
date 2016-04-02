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
 * @file TrkUserEventInformation.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "TrkUserEventInformation.hh"
#include "TrkDetectorConstruction.hh"

TrkUserEventInformation::TrkUserEventInformation(uint64_t maxBunchSize_,
                                                 I3CLSimStepStorePtr stepStore_,
                                                 boost::shared_ptr<std::deque<boost::tuple<I3CLSimLightSourceConstPtr, uint32_t, const I3CLSimLightSourceParameterization> > > sendToParameterizationQueue_,
                                                 const I3CLSimLightSourceParameterizationSeries &parameterizationAvailable_,
                                                 boost::shared_ptr<I3CLSimQueue<I3CLSimLightSourceToStepConverterGeant4::FromGeant4Pair_t> > queueFromGeant4_,
                                                 boost::this_thread::disable_interruption &threadDisabledInterruptionState_,
                                                 uint32_t currentExternalParticleID_,
                                                 double maxRefractiveIndex_)
:
abortRequested(false),
maxBunchSize(maxBunchSize_),
stepStore(stepStore_),
sendToParameterizationQueue(sendToParameterizationQueue_),
parameterizationAvailable(parameterizationAvailable_),
queueFromGeant4(queueFromGeant4_),
threadDisabledInterruptionState(threadDisabledInterruptionState_),
currentExternalParticleID(currentExternalParticleID_),
maxRefractiveIndex(maxRefractiveIndex_)
{
}

TrkUserEventInformation::~TrkUserEventInformation()
{
}

