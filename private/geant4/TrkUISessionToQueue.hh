#ifndef TRKUISESSIONTOQUEUE_H_INCLUDED
#define TRKUISESSIONTOQUEUE_H_INCLUDED
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
 * @file TrkUISessionToQueue.hh
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef TrkUISessionToQueue_H
#define TrkUISessionToQueue_H 1

#include "globals.hh"
#include "G4UIsession.hh"
#include "G4Version.hh"

#include <boost/shared_ptr.hpp>
#include "clsim/I3CLSimQueue.h"
#include <string>

class TrkUISessionToQueue : public G4UIsession
{
public:
    TrkUISessionToQueue(boost::shared_ptr<I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > > > queueFromGeant4Messages);
    virtual ~TrkUISessionToQueue();

    // These two methods will be invoked by G4strstreambuf.
#if G4VERSION_NUMBER >= 960
    // yay! Geant4 changed the interface in version 4.9.6
    virtual G4int ReceiveG4cout(const G4String& coutString);
    virtual G4int ReceiveG4cerr(const G4String& cerrString);
#else
    virtual G4int ReceiveG4cout(G4String coutString);
    virtual G4int ReceiveG4cerr(G4String cerrString);
#endif

private:
    boost::shared_ptr<I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > > > queueFromGeant4Messages_;
};

#endif

#endif  // TRKUISESSIONTOQUEUE_H_INCLUDED
