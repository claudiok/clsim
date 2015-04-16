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
 * @file TrkUISessionToQueue.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "TrkUISessionToQueue.hh"


TrkUISessionToQueue::TrkUISessionToQueue(boost::shared_ptr<I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > > > queueFromGeant4Messages)
:
queueFromGeant4Messages_(queueFromGeant4Messages)
{
    
}

TrkUISessionToQueue::~TrkUISessionToQueue()
{
    
}

// These two methods will be invoked by G4strstreambuf.
#if G4VERSION_NUMBER >= 960
G4int TrkUISessionToQueue::ReceiveG4cout(const G4String& coutString)
#else
G4int TrkUISessionToQueue::ReceiveG4cout(G4String coutString)
#endif
{
    if (queueFromGeant4Messages_)
    {
        boost::shared_ptr<std::pair<const std::string, bool> > strOut(new std::pair<const std::string, bool>(std::string(coutString), false));
        queueFromGeant4Messages_->Put(strOut);
    }
    else
    {
        std::cout << coutString << std::flush;
    }
    
    return 0;   
}

#if G4VERSION_NUMBER >= 960
G4int TrkUISessionToQueue::ReceiveG4cerr(const G4String& cerrString)
#else
G4int TrkUISessionToQueue::ReceiveG4cerr(G4String cerrString)
#endif
{
    if (queueFromGeant4Messages_)
    {
        boost::shared_ptr<std::pair<const std::string, bool> > strOut(new std::pair<const std::string, bool>(std::string(cerrString), true));
        queueFromGeant4Messages_->Put(strOut);
    }
    else
    {
        std::cout << cerrString << std::flush;
    }
    
    return 0;
}

