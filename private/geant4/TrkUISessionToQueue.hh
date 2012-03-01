#ifndef TrkUISessionToQueue_H
#define TrkUISessionToQueue_H 1

#include "globals.hh"
#include "G4UIsession.hh"

#include <boost/shared_ptr.hpp>
#include "clsim/I3CLSimQueue.h"
#include <string>

class TrkUISessionToQueue : public G4UIsession
{
public:
    TrkUISessionToQueue(boost::shared_ptr<I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > > > queueFromGeant4Messages);
    ~TrkUISessionToQueue();

    // These two methods will be invoked by G4strstreambuf.
    virtual G4int ReceiveG4cout(G4String coutString);
    virtual G4int ReceiveG4cerr(G4String cerrString);

private:
    boost::shared_ptr<I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > > > queueFromGeant4Messages_;
};

#endif
