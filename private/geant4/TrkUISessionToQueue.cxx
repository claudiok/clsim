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
G4int TrkUISessionToQueue::ReceiveG4cout(G4String coutString)
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

G4int TrkUISessionToQueue::ReceiveG4cerr(G4String cerrString)
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

