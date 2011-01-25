#ifndef TrkStackingAction_H
#define TrkStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class TrkStackingAction : public G4UserStackingAction
{
    public:
	TrkStackingAction();
	~TrkStackingAction();
	
	virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
	virtual void NewStage();
	virtual void PrepareNewEvent();
	
    private:
};

#endif
