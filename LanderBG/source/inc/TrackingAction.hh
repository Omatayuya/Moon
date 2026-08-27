#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class EventAction;

class TrackingAction : public G4UserTrackingAction
{
public:
	TrackingAction(EventAction *eventAction);
	virtual ~TrackingAction();

	virtual void PreUserTrackingAction(const G4Track *track);

private:
	EventAction *fEventAction;
};

#endif
