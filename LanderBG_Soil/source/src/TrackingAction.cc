#include "TrackingAction.hh"
#include "EventAction.hh"
#include "G4Track.hh"

TrackingAction::TrackingAction(EventAction *eventAction)
: G4UserTrackingAction(), fEventAction(eventAction)
{}

TrackingAction::~TrackingAction()
{}

void TrackingAction::PreUserTrackingAction(const G4Track *track)
{
	fEventAction->RecordTrackParticle(track->GetTrackID(), track->GetDefinition()->GetParticleName());
}
