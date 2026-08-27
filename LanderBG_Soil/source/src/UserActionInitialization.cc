#include "UserActionInitialization.hh"
#include "PrimaryGenerator.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"

UserActionInitialization::UserActionInitialization()
: G4VUserActionInitialization()
{}

UserActionInitialization::~UserActionInitialization()
{}

void UserActionInitialization::Build() const
{
	SetUserAction(new PrimaryGenerator());

	auto runAction = new RunAction();
	SetUserAction(runAction);

	auto eventAction = new EventAction(runAction);
	SetUserAction(eventAction);

	SetUserAction(new TrackingAction(eventAction));
}

void UserActionInitialization::BuildForMaster() const
{
	SetUserAction(new RunAction());
}

