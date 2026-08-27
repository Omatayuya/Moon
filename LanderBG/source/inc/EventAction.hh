#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <map>

#include "TargetHit.hh"

class RunAction;

class EventAction : public G4UserEventAction
{
public:
	EventAction(RunAction *runAction);
	virtual ~EventAction();

	virtual void BeginOfEventAction(const G4Event *);
	virtual void EndOfEventAction(const G4Event *);

	void RecordTrackParticle(G4int trackID, const G4String &name) { fTrackIDToParticleNameMap[trackID] = name; };
	G4String GetTrackParticleName(G4int trackID) const
	{
		auto it = fTrackIDToParticleNameMap.find(trackID);
		return (it != fTrackIDToParticleNameMap.end()) ? it->second : "Primary";
	};

	G4int ParticleID(G4String particle);
	G4int ProcessID(G4String process);

private:
	HitsCollection *GetHitsCollection(G4int hcID, const G4Event *event) const;
	void PrintEventStatistics() const;

	RunAction *fRunAction;
	G4int fTargetHCID;

	std::map<G4int, G4String> fTrackIDToParticleNameMap;
};

#endif
