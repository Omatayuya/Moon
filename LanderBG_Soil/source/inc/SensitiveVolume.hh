#ifndef SensitiveVolume_h
#define SensitiveVolume_h 1

#include "TargetHit.hh"
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

class SensitiveVolume : public G4VSensitiveDetector
{
public:
	SensitiveVolume(G4String name, G4String hitsCollectionName);
	~SensitiveVolume();

	void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
	void EndOfEvent(G4HCofThisEvent*);

private:
	HitsCollection* fHitsCollection;
};

#endif

