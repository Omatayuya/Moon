#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;

class RunAction : public G4UserRunAction
{
public:
	RunAction();
	virtual ~RunAction();

	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);

	void ClearScorerData();

	void AddTrackID(G4int id) { fTrackID.emplace_back(id); };
	void AddNeutronEntranceEnergy(G4double energy) { fNeutronEntranceEnergy.emplace_back(energy); };
	void AddNeutronEntrancePos(G4double x, G4double y, G4double z)
	{
		fNeutronEntrancePosX.emplace_back(x); fNeutronEntrancePosY.emplace_back(y); fNeutronEntrancePosZ.emplace_back(z);
	};

	void AddParentEntranceEnergy(G4double e) { fParentEntranceEnergy.emplace_back(e); };
	void AddParentID(G4int id) { fParentID.emplace_back(id); };
	void AddParentParticleID(G4int id) { fParentParticleID.emplace_back(id); };
	void AddParticleID(G4int id) { fParticleID.emplace_back(id); };
	void AddStepProcessID(G4int id) { fStepProcessID.emplace_back(id); };
	void AddCreatorProcessID(G4int id) { fCreatorProcessID.emplace_back(id); };
	void AddPreEnergy(G4double e) { fPreEnergy.emplace_back(e); };
	void AddPostEnergy(G4double e) { fPostEnergy.emplace_back(e); };
	void AddEdep(G4double e) { fEdep.emplace_back(e); };
	void AddGTime(G4double t) { fGTime.emplace_back(t); };
	void AddPrePos(G4double x, G4double y, G4double z)
	{
		fPrePosX.emplace_back(x); fPrePosY.emplace_back(y); fPrePosZ.emplace_back(z);
	};
	void AddPostPos(G4double x, G4double y, G4double z)
	{
		fPostPosX.emplace_back(x); fPostPosY.emplace_back(y); fPostPosZ.emplace_back(z);
	};
	void AddPreMom(G4double x, G4double y, G4double z)
	{
		fPreMomX.emplace_back(x); fPreMomY.emplace_back(y); fPreMomZ.emplace_back(z);
	};

private:
	std::vector<G4int> fTrackID;
	std::vector<G4double> fNeutronEntranceEnergy;
	std::vector<G4double> fNeutronEntrancePosX, fNeutronEntrancePosY, fNeutronEntrancePosZ;
	
	std::vector<G4double> fParentEntranceEnergy;
	std::vector<G4int> fParentID;
	std::vector<G4int> fParentParticleID;
	std::vector<G4int> fParticleID;
	std::vector<G4int> fStepProcessID;
	std::vector<G4int> fCreatorProcessID;
	std::vector<G4double> fPreEnergy;
	std::vector<G4double> fPostEnergy;
	std::vector<G4double> fEdep;
	std::vector<G4double> fGTime;
	std::vector<G4double> fPrePosX, fPrePosY, fPrePosZ;
	std::vector<G4double> fPostPosX, fPostPosY, fPostPosZ;
	std::vector<G4double> fPreMomX, fPreMomY, fPreMomZ;
};

#endif

