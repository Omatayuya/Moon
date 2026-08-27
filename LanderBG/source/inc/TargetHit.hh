#ifndef TargetHit_h
#define TargetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class TargetHit : public G4VHit
{
public:
	TargetHit();
	TargetHit(const TargetHit&);
	virtual ~TargetHit();

	// operators
	const TargetHit& operator=(const TargetHit&);
	G4int operator==(const TargetHit&) const;

	inline void* operator new(size_t);
	inline void operator delete(void*);

	// methods from base class
	virtual void Draw();
	virtual void Print();

	// Set methods
	void SetTrackID(G4int track){ fTrackID = track; };
	void SetParentID(G4int parent){ fParentID = parent; };
	void SetStepNumber(G4int step){ fStepNb = step; };
	void SetParticleName(G4String particle){ fParticleName = particle; };
	void SetParticleCharge(G4double charge){ fParticleCharge = charge; };
	void SetStepProcessName(G4String stepProcess){ fStepProcessName = stepProcess; };
	//void SetGlobalTime(G4double globalTime){ fGlobalTime = globalTime; };
	void SetCreatorProcessName(G4String creatorProcess){ fCreatorProcessName = creatorProcess; };
	void SetEdep(G4double de){ fEdep = de; };
	void SetStepLength(G4double length){ fStepLength = length; };
	void SetPreKineticEnergy(G4double ke){ fPreKineticEnergy = ke; };
	void SetPostKineticEnergy(G4double ke){ fPostKineticEnergy = ke; };
	void SetPreGlobalTime(G4double globalTime){ fPreGlobalTime = globalTime; };
	void SetPostGlobalTime(G4double globalTime){ fPostGlobalTime = globalTime; };
	void SetPreStepPosition(G4ThreeVector xyz){ fPreStepPosition = xyz; };
	void SetPostStepPosition(G4ThreeVector xyz){ fPostStepPosition = xyz; };
	void SetPreMomentumDirection(G4ThreeVector xyz){ fPreMomentumDirection = xyz; };
	void SetPostMomentumDirection(G4ThreeVector xyz){ fPostMomentumDirection = xyz; };
	void SetEnterVolume(G4bool bound){ fEnterVolume = bound; };
	void SetLeaveVolume(G4bool bound){ fLeaveVolume = bound; };
	void SetPreStepVolumeName(G4String volName){ fPreStepVolumeName = volName; };
	void SetPostStepVolumeName(G4String volName){ fPostStepVolumeName = volName; };
	void SetPreStepVolumeCopyNb(G4int copyNb){ fPreStepVolumeCopyNb = copyNb; };
	void SetPostStepVolumeCopyNb(G4int copyNb){ fPostStepVolumeCopyNb = copyNb; };
	void SetPreStepLocalPosition(G4ThreeVector xyz){ fPreStepLocalPosition = xyz; };
	void SetPostStepLocalPosition(G4ThreeVector xyz){ fPostStepLocalPosition = xyz; };

	// Get methods
	G4int GetTrackID() const { return fTrackID; };
	G4int GetParentID() const { return fParentID; };
	G4int GetStepNumber() const { return fStepNb; };
	G4String GetParticleName() const { return fParticleName; };
	G4double GetParticleCharge() const { return fParticleCharge; };
	G4String GetStepProcessName() const { return fStepProcessName; };
	//G4double GetGlobalTime() const { return fGlobalTime; };
	G4String GetCreatorProcessName() const { return fCreatorProcessName; };
	G4double GetEdep() const { return fEdep; };
	G4double GetStepLength() const { return fStepLength; };
	G4double GetPreKineticEnergy() const { return fPreKineticEnergy; };
	G4double GetPostKineticEnergy() const { return fPostKineticEnergy; };
	G4double GetPreGlobalTime() const { return fPreGlobalTime; };
	G4double GetPostGlobalTime() const { return fPostGlobalTime; };
	G4ThreeVector GetPreStepPosition() const { return fPreStepPosition; };
	G4ThreeVector GetPostStepPosition() const { return fPostStepPosition; };
	G4ThreeVector GetPreMomentumDirection() const { return fPreMomentumDirection; };
	G4ThreeVector GetPostMomentumDirection() const { return fPostMomentumDirection; };
	G4bool isEnteringVolume() const { return fEnterVolume; };
	G4bool isLeavingVolume() const { return fLeaveVolume; };
	G4String GetPreStepVolumeName() const { return fPreStepVolumeName; };
	G4String GetPostStepVolumeName() const { return fPostStepVolumeName; };
	G4int GetPreStepVolumeCopyNb() const { return fPreStepVolumeCopyNb; };
	G4int GetPostStepVolumeCopyNb() const { return fPostStepVolumeCopyNb; };
	G4ThreeVector GetPreStepLocalPosition() const { return fPreStepLocalPosition; };
	G4ThreeVector GetPostStepLocalPosition() const { return fPostStepLocalPosition; };

private:
	G4int fTrackID;
	G4int fParentID;
	G4int fStepNb;
	G4String fParticleName;
	G4double fParticleCharge;
	G4String fStepProcessName;
	//G4double fGlobalTime;
	G4double fEdep;
	G4double fStepLength;
	G4String fCreatorProcessName;
	G4double fPreKineticEnergy;
	G4double fPostKineticEnergy;
	G4double fPreGlobalTime;
	G4double fPostGlobalTime;
	G4ThreeVector fPreStepPosition;
	G4ThreeVector fPostStepPosition;
	G4ThreeVector fPreMomentumDirection;
	G4ThreeVector fPostMomentumDirection;
	G4bool fEnterVolume;
	G4bool fLeaveVolume;
	G4String fPreStepVolumeName;
	G4String fPostStepVolumeName;
	G4int fPreStepVolumeCopyNb;
	G4int fPostStepVolumeCopyNb;
	G4ThreeVector fPreStepLocalPosition;
	G4ThreeVector fPostStepLocalPosition;
};

typedef G4THitsCollection<TargetHit> HitsCollection;

extern G4ThreadLocal G4Allocator<TargetHit>* HitAllocator;

inline void* TargetHit::operator new(size_t){
	if(!HitAllocator) HitAllocator = new G4Allocator<TargetHit>;

	return (void*)HitAllocator->MallocSingle();
}

inline void TargetHit::operator delete(void *hit)
{
	HitAllocator->FreeSingle((TargetHit*) hit);
}

#endif

