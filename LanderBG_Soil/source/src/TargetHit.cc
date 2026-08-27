#include "TargetHit.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<TargetHit>* HitAllocator = 0;

TargetHit::TargetHit()
: G4VHit(),
	fTrackID(0),
	fParentID(0),
	fStepNb(0),
	fParticleName(G4String()),
	fParticleCharge(0.),
	fStepProcessName(G4String()),
	//fGlobalTime(0.),
	fEdep(0.),
	fStepLength(0.),
	fCreatorProcessName(G4String()),
	fPreKineticEnergy(0.),
	fPostKineticEnergy(0.),
	fPreGlobalTime(0.),
	fPostGlobalTime(0.),
	fPreStepPosition(G4ThreeVector()),
	fPostStepPosition(G4ThreeVector()),
	fPreMomentumDirection(G4ThreeVector()),
	fPostMomentumDirection(G4ThreeVector()),
	fEnterVolume(false),
	fLeaveVolume(false),
	fPreStepVolumeName(G4String()),
	fPostStepVolumeName(G4String()),
	fPreStepVolumeCopyNb(-1),
	fPostStepVolumeCopyNb(-1),
	fPreStepLocalPosition(G4ThreeVector()),
	fPostStepLocalPosition(G4ThreeVector())
{}

TargetHit::~TargetHit()
{}

TargetHit::TargetHit(const TargetHit& /*right*/)
: G4VHit()
{}

const TargetHit& TargetHit::operator=(const TargetHit& /*right*/)
{
	return *this;
}

G4int TargetHit::operator==(const TargetHit& right) const
{
	return (this==&right) ? 1:0;
}

void TargetHit::Draw()
{}

void TargetHit::Print()
{
	G4cout << "-------------------------------" << G4endl;
	G4cout << "trackID: " << fTrackID
			<< ", particleName: " << fParticleName
			<< ", particleCharge: " << fParticleCharge
			<< ", parentID: " << fParentID
			<< ", creatorProcess: " << fCreatorProcessName
			<< G4endl;
	G4cout << "stepNumber: " << fStepNb
			<< ", stepProcess: " << fStepProcessName
			<< ", stepLength: " << std::setw(6) << G4BestUnit(fStepLength, "Length")
			<< ", Edep: " << std::setw(6) << G4BestUnit(fEdep, "Energy")
			<< G4endl;

	G4cout << "preKineticEnergy: " << std::setw(6) << G4BestUnit(fPreKineticEnergy, "Energy")
			<< ", preGlobalTime: " << std::setw(6) << G4BestUnit(fPreGlobalTime, "Time")
			<< G4endl;
	G4cout << "preStepPoint: " << std::setw(6) << G4BestUnit(fPreStepPosition, "Length")
			<< ", preStepLocalPoint: " << std::setw(6) << G4BestUnit(fPreStepLocalPosition, "Length")
			<< G4endl; 
	G4cout << "preMomentumDirection: " << fPreMomentumDirection
			<< G4endl;
	G4cout << "isEnteringVolume: " << fEnterVolume
			<< ", preStepVolumeName: " << fPreStepVolumeName
			<< ", preStepVolumeCopyNumber: " << fPreStepVolumeCopyNb
			<< G4endl;

	G4cout << "postKineticEnergy: " << std::setw(6) << G4BestUnit(fPostKineticEnergy, "Energy")
			<< ", postGlobalTime: " << std::setw(6) << G4BestUnit(fPostGlobalTime, "Time")
			<< G4endl;
	G4cout << "postStepPoint: " << std::setw(6) << G4BestUnit(fPostStepPosition, "Length")
			<< ", postStepLocalPoint: " << std::setw(6) << G4BestUnit(fPostStepLocalPosition, "Length")
			<< G4endl;
	G4cout << "postMomentumDirection: " << fPostMomentumDirection
			<< G4endl;
	G4cout << "isLeavingVolume: " << fLeaveVolume
			<< ", postStepVolumeName: " << fPostStepVolumeName
			<< ", postStepVolumeCopyNumber: " << fPostStepVolumeCopyNb
			<< G4endl;
	G4cout << "-------------------------------" << G4endl;
}

