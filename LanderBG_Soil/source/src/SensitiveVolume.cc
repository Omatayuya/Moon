#include "SensitiveVolume.hh"
#include "TargetHit.hh"
#include "Analysis.hh"

#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"

SensitiveVolume::SensitiveVolume(G4String name, G4String hitsCollectionName)
: G4VSensitiveDetector(name), fHitsCollection(NULL)
{
	collectionName.insert(hitsCollectionName);
}

SensitiveVolume::~SensitiveVolume()
{}

void SensitiveVolume::Initialize(G4HCofThisEvent* hce)
{
	fHitsCollection = new HitsCollection(SensitiveDetectorName, collectionName[0]);

	auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	hce->AddHitsCollection(hcID, fHitsCollection);
}

void SensitiveVolume::EndOfEvent(G4HCofThisEvent*)
{}

G4bool SensitiveVolume::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	//auto analysisManager = G4AnalysisManager::Instance();

	auto aHit = new TargetHit();

	aHit->SetTrackID(aStep->GetTrack()->GetTrackID());
	aHit->SetParentID(aStep->GetTrack()->GetParentID());
	aHit->SetStepNumber(aStep->GetTrack()->GetCurrentStepNumber());
	aHit->SetParticleName(aStep->GetTrack()->GetDefinition()->GetParticleName());
	aHit->SetParticleCharge(aStep->GetTrack()->GetDefinition()->GetPDGCharge());
	//aHit->SetGlobalTime(aStep->GetTrack()->GetGlobalTime());
	aHit->SetEdep(aStep->GetTotalEnergyDeposit());
	aHit->SetStepLength(aStep->GetStepLength());

	if(aHit->GetTrackID()==1){
		aHit->SetCreatorProcessName("Primary");
	}else{
		aHit->SetCreatorProcessName(aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
	}

	// Get pre- and post- step points
	auto preStepPoint = aStep->GetPreStepPoint();
	auto postStepPoint = aStep->GetPostStepPoint();
	auto preTouch = preStepPoint->GetTouchableHandle();
	auto postTouch = postStepPoint->GetTouchableHandle();

	aHit->SetPreKineticEnergy(preStepPoint->GetKineticEnergy());
	aHit->SetPostKineticEnergy(postStepPoint->GetKineticEnergy());
	aHit->SetPreGlobalTime(preStepPoint->GetGlobalTime());
	aHit->SetPostGlobalTime(postStepPoint->GetGlobalTime());
	aHit->SetPreStepPosition(preStepPoint->GetPosition());
	aHit->SetPostStepPosition(postStepPoint->GetPosition());
	aHit->SetPreMomentumDirection(preStepPoint->GetMomentumDirection());
	aHit->SetPostMomentumDirection(postStepPoint->GetMomentumDirection());
	if(preStepPoint->GetStepStatus()==fGeomBoundary) aHit->SetEnterVolume(true);
	if(postStepPoint->GetStepStatus()==fGeomBoundary) aHit->SetLeaveVolume(true);
	aHit->SetStepProcessName(postStepPoint->GetProcessDefinedStep()->GetProcessName());

	aHit->SetPreStepVolumeName(preTouch->GetVolume()->GetName());
	aHit->SetPostStepVolumeName(postTouch->GetVolume()->GetName());
	aHit->SetPreStepVolumeCopyNb(preTouch->GetCopyNumber());
	aHit->SetPostStepVolumeCopyNb(postTouch->GetCopyNumber());

	aHit->SetPreStepLocalPosition(preTouch->GetHistory()->GetTopTransform().TransformPoint(preStepPoint->GetPosition()));
	aHit->SetPostStepLocalPosition(postTouch->GetHistory()->GetTopTransform().TransformPoint(postStepPoint->GetPosition()));

	fHitsCollection->insert(aHit);

	//aHit->Print();

	return true;
}

