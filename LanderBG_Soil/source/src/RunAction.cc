#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

RunAction::RunAction()
	: G4UserRunAction()
{
	// Create analysis manager
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetDefaultFileType("root");
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetNtupleMerging(true);
	analysisManager->SetNtupleRowWise(false, true);

	// Create histograms
	analysisManager->CreateH1("PrimaryEnergy", "Primary energy; Energy (MeV); Counts", 800, 1, 1e+8, "none", "none", "log");
	analysisManager->CreateH1("PrimaryX", "Primary position X; X (m); Counts", 320, -16, 16);
	analysisManager->CreateH1("PrimaryY", "Primary position Y; Y (m); Counts", 320, -16, 16);
	analysisManager->CreateH1("PrimaryZ", "Primary position Z; Z (m); Counts", 120, -6, 6);
	analysisManager->CreateH1("PrimaryTheta", "Primary polar angle; Theta (deg); Counts", 180, 0, 180);
	analysisManager->CreateH1("PrimaryPhi", "Primary azimuth angle; Phi (deg); Counts", 360, -180, 180);
	analysisManager->CreateH2("PrimaryXY", "Primary XY; X (m); Y (m); Counts", 320, -16, 16, 320, -16, 16);
	analysisManager->CreateH2("PrimaryXZ", "Primary XZ; X (m); Z (m); Counts", 320, -16, 16, 120, -6, 6);

	// neutronEntrance: EJ-270に入射した中性子(反応の有無によらず全件)
	G4int ntupleIdx = 0;
	analysisManager->CreateNtuple("neutronEntrance", "neutronEntrance");
	analysisManager->CreateNtupleIColumn(ntupleIdx, "eventID");
	analysisManager->CreateNtupleIColumn(ntupleIdx, "trackID", fTrackID);
	analysisManager->CreateNtupleDColumn(ntupleIdx, "energy", fNeutronEntranceEnergy);
	analysisManager->CreateNtupleDColumn(ntupleIdx, "posX", fNeutronEntrancePosX);
	analysisManager->CreateNtupleDColumn(ntupleIdx, "posY", fNeutronEntrancePosY);
	analysisManager->CreateNtupleDColumn(ntupleIdx, "posZ", fNeutronEntrancePosZ);
	analysisManager->FinishNtuple(ntupleIdx++);

	// Create ntuple (1 row = 1 event, per-hit values stored as vectors)
	analysisManager->CreateNtuple("scorer", "scorer");
	analysisManager->CreateNtupleIColumn(ntupleIdx, "eventID");
	analysisManager->CreateNtupleIColumn(ntupleIdx, "parentID", fParentID);                   // 3
	analysisManager->CreateNtupleDColumn(ntupleIdx, "parentEntranceEnergy", fParentEntranceEnergy); // 4
	analysisManager->CreateNtupleIColumn(ntupleIdx, "parentParticleID", fParentParticleID);   // 5
	analysisManager->CreateNtupleIColumn(ntupleIdx, "particleID", fParticleID);               // 6
	analysisManager->CreateNtupleIColumn(ntupleIdx, "stepProcessID", fStepProcessID);         // 7
	analysisManager->CreateNtupleIColumn(ntupleIdx, "creatorProcessID", fCreatorProcessID);   // 8
	analysisManager->CreateNtupleDColumn(ntupleIdx, "preEnergy", fPreEnergy);                 // 9
	analysisManager->CreateNtupleDColumn(ntupleIdx, "postEnergy", fPostEnergy);               // 10
	analysisManager->CreateNtupleDColumn(ntupleIdx, "edep", fEdep);                           // 11
	analysisManager->CreateNtupleDColumn(ntupleIdx, "gTime", fGTime);                         // 12
	analysisManager->CreateNtupleDColumn(ntupleIdx, "prePosX", fPrePosX);                     // 13
	analysisManager->CreateNtupleDColumn(ntupleIdx, "prePosY", fPrePosY);                     // 14
	analysisManager->CreateNtupleDColumn(ntupleIdx, "prePosZ", fPrePosZ);                     // 15
	analysisManager->CreateNtupleDColumn(ntupleIdx, "postPosX", fPostPosX);                   // 16
	analysisManager->CreateNtupleDColumn(ntupleIdx, "postPosY", fPostPosY);                   // 17
	analysisManager->CreateNtupleDColumn(ntupleIdx, "postPosZ", fPostPosZ);                   // 18
	analysisManager->CreateNtupleDColumn(ntupleIdx, "preMomX", fPreMomX);                     // 19
	analysisManager->CreateNtupleDColumn(ntupleIdx, "preMomY", fPreMomY);                     // 20
	analysisManager->CreateNtupleDColumn(ntupleIdx, "preMomZ", fPreMomZ);                     // 21
	analysisManager->FinishNtuple(ntupleIdx++);
}

RunAction::~RunAction()
{
}

void RunAction::ClearScorerData()
{
	fTrackID.clear();
	fNeutronEntranceEnergy.clear();
	fNeutronEntrancePosX.clear(); fNeutronEntrancePosY.clear(); fNeutronEntrancePosZ.clear();
	
	fParentEntranceEnergy.clear();
	fParentID.clear();
	fParentParticleID.clear();
	fParticleID.clear();
	fStepProcessID.clear();
	fCreatorProcessID.clear();
	fPreEnergy.clear();
	fPostEnergy.clear();
	fEdep.clear();
	fGTime.clear();
	fPrePosX.clear(); fPrePosY.clear(); fPrePosZ.clear();
	fPostPosX.clear(); fPostPosY.clear(); fPostPosZ.clear();
	fPreMomX.clear(); fPreMomY.clear(); fPreMomZ.clear();
}

void RunAction::BeginOfRunAction(const G4Run *)
{
	// Get analysis manager and open an output file
	auto analysisManager = G4AnalysisManager::Instance();

	G4String fileName = "results";
	analysisManager->OpenFile(fileName);
}

void RunAction::EndOfRunAction(const G4Run *)
{
	// Save histograms & ntuple
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();

	G4cout << "End of Runaction" << G4endl;
}
