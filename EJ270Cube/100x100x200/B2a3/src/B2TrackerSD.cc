//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2TrackerSD.cc
/// \brief Implementation of the B2TrackerSD class

#include "B2TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"   // yk
#include "G4RunManager.hh" // yk
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::B2TrackerSD(const G4String &name,
                         const G4String &hitsCollectionName)
    : G4VSensitiveDetector(name),
      fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::~B2TrackerSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::Initialize(G4HCofThisEvent *hce)
{
  // Create hits collection

  fHitsCollection = new B2TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B2TrackerSD::ProcessHits(G4Step *aStep,
                                G4TouchableHistory *)
{
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  // if (edep==0.) return false; //yk

  B2TrackerHit *newHit = new B2TrackerHit();

  newHit->SetEventID(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()); // yk
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
  newHit->SetEdep(edep);
  newHit->SetStepLength(aStep->GetStepLength());
  newHit->SetPrePos(aStep->GetPreStepPoint()->GetPosition()); // yk
  newHit->SetPostPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetGlobalTime(aStep->GetTrack()->GetGlobalTime());         // yk
  newHit->SetLocalTime(aStep->GetPostStepPoint()->GetLocalTime());   // yk
  newHit->SetProperTime(aStep->GetPostStepPoint()->GetProperTime()); // yk
  newHit->SetIsEntering(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary);
  newHit->SetDeltaTime(aStep->GetDeltaTime()); // yk

  G4Track *track = aStep->GetTrack();                                 // yk
  const G4DynamicParticle *dynParticle = track->GetDynamicParticle(); // yk
  G4ParticleDefinition *particle = dynParticle->GetDefinition();      // yk

  G4String ParticleName = particle->GetParticleName(); // yk
  G4int ParticleID = track->GetTrackID();              // yk
  G4int ParentParticleID = track->GetParentID();       // yk

  newHit->SetParticleName(ParticleName);         // yk
  newHit->SetParticleID(ParticleID);             // yk
  newHit->SetParentParticleID(ParentParticleID); // yk

  // yk: Process Name of PreStepPoint
  const G4VProcess *aPreProcess = aStep->GetPreStepPoint()->GetProcessDefinedStep();
  if (aPreProcess == 0)
  {
    newHit->SetPreProcessName("1stStep");
  }
  else
  {
    G4String PreProcessName = aPreProcess->GetProcessName();
    newHit->SetPreProcessName(PreProcessName);
  }

  // yk: Process Name of PostStepPoint
  const G4VProcess *aPostProcess = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  G4String PostProcessName = aPostProcess->GetProcessName();
  newHit->SetPostProcessName(PostProcessName);

  // yk: Creator Process Name
  if (ParentParticleID == 0)
  {
    newHit->SetCreatorProcessName("Primary");
  }
  else
  {
    const G4VProcess *anotherProcess = track->GetCreatorProcess();
    G4String CreatorProcessName = anotherProcess->GetProcessName();
    newHit->SetCreatorProcessName(CreatorProcessName);
  }

  // // yk: Creator Model Name
  // if (ParentParticleID == 0)
  // {
  //   newHit->SetCreatorModelName("None");
  // }
  // else
  // {
  //   G4String CreatorModelName = track->GetCreatorModelName();
  //   newHit->SetCreatorModelName(CreatorModelName);
  // }

  G4double PreKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy(); // yk
  newHit->SetPreKineticEnergy(PreKineticEnergy);                            // yk

  G4double PostKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy(); // yk
  newHit->SetPostKineticEnergy(PostKineticEnergy);                            // yk

  G4String VertexVolume = track->GetLogicalVolumeAtVertex()->GetName(); // yk
  newHit->SetLogicalVolumeAtVertexName(VertexVolume);                   // yk

  fHitsCollection->insert(newHit);

  // newHit->Print();

  return true;
}

void B2TrackerSD::FillTree(B2TrackerHit *hit)
{
  auto analysisManager = G4AnalysisManager::Instance();
  const G4int ntupleID = 0; // Hit
  G4int col = 0;

  auto currentEvent = G4RunManager::GetRunManager()->GetCurrentEvent();
  G4double primEnergy = currentEvent->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();

  analysisManager->FillNtupleIColumn(ntupleID, col++, hit->GetEventID());
  analysisManager->FillNtupleIColumn(ntupleID, col++, hit->GetIsEntering() ? 1 : 0);
  analysisManager->FillNtupleDColumn(ntupleID, col++, primEnergy / MeV);
  // analysisManager->FillNtupleSColumn(ntupleID, col++, Form("target%d", hit->GetChamberNb()));
  analysisManager->FillNtupleSColumn(ntupleID, col++, collectionName[0]);
  analysisManager->FillNtupleIColumn(ntupleID, col++, hit->GetParticleID());
  analysisManager->FillNtupleIColumn(ntupleID, col++, hit->GetParentParticleID());
  analysisManager->FillNtupleSColumn(ntupleID, col++, hit->GetParticleName());
  analysisManager->FillNtupleSColumn(ntupleID, col++, hit->GetPreProcessName());
  analysisManager->FillNtupleSColumn(ntupleID, col++, hit->GetPostProcessName());
  analysisManager->FillNtupleSColumn(ntupleID, col++, hit->GetCreatorProcessName());
  // analysisManager->FillNtupleSColumn(ntupleID, col++, hit->GetCreatorModelName());
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetPreKineticEnergy() / MeV);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetPostKineticEnergy() / MeV);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetEdep() / MeV);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetGlobalTime() / ns);
  // analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetLocalTime()/ns);
  // analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetProperTime()/ns);
  // analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetDeltaTime()/ns);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetPrePos().x() / mm);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetPrePos().y() / mm);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetPrePos().z() / mm);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetPostPos().x() / mm);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetPostPos().y() / mm);
  analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetPostPos().z() / mm);
  // analysisManager->FillNtupleDColumn(ntupleID, col++, hit->GetStepLength()/mm);

  analysisManager->AddNtupleRow(ntupleID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::EndOfEvent(G4HCofThisEvent *)
{
  G4int nofHits = fHitsCollection->entries();
  if (nofHits > 0)
  {
    for (G4int i = 0; i < nofHits; i++)
      FillTree((*fHitsCollection)[i]);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
