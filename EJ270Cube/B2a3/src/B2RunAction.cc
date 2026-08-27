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
/// \file B2RunAction.cc
/// \brief Implementation of the B2RunAction class

#include "B2RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "B2PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2RunAction::B2RunAction()
    : G4UserRunAction()
{
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetNtupleMerging(true);

  // Ntuple 0: 初期粒子情報
  const G4int nThetaBins = 18;
  const G4double thetaBinWidth = 90.0 / nThetaBins;
  const G4int nEnergyBins = 500;
  const G4double energyMin = 1e-10; // MeV
  const G4double energyMax = 1e4;   // MeV
  for (G4int i = 0; i < nThetaBins; ++i)
  {
    G4String name = "hPrimEnergyByTheta_" + std::to_string(i);
    G4String title = std::to_string(i * thetaBinWidth) + " - " + std::to_string((i + 1) * thetaBinWidth) + " deg";
    analysisManager->CreateH1(name, title, nEnergyBins, energyMin, energyMax, "none", "none", "log");
  }
  analysisManager->CreateH1("hPrimEnergy", "Primary Energy ; energy [MeV]; Entries", nEnergyBins, energyMin, energyMax, "none", "none", "log");
  analysisManager->CreateH1("hPrimTheta", "Primary Theta ; theta [deg]; Entries", 90, 0., 90.);
  analysisManager->CreateH1("hPrimPhi", "Primary Phi ; phi [deg]; Entries", 360, 0., 360.);
  analysisManager->CreateH1("hPrimPosZ", "Primary Position Z ; position [mm]; Entries", 100, -50., 50.);
  analysisManager->CreateH2("hPrimPosXY", "Primary Position XY ; position X [mm]; position Y [mm]; Entries", 100, -50., 50., 100, -50., 50.);

  // Ntuple 1: ヒット情報（本体+並行世界を共有）
  analysisManager->CreateNtuple("Hit", "tracker/scorer hit info");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleIColumn("IsEntering");
  analysisManager->CreateNtupleDColumn("primEnergy");
  // analysisManager->CreateNtupleSColumn("chamberNb");
  analysisManager->CreateNtupleSColumn("collection");
  analysisManager->CreateNtupleIColumn("PID");
  analysisManager->CreateNtupleIColumn("PPID");
  analysisManager->CreateNtupleSColumn("Pname");
  analysisManager->CreateNtupleSColumn("PreProc");
  analysisManager->CreateNtupleSColumn("PostProc");
  analysisManager->CreateNtupleSColumn("CProc");
  // analysisManager->CreateNtupleSColumn("CModel");
  analysisManager->CreateNtupleDColumn("PreKinE");
  analysisManager->CreateNtupleDColumn("PostKinE");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("GTime");
  // analysisManager->CreateNtupleDColumn("LTime");
  // analysisManager->CreateNtupleDColumn("PTime");
  // analysisManager->CreateNtupleDColumn("DTime");
  analysisManager->CreateNtupleDColumn("PrePosX");
  analysisManager->CreateNtupleDColumn("PrePosY");
  analysisManager->CreateNtupleDColumn("PrePosZ");
  analysisManager->CreateNtupleDColumn("PostPosX");
  analysisManager->CreateNtupleDColumn("PostPosY");
  analysisManager->CreateNtupleDColumn("PostPosZ");
  // analysisManager->CreateNtupleDColumn("StepLength");
  analysisManager->FinishNtuple();

  // Ntuple 2: Run情報
  fRunInfoNtupleId = analysisManager->CreateNtuple("RunInfo", "Run information");
  analysisManager->CreateNtupleDColumn(fRunInfoNtupleId, "TotalFlux");
  analysisManager->FinishNtuple(fRunInfoNtupleId);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2RunAction::~B2RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2RunAction::BeginOfRunAction(const G4Run *run)
{
  G4int currentRunID = run->GetRunID();
  G4cout << "### Starting Run " << currentRunID << " ###" << G4endl;
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile("results");

  auto *gen = dynamic_cast<const B2PrimaryGeneratorAction *>(
      G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if (gen)
  {
    analysisManager->FillNtupleDColumn(fRunInfoNtupleId, 0, gen->GetTotalFlux());
    analysisManager->AddNtupleRow(fRunInfoNtupleId);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2RunAction::EndOfRunAction(const G4Run *)
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
