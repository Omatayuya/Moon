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
/// \file B2PrimaryGeneratorAction.cc
/// \brief Implementation of the B2PrimaryGeneratorAction class

#include "B2PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
// #include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh" //20211222 HN
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh" // G4RunManager::GetCurrentRun() を使うため
#include "G4Event.hh"      // G4Event::GetEventID() を使うため
#include "G4Run.hh"        // G4Run::GetRunID() を使うため

#include "Randomize.hh"
#include <fstream>
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::B2PrimaryGeneratorAction()
    //: G4VUserPrimaryGeneratorAction()
    : fpParticleGPS(0) // 20211222 HN
{

  // G4int nofParticles = 1;
  // fParticleGun = new G4ParticleGun(nofParticles);
  fpParticleGPS = new G4GeneralParticleSource(); // 20211222 HN

  // default particle kinematic

  // G4ParticleDefinition* particleDefinition
  //   = G4ParticleTable::GetParticleTable()->FindParticle("proton");//yk

  // fParticleGun->SetParticleDefinition(particleDefinition);
  // fGeneralParticleSource->SetParticleDefinition(particleDefinition);      //20211222 HN
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  // fGeneralParticleSource->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.)); //20211222 HN
  // fParticleGun->SetParticleEnergy(10.0*MeV);//yk
  // fGeneralParticleSource->SetParticleEnergy(10.0*MeV);                              //20211222 HN
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::~B2PrimaryGeneratorAction()
{
  // delete fParticleGun;
  delete fpParticleGPS; // 20211222 HN
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore.

  G4double worldZHalfLength = 0;
  G4LogicalVolume *worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box *worldBox = NULL;
  if (worldLV)
    worldBox = dynamic_cast<G4Box *>(worldLV->GetSolid());
  if (worldBox)
    worldZHalfLength = worldBox->GetZHalfLength();
  else
  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  // fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., worldZHalfLength));//yk: z=-worldZHalfLength
  // fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.5 * cm));//HN
  // fGeneralParticleSource->SetParticlePosition(G4ThreeVector(0., 0., 0.5 * cm));//20211222 HN

  // fParticleGun->GeneratePrimaryVertex(anEvent);
  fpParticleGPS->GeneratePrimaryVertex(anEvent); // 20211222 HN

  const G4String outputFileName = "initial_particle.csv";

  static bool isFirstEvent = true;

  std::ofstream outFile(outputFileName, std::ios::app);

  if (!outFile.is_open())
  {
    G4cerr << "Error opening file: " << outputFileName << G4endl;
    return;
  }
  // if (isFirstEvent)
  // {
  //   // outFile.seekp(0, std::ios::end);
  //   outFile << "EventID,ParticleName,Energy(MeV),PositionX(mm),PositionY(mm),PositionZ(mm)\n";
  //   outFile << "Matsunagahayakumathinguapurisiro";
  //   isFirstEvent = false;
  // }

  // イベントIDを取得
  G4int eventID = anEvent->GetEventID();
  G4int currentRunID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();

  // 入射中性子の運動エネルギーを取得
  G4double incidentEnergy = fpParticleGPS->GetCurrentSource()->GetParticleEnergy();

  // 入射位置を取得（GeneratePrimaryVertex 後なので頂点から取得）
  G4PrimaryVertex *primaryVertex = anEvent->GetPrimaryVertex();
  G4double posX = primaryVertex->GetX0() / mm;
  G4double posY = primaryVertex->GetY0() / mm;
  G4double posZ = primaryVertex->GetZ0() / mm;

  // ファイルに書き出し (EventID, Energy[MeV], X[mm], Y[mm], Z[mm])
  outFile
      << eventID
      << "," << incidentEnergy / MeV
      << "," << posX << "," << posY << "," << posZ << std::endl;

  outFile.close(); // ファイルを閉じる
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
