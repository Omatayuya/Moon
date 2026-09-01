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
/// \file B2aDetectorConstruction.cc
/// \brief Implementation of the B2aDetectorConstruction class

#include "B2aDetectorConstruction.hh"
#include "B2aDetectorMessenger.hh"
#include "B2TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger *B2aDetectorConstruction::fMagFieldMessenger = 0;

B2aDetectorConstruction::B2aDetectorConstruction()
    : G4VUserDetectorConstruction(),
      fNbOfChambers(0),
      fLogicTarget(NULL), fLogicTarget_0(NULL), fLogicTarget_1(NULL), fLogicTarget_2(NULL),
      fLogicTarget_3(NULL), fLogicTarget_4(NULL), fLogicTarget_5(NULL),
      fLogicTarget_6(NULL), fLogicTarget_7(NULL),
      fLogicFilm(NULL), fLogicChamber(NULL), fLogicSphere(NULL),
      fLogicAbsorber(NULL), fLogicAbsorber_4(NULL), fLogicAbsorber_5(NULL), fLogicAbsorber_6(NULL), fLogicAbsorber_7(NULL), // byHN
      fTargetMaterial(NULL), fTargetMaterial_GAGG(NULL), fFilmMaterial(NULL), fChamberMaterial(NULL), fSphereMaterial(NULL),
      fAbsorberMaterial(NULL), // byHN
      fStepLimit(NULL),
      fCheckOverlaps(true)
{
  fMessenger = new B2aDetectorMessenger(this);

  fNbOfChambers = 5;
  fLogicChamber = new G4LogicalVolume *[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2aDetectorConstruction::~B2aDetectorConstruction()
{
  delete[] fLogicChamber;
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B2aDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::DefineMaterials()
{
  /*********************************************************/
  constexpr G4double roomTemp = 293.15 * kelvin;
  /*********************************************************/

  // Material definition

  G4NistManager *nist = G4NistManager::Instance();

  // definition using NIST Manager
  nist->FindOrBuildMaterial("G4_AIR");
  nist->FindOrBuildMaterial("G4_Galactic");                // yk
  nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); // yk: G4_SODIUM_IODIDE, G4_PLASTIC_SC_VINYLTOLUENE
  nist->FindOrBuildMaterial("G4_PLEXIGLASS");              // yk: G4_PLEXIGLASS (acryl)

  G4bool isotopes = false;
  G4Element *H = nist->FindOrBuildElement("H", isotopes);   // Z =  1
  G4Element *He = nist->FindOrBuildElement("He", isotopes); // Z =  2
  G4Element *Li = nist->FindOrBuildElement("Li", isotopes); // Z =  3
  G4Element *C = nist->FindOrBuildElement("C", isotopes);   // Z =  6
  G4Element *N = nist->FindOrBuildElement("N", isotopes);   // Z =  7
  G4Element *O = nist->FindOrBuildElement("O", isotopes);   // Z =  8
  G4Element *Na = nist->FindOrBuildElement("Na", isotopes); // Z = 11
  G4Element *Cl = nist->FindOrBuildElement("Cl", isotopes); // Z = 17
  G4Element *Br = nist->FindOrBuildElement("Br", isotopes); // Z = 35
  G4Element *Y = nist->FindOrBuildElement("Y", isotopes);   // Z = 39
  G4Element *I = nist->FindOrBuildElement("I", isotopes);   // Z = 53
  G4Element *Cs = nist->FindOrBuildElement("Cs", isotopes); // Z = 55
  G4Element *La = nist->FindOrBuildElement("La", isotopes); // Z = 57
  G4Element *Tl = nist->FindOrBuildElement("Tl", isotopes); // Z = 81
  G4Element *Gd = nist->FindOrBuildElement("Gd", isotopes); // Z = 64
  G4Element *Cd = nist->FindOrBuildElement("Cd", isotopes); // Z = 48
  G4Element *Si = nist->FindOrBuildElement("Si", isotopes); // Z = 14
  G4Element *Mg = nist->FindOrBuildElement("Mg", isotopes); // Z = 12
  G4Element *Al = nist->FindOrBuildElement("Al", isotopes); // Z = 13
  G4Element *Ce = nist->FindOrBuildElement("Ce", isotopes); // Z = 58
  G4Element *Ga = nist->FindOrBuildElement("Ga", isotopes); // Z = 31
  G4Element *B = nist->FindOrBuildElement("B");             // Z = 5
  G4Element *Fe = nist->FindOrBuildElement("Fe");           // Z = 26

  // Li in EJ270: 95.5 atomic percent of lithium-6
  G4Isotope *Li6 = new G4Isotope("Li6", 3, 6, 6.02 * g / mole);
  G4Isotope *Li7 = new G4Isotope("Li7", 3, 7, 7.02 * g / mole);
  G4Element *IsoLi = new G4Element("isoLi", "Li", 2);
  IsoLi->AddIsotope(Li6, 95.5 * perCent);
  IsoLi->AddIsotope(Li7, 4.5 * perCent);

  // EJ-270 properties https://eljentechnology.com/products/plastic-scintillators/ej-270 2025/12/10
  double EJ270_H_frac = 5.20 * 1;
  double EJ270_C_frac = 4.62 * 12;
  double EJ270_Li_frac = 0.0575 * 6;
  double EJ270_O_frac = 0.479 * 16;
  double EJ270_N_frac = 0.0931 * 14;
  double EJ270_tot_frac = EJ270_H_frac + EJ270_C_frac + EJ270_Li_frac + EJ270_O_frac + EJ270_N_frac;
  EJ270 = new G4Material("G4_EJ270", 1.14 * g / cm3, 5);
  EJ270->AddElement(H, EJ270_H_frac / EJ270_tot_frac * 100 * perCent);
  EJ270->AddElement(C, EJ270_C_frac / EJ270_tot_frac * 100 * perCent);
  EJ270->AddElement(IsoLi, EJ270_Li_frac / EJ270_tot_frac * 100 * perCent);
  EJ270->AddElement(O, EJ270_O_frac / EJ270_tot_frac * 100 * perCent);
  EJ270->AddElement(N, EJ270_N_frac / EJ270_tot_frac * 100 * perCent);

  // B4C(Mirrobor)
  std::vector<std::pair<G4Element *, G4double>> vEl_Mirrobor{
      {H, 0.01},
      {B, 0.6321},
      {C, 0.3061},
      {O, 0.05},
      {Fe, 0.0018}};
  Mirrobor = new G4Material("Mirrobor", 1.36 * g / cm3, vEl_Mirrobor.size(), kStateSolid, roomTemp);
  for (const auto &el : vEl_Mirrobor)
    Mirrobor->AddElement(el.first, el.second);

  //Gdsheet (ニュートロン・ストップSY)
  std::vector<std::pair<G4Element *, G4double>> vEl_Gdsheet{
      {H, 0.049},
      {C, 0.194},
      {O, 0.183},
      {Si, 0.227},
      {Gd, 0.347}};
  Gdsheet = new G4Material("Gdsheet", 1.50 * g / cm3, vEl_Gdsheet.size(), kStateSolid, roomTemp);
  for (const auto &el : vEl_Gdsheet)
    Gdsheet->AddElement(el.first, el.second); 

  // polyethylene
  Polyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

  // G10(FR4)
  std::vector<std::pair<G4Element *, G4double>> vEl_G10{
      {H, 0.030500},
      {C, 0.363446},
      {O, 0.322763},
      {Si, 0.283291}};
  G10 = new G4Material("G10", 1.7 * g / cm3, vEl_G10.size());
  for (const auto &el : vEl_G10)
    G10->AddElement(el.first, el.second);

  // PTFE
  PTFE = nist->FindOrBuildMaterial("G4_TEFLON");

  // Cu
  Cu = nist->FindOrBuildMaterial("G4_Cu");

  // Air
  Air = nist->FindOrBuildMaterial("G4_AIR");

  // Pb Block
  PbBlock = nist->FindOrBuildMaterial("G4_Pb");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B2aDetectorConstruction::DefineVolumes()
{
  G4double worldLength = 10 * m; // yk: 1.2 * (2*targetLength + trackerLength)=588cm

  /*world*/
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / mm
         << " mm" << G4endl;

  G4Box *worldS = new G4Box("world",                                            // its name
                            worldLength / 2, worldLength / 2, worldLength / 2); // its size
  G4LogicalVolume *worldLV = new G4LogicalVolume(
      worldS,         // its solid
      Air, // its material
      "World");       // its name

  worldPV = new G4PVPlacement(
      0,               // no rotation
      G4ThreeVector(), // at (0,0,0)
      worldLV,         // its logical volume
      "World",         // its name
      0,               // its mother  volume
      false,           // no boolean operations
      0,               // copy number
      fCheckOverlaps); // checking overlaps

  ConstructSensor();
  // ConstructLab();

  // Visualization attributes
  // worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  return worldPV;
}

void B2aDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD *aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4cout << G4endl;
  G4cout << "----- Sensitive detectors -----" << G4endl;
  for (const auto &el : fSensitiveLogVolName)
  {
    SetSensitiveDetector(el, aTrackerSD, true);
  }
  G4cout << "-------------------------------" << G4endl;
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
}

void B2aDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager *nistManager = G4NistManager::Instance();

  G4Material *pttoMaterial =
      nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial_GAGG != pttoMaterial)
  {
    if (pttoMaterial)
    {
      fTargetMaterial_GAGG = pttoMaterial;
      if (fLogicTarget)
        fLogicTarget->SetMaterial(fTargetMaterial_GAGG);
      G4cout
          << G4endl
          << "----> The target is made of " << materialName << G4endl;
    }
    else
    {
      G4cout
          << G4endl
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
    }
  }
}

void B2aDetectorConstruction::SetTargetMaterial_0(G4String materialName)
{
  G4NistManager *nistManager = G4NistManager::Instance();

  G4Material *pttoMaterial =
      nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial)
  {
    if (pttoMaterial)
    {
      fTargetMaterial = pttoMaterial;
      if (fLogicTarget_0)
        fLogicTarget_0->SetMaterial(fTargetMaterial);
      G4cout
          << G4endl
          << "----> The target is made of " << materialName << G4endl;
    }
    else
    {
      G4cout
          << G4endl
          << "-->  WARNING from SetTargetMaterial_0 : "
          << materialName << " not found" << G4endl;
    }
  }
}

void B2aDetectorConstruction::SetTargetMaterial_1(G4String materialName)
{
  G4NistManager *nistManager = G4NistManager::Instance();

  G4Material *pttoMaterial =
      nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial)
  {
    if (pttoMaterial)
    {
      fTargetMaterial = pttoMaterial;
      if (fLogicTarget_1)
        fLogicTarget_1->SetMaterial(fTargetMaterial);
      G4cout
          << G4endl
          << "----> The target is made of " << materialName << G4endl;
    }
    else
    {
      G4cout
          << G4endl
          << "-->  WARNING from SetTargetMaterial_1 : "
          << materialName << " not found" << G4endl;
    }
  }
}

// void B2aDetectorConstruction::SetTargetMaterial_2(G4String materialName)
// {
//   G4NistManager* nistManager = G4NistManager::Instance();

//   G4Material* pttoMaterial =
//               nistManager->FindOrBuildMaterial(materialName);

//   if (fTargetMaterial != pttoMaterial) {
//      if ( pttoMaterial ) {
//         fTargetMaterial = pttoMaterial;
//         if (fLogicTarget_2) fLogicTarget_2->SetMaterial(fTargetMaterial);
//         G4cout
//           << G4endl
//           << "----> The target is made of " << materialName << G4endl;
//      } else {
//         G4cout
//           << G4endl
//           << "-->  WARNING from SetTargetMaterial_2 : "
//           << materialName << " not found" << G4endl;
//      }
//   }
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetFilmMaterial(G4String materialName)
{
  G4NistManager *nistManager = G4NistManager::Instance();

  G4Material *pttoMaterial =
      nistManager->FindOrBuildMaterial(materialName);

  if (fFilmMaterial != pttoMaterial)
  {
    if (pttoMaterial)
    {
      fFilmMaterial = pttoMaterial;
      if (fLogicFilm)
        fLogicFilm->SetMaterial(fFilmMaterial);
      G4cout
          << G4endl
          << "----> The film is made of " << materialName << G4endl;
    }
    else
    {
      G4cout
          << G4endl
          << "-->  WARNING from SetFilmMaterial : "
          << materialName << " not found" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager *nistManager = G4NistManager::Instance();

  G4Material *pttoMaterial =
      nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial)
  {
    if (pttoMaterial)
    {
      fChamberMaterial = pttoMaterial;
      for (G4int copyNo = 0; copyNo < fNbOfChambers; copyNo++)
      {
        if (fLogicChamber[copyNo])
          fLogicChamber[copyNo]->SetMaterial(fChamberMaterial);
      }
      G4cout
          << G4endl
          << "----> The chambers are made of " << materialName << G4endl;
    }
    else
    {
      G4cout
          << G4endl
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit) && (maxStep > 0.))
    fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}
