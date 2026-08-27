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
#include "G4SubtractionSolid.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

#include "CADMesh.hh"

void B2aDetectorConstruction::ConstructLab()
{
  G4int copyNb = 3000;

  /*virtual scoring detector*/
  constexpr G4double scorerX = 7 * cm;
  constexpr G4double scorerY = 7 * cm;
  constexpr G4double scorerZ = 10 * um;
  constexpr G4double scorerOffset = 10 * cm;
  {
    auto box = new G4Box("SL_scorer2", scorerX / 2, scorerY / 2, scorerZ / 2);
    auto lv = new G4LogicalVolume(box, Air, "LV_scorer2");
    lv->SetVisAttributes(G4VisAttributes(G4Colour::White()));
    AddSensitiveDetector(lv->GetName());

    auto pos = G4ThreeVector(0, 0, scorerOffset);
    new G4PVPlacement(0, pos, lv, "PV_scorer2", worldPV->GetLogicalVolume(), false, ++copyNb, fCheckOverlaps);
  }

  /*Moderator*/

  
}