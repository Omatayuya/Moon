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
#include "G4SubtractionSolid.hh"

#include "CADMesh.hh"

void B2aDetectorConstruction::ConstructSensor()
{
  /*EJ-270 cube (10*10*20 cm^3)*/
  G4int copyNb = 0;
  {
    constexpr G4double SizeX = 7 * cm;
    constexpr G4double SizeY = 7 * cm;
    constexpr G4double SizeZ = 80 * cm;

    constexpr G4double offsetZ = 50 * cm;

    auto box = new G4Box("Box_EJ270cube", SizeX / 2, SizeY / 2, SizeZ / 2);
    auto lv = new G4LogicalVolume(box, EJ270, "LV_EJ270cube");
    lv->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
    AddSensitiveDetector(lv->GetName());

    auto tf = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, offsetZ));
    new G4PVPlacement(tf, "PV_EJ270cube", lv, worldPV, false, ++copyNb, true);

    fSensorCenter = G4ThreeVector(0, 0, offsetZ);
    fSensorHalfSize = G4ThreeVector(SizeX / 2, SizeY / 2, SizeZ / 2);
  }

  /*socrer*/
  // G4int copyNb_scorer = 1000;
  // constexpr G4double scorerT = 10 * um;

  // constexpr G4double scorerbackhX = scorerT;
  // constexpr G4double scorerbackhY = 7 * cm;
  // constexpr G4double scorerbackhZ = 7 * cm;

  // G4double hX = DetHalfSize.x();
  // G4double hY = DetHalfSize.y();
  // G4double hZ = DetHalfSize.z();
  // G4double hT = scorerT / 2.0;

  // struct ScorerParam
  // {
  //   G4String name;
  //   G4ThreeVector halfSize;
  //   G4ThreeVector offset;
  //   G4Colour color;
  // };

  // std::vector<ScorerParam> vScorer = {
  //     {"Back", G4ThreeVector(scorerbackhX / 2, scorerbackhY / 2, scorerbackhZ / 2), G4ThreeVector(-(hX + hT), 0, 0), G4Colour::Red()},
  //     {"Top", G4ThreeVector(hX, hT, hZ), G4ThreeVector(0, hY + hT + GdsheetT, 0), G4Colour::Red()},
  //     {"Bot", G4ThreeVector(hX, hT, hZ), G4ThreeVector(0, -(hY + hT + GdsheetT), 0), G4Colour::Red()},
  //     {"Right", G4ThreeVector(hX, hY, hT), G4ThreeVector(0, 0, hZ + hT + GdsheetT), G4Colour::Red()},
  //     {"Left", G4ThreeVector(hX, hY, hT), G4ThreeVector(0, 0, -(hZ + hT + GdsheetT)), G4Colour::Red()},
  //     {"Front", G4ThreeVector(hT, hY, hZ), G4ThreeVector(hX + hT + GdsheetT, 0, 0), G4Colour::Red()},
  // };

  // for (const auto &param : vScorer)
  // {
  //   auto box = new G4Box("Box_Scorer" + param.name,
  //                        param.halfSize.x(),
  //                        param.halfSize.y(),
  //                        param.halfSize.z());

  //   auto lv = new G4LogicalVolume(box, Air, "LV_Scorer" + param.name);
  //   lv->SetVisAttributes(G4VisAttributes(param.color));
  //   AddSensitiveDetector(lv->GetName());

  //   G4ThreeVector pos = DetCenter + param.offset;
  //   new G4PVPlacement(0, pos, lv, "PV_Scorer" + param.name,
  //                     worldPV->GetLogicalVolume(), false, ++copyNb_scorer, true);
  // }
}