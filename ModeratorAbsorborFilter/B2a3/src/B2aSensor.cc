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
  /*multi-layer based scintillator*/
  // G4String cadModelPath = "/home/yomata/work/Nspectrum/exp_202603/Sensor_CADmodel/";
  G4String cadModelPath = "/home/yomata/work/Nspectrum/exp_202603/Sensor_CADmodel/";

  std::vector<CadParam> vCadVolume = {
      {"EJ-270_1_5mmt.stl", EJ270, G4Color::Blue(), true},
      {"EJ-270_1_5mmt_MPPC.stl", G10, G4Color::Green()},
      {"EJ-270_1_5mmt_Reflector.stl", PTFE, G4VisAttributes::GetInvisible()},
      {"EJ-270_2_5mmt.stl", EJ270, G4Color::Blue(), true},
      {"EJ-270_2_5mmt_MPPC.stl", G10, G4Color::Green()},
      {"EJ-270_2_5mmt_Reflector.stl", PTFE, G4VisAttributes::GetInvisible()},
      {"EJ-270_3_10mmt.stl", EJ270, G4Color::Blue(), true},
      {"EJ-270_3_10mmt_MPPC.stl", G10, G4Color::Green()},
      {"EJ-270_3_10mmt_Reflector.stl", PTFE, G4VisAttributes::GetInvisible()},
      {"EJ-270_4_10mmt.stl", EJ270, G4Color::Blue(), true},
      {"EJ-270_4_10mmt_MPPC.stl", G10, G4Color::Green()},
      {"EJ-270_4_10mmt_Reflector.stl", PTFE, G4VisAttributes::GetInvisible()},
      {"EJ-270_5_10mmt.stl", EJ270, G4Color::Blue(), true},
      {"EJ-270_5_10mmt_MPPC.stl", G10, G4Color::Green()},
      {"EJ-270_5_10mmt_Reflector.stl", PTFE, G4VisAttributes::GetInvisible()},
      {"EJ-270_6_10mmt.stl", EJ270, G4Color::Blue(), true},
      {"EJ-270_6_10mmt_MPPC.stl", G10, G4Color::Green()},
      {"EJ-270_6_10mmt_Reflector.stl", PTFE, G4VisAttributes::GetInvisible()},
      {"EJ-270_7_10mmt.stl", EJ270, G4Color::Blue(), true},
      {"EJ-270_7_10mmt_MPPC.stl", G10, G4Color::Green()},
      {"EJ-270_7_10mmt_Reflector.stl", PTFE, G4VisAttributes::GetInvisible()},
      {"EJ-270_8_30mmt.stl", EJ270, G4Color::Blue(), true},
      {"EJ-270_8_30mmt_MPPC.stl", G10, G4Color::Green()},
      {"EJ-270_8_30mmt_Reflector.stl", PTFE, G4VisAttributes::GetInvisible()},
  };

  constexpr G4double scale = 1;
  constexpr G4double rotX = 90 * deg;
  constexpr G4double rotY = 0 * deg;
  constexpr G4double rotZ = -90 * deg;
  constexpr G4double offsetX = 0 * mm;

  G4int copyNb = 2000;
  /**************************************************/

  G4ThreeVector detMin(kInfinity, kInfinity, kInfinity);
  G4ThreeVector detMax(-kInfinity, -kInfinity, -kInfinity);
  G4ThreeVector sensMin(kInfinity, kInfinity, kInfinity);
  G4ThreeVector sensMax(-kInfinity, -kInfinity, -kInfinity);

  // CADMesh::STL
  for (const auto &el : vCadVolume)
  {
    G4String basename = el.filename.substr(0, el.filename.rfind('.'));

    auto mesh = CADMesh::TessellatedMesh::FromSTL(cadModelPath + el.filename);
    mesh->SetScale(scale);
    mesh->SetOffset(0, 0, 0);

    auto sl = mesh->GetSolid();
    sl->SetName("SL_" + basename);

    auto lv = new G4LogicalVolume(sl, el.material, "LV_" + basename);
    lv->SetVisAttributes(el.vis);
    if (el.sensVolFlag)
      AddSensitiveDetector(lv->GetName());

    auto rot = G4RotationMatrix();
    rot.rotateX(rotX);
    rot.rotateY(rotY);
    rot.rotateZ(rotZ);
    auto tf = G4Transform3D(rot, G4ThreeVector(offsetX, 0, 0));

    G4ThreeVector pMin, pMax;
    sl->BoundingLimits(pMin, pMax); // STLのローカルな最小・最大座標を取得

    G4Point3D corners[8] = {
        G4Point3D(pMin.x(), pMin.y(), pMin.z()),
        G4Point3D(pMax.x(), pMin.y(), pMin.z()),
        G4Point3D(pMin.x(), pMax.y(), pMin.z()),
        G4Point3D(pMax.x(), pMax.y(), pMin.z()),
        G4Point3D(pMin.x(), pMin.y(), pMax.z()),
        G4Point3D(pMax.x(), pMin.y(), pMax.z()),
        G4Point3D(pMin.x(), pMax.y(), pMax.z()),
        G4Point3D(pMax.x(), pMax.y(), pMax.z())};

    for (int i = 0; i < 8; ++i)
    {
      G4Point3D globalPt = tf * corners[i]; // グローバル座標系に変換
      detMin.setX(std::min(detMin.x(), globalPt.x()));
      detMin.setY(std::min(detMin.y(), globalPt.y()));
      detMin.setZ(std::min(detMin.z(), globalPt.z()));
      detMax.setX(std::max(detMax.x(), globalPt.x()));
      detMax.setY(std::max(detMax.y(), globalPt.y()));
      detMax.setZ(std::max(detMax.z(), globalPt.z()));
    }
    if (el.sensVolFlag)
      for (int i = 0; i < 8; ++i)
      {
        G4Point3D globalPt = tf * corners[i]; // グローバル座標系に変換
        sensMin.setX(std::min(sensMin.x(), globalPt.x()));
        sensMin.setY(std::min(sensMin.y(), globalPt.y()));
        sensMin.setZ(std::min(sensMin.z(), globalPt.z()));
        sensMax.setX(std::max(sensMax.x(), globalPt.x()));
        sensMax.setY(std::max(sensMax.y(), globalPt.y()));
        sensMax.setZ(std::max(sensMax.z(), globalPt.z()));
      }

    new G4PVPlacement(tf, "PV_" + basename, lv, worldPV, false, ++copyNb, true);
  }

  G4ThreeVector DetCenter = (detMin + detMax) / 2;
  G4ThreeVector DetHalfSize = (detMax - detMin) / 2;
  G4ThreeVector SensCenter = (sensMin + sensMax) / 2;
  G4ThreeVector SensHalfSize = (sensMax - sensMin) / 2;

  /*B4C mirrorbor*/
  // G4ThreeVector DetCenter = (sensMax + sensMin) / 2;
  // G4ThreeVector DetHalfSize = (sensMax - sensMin) / 2;
  // constexpr G4double mirrorborT = 0.6 * cm;
  // constexpr G4double mirrorborWindowY = 3.2 * cm;
  // constexpr G4double mirrorborWindowZ = 3.2 * cm;
  // {
  //   auto box = new G4Box("Box_mirrorbor", DetHalfSize.x() + mirrorborT, DetHalfSize.y() + mirrorborT, DetHalfSize.z() + mirrorborT);
  //   auto hole = new G4Box("Hole_mirrorbor", DetHalfSize.x(), DetHalfSize.y(), DetHalfSize.z());
  //   auto sl_1 = new G4SubtractionSolid("Subt_mirrorbor_1", box, hole);

  //   auto widow = new G4Box("Window_mirrorbor", mirrorborT, mirrorborWindowY / 2, mirrorborWindowZ / 2);
  //   G4ThreeVector windowPos = G4ThreeVector(-(DetHalfSize.x() + mirrorborT / 2), 0, 0);

  //   auto sl = new G4SubtractionSolid("Subt_mirrorbor", sl_1, widow, 0, windowPos);
  //   auto lv = new G4LogicalVolume(sl, Mirrobor, "LV_mirrorbor");
  //   lv->SetVisAttributes(G4VisAttributes(G4Colour::Gray()));

  //   auto tf = G4Transform3D(G4RotationMatrix(), DetCenter);
  //   new G4PVPlacement(tf, "PV_mirrorbor", lv, worldPV, false, ++copyNb, true);
  // }

  /*Gdシート(ニュートロン・ストップ SY)*/
  constexpr G4double GdsheetT = 0.5 * cm;
  constexpr G4double GdsheetWindowY = 7 * cm;
  constexpr G4double GdsheetWindowZ = 7 * cm;
  {
    auto box = new G4Box("Box_Gdsheet", DetHalfSize.x() + GdsheetT / 2, DetHalfSize.y() + GdsheetT / 2, DetHalfSize.z() + GdsheetT / 2);
    auto hole = new G4Box("Hole_Gdsheet", DetHalfSize.x(), DetHalfSize.y(), DetHalfSize.z());
    auto sl_1 = new G4SubtractionSolid("Subt_Gdsheet_1", box, hole);

    auto widow = new G4Box("Window_Gdsheet", GdsheetT, GdsheetWindowY / 2, GdsheetWindowZ / 2);
    G4ThreeVector windowPos = G4ThreeVector(-(DetHalfSize.x() + GdsheetT / 2), 0, 0);

    auto sl = new G4SubtractionSolid("Subt_Gdsheet", sl_1, widow, 0, windowPos);
    auto lv = new G4LogicalVolume(sl, Gdsheet, "LV_Gdsheet");
    lv->SetVisAttributes(G4VisAttributes(G4Colour::White()));

    auto tf = G4Transform3D(G4RotationMatrix(), DetCenter);
    new G4PVPlacement(tf, "PV_Gdsheet", lv, worldPV, false, ++copyNb, true);
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

  // fSensorCenter = SensCenter;
  // fSensorHalfSize = SensHalfSize;
  // fGdsheetT = GdsheetT;
}