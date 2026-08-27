#include "B2ParallelWorld.hh"
#include "B2aDetectorConstruction.hh"
#include "B2TrackerSD.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"

#include <vector>

MyParallelWorld::MyParallelWorld(G4String worldName, const B2aDetectorConstruction *detector)
    : G4VUserParallelWorld(worldName), fDetector(detector) {}

void MyParallelWorld::Construct()
{
    if (!fDetector)
    {
        G4Exception("MyParallelWorld::Construct", "B2aParallelWorld001",
                    FatalException, "B2aDetectorConstruction pointer is null.");
    }

    // 1. 並行世界の一番外側の枠（ゴーストワールド）を取得する
    G4VPhysicalVolume *ghostWorld = GetWorld();
    G4LogicalVolume *ghostLogical = ghostWorld->GetLogicalVolume();

    const auto &SensorCenter = fDetector->GetSensorCenter();
    const auto &SensorHalfSize = fDetector->GetSensorHalfSize();
    const auto GdsheetT = fDetector->GetGdsheetT();

    constexpr G4double scorerT = 10 * um;
    constexpr G4double scorerbackhX = scorerT;
    constexpr G4double scorerbackhY = 7 * cm;
    constexpr G4double scorerbackhZ = 7 * cm;

    G4double hX = SensorHalfSize.x();
    G4double hY = SensorHalfSize.y();
    G4double hZ = SensorHalfSize.z();
    G4double hT = scorerT / 2.0;

    struct ScorerParam
    {
        G4String name;
        G4ThreeVector halfSize;
        G4ThreeVector offset;
    };

    std::vector<ScorerParam> vScorer = {
        {"Back", G4ThreeVector(scorerbackhX / 2, scorerbackhY / 2, scorerbackhZ / 2), G4ThreeVector(-(hX + hT), 0, 0)},
        {"Right", G4ThreeVector(hX, hT, hZ), G4ThreeVector(0, hY, 0)},
        {"Left", G4ThreeVector(hX, hT, hZ), G4ThreeVector(0, -hY, 0)},
        {"Top", G4ThreeVector(hX, hY, hT), G4ThreeVector(0, 0, hZ)},
        {"Bottom", G4ThreeVector(hX, hY, hT), G4ThreeVector(0, 0, -hZ)},
        {"Front", G4ThreeVector(hT, hY, hZ), G4ThreeVector(hX, 0, 0)},
    };

    G4int copyNb_scorer = 1000;

    for (const auto &param : vScorer)
    {
        auto box = new G4Box("Box_Scorer" + param.name,
                             param.halfSize.x(),
                             param.halfSize.y(),
                             param.halfSize.z());

        auto lv = new G4LogicalVolume(box, nullptr, "LV_Scorer" + param.name);
        lv->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

        G4ThreeVector pos = SensorCenter + param.offset;
        new G4PVPlacement(nullptr, pos, lv, "PV_Scorer" + param.name,
                          ghostLogical, false, ++copyNb_scorer, true);
    }
}

void MyParallelWorld::ConstructSD()
{
    G4SDManager *sdMan = G4SDManager::GetSDMpointer();

    auto scorerSD = new B2TrackerSD("ParallelWorld/ScorerSD",
                                    "ParallelWorldHitsCollection");
    sdMan->AddNewDetector(scorerSD);

    const std::vector<G4String> scorerNames = {"Back", "Top", "Bottom", "Right", "Left", "Front"};
    for (const auto &name : scorerNames)
    {
        SetSensitiveDetector("LV_Scorer" + name, scorerSD,true);
    }
}
