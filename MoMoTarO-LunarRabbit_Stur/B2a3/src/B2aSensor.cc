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
	G4String cadModelPath = "/home/yomata/work/Moon/MoMoTarO-LunarRabbit_Stur/Sensor_CAD_ver20260906/";
	
	std::vector<CadParam> vCadVolume{// MoMoTarO_B_SoilExp202412_v2
									 {"OuterBox_bottom.stl", TS_A5052, G4Color::Gray()},
									 {"OuterBox_side_1.stl", TS_A5052, G4Color::Gray()},
									 {"OuterBox_side_2.stl", TS_A5052, G4Color::Gray()},
									 {"OuterBox_side_3.stl", TS_A5052, G4Color::Gray()},
									 {"OuterBox_side_4.stl", TS_A5052, G4Color::Gray()},
									 {"EJ-270_1.stl", EJ270, G4Color::Blue(), true},
									 {"EJ-270_1_MPPC.stl", G10, G4Color::Green()},
									 {"EJ-270_1_Reflector.stl", PTFE, G4Color::Cyan()},
									 {"EJ-270_1_Box_top.stl", TS_A5052, G4Color::Gray()},
									 {"EJ-270_1_Box_side.stl", TS_A5052, G4Color::Gray()},
									 {"EJ-270_1_Box_bottom.stl", TS_A5052, G4Color::Gray()},
									 {"EJ-270_2-1.stl", EJ270, G4Color::Blue(), true},
									 {"EJ-270_2-1_MPPC.stl", G10, G4Color::Green()},
									 {"EJ-270_2-1_Reflector.stl", PTFE, G4Color::Cyan()},
									 {"EJ-270_2-2.stl", EJ270, G4Color::Blue(), true},
									 {"EJ-270_2-2_MPPC.stl", G10, G4Color::Green()},
									 {"EJ-270_2-2_Reflector.stl", PTFE, G4Color::Cyan()},
									 {"EJ-270_2-3.stl", EJ270, G4Color::Blue(), true},
									 {"EJ-270_2-3_MPPC.stl", G10, G4Color::Green()},
									 {"EJ-270_2-3_Reflector.stl", PTFE, G4Color::Cyan()},
									 {"EJ-270_2_Cd.stl", Cd, G4Color::Red()},
									 {"EJ-270_2_Box_top.stl", TS_A5052, G4Color::Gray()},
									 {"EJ-270_2_Box_side.stl", TS_A5052, G4Color::Gray()},
									 {"EJ-270_2_Box_bottom.stl", TS_A5052, G4Color::Gray()},
									 {"GAGG.stl", GAGG, G4Color::Yellow(), true},
									 {"GAGG_MPPC.stl", G10, G4Color::Green()},
									 {"GAGG_Reflector.stl", PTFE, G4Color::Cyan()},
									 {"GAGG_Box_top.stl", TS_A5052, G4Color::Gray()},
									 {"GAGG_Box_side.stl", TS_A5052, G4Color::Gray()},
									 {"GAGG_Box_bottom.stl", TS_A5052, G4Color::Gray()},
									 {"DAQ_board.stl", G10, G4Color::Green()},
									 {"Bolt.stl", TS_SUS304, G4Color::Brown()},
									 {"Spacer_Box.stl", TS_A5052, G4Color::Gray()},
									 {"Spacer_CircuitBox.stl", TS_A5052, G4Color::Gray()}};

	constexpr G4double scale = 1;
	constexpr G4double rotX = 90 * deg;
	constexpr G4double offsetZ = 2.51 * mm;

	G4int copyNb = 2000;
	/**************************************************/

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

		auto tf = G4Transform3D(rot, G4ThreeVector(0, 0, offsetZ));
		new G4PVPlacement(tf, "PV_" + basename, lv, worldPV, false, ++copyNb, true);
	}
}