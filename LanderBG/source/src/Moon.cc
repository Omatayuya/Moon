#include "Geometry.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SubtractionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

/* Lunar soil */
void Geometry::ConstructMoon()
{
	/*********************************************************/
	G4int copyNb = 1000;
	/*********************************************************/

	/* lunar soil */
	{
		auto sl = new G4Box("SL_Moon", 0.5*30*m, 0.5*30*m, 0.5*10*m);
		auto lv = new G4LogicalVolume(sl, FHT, "LV_Moon");
		lv->SetVisAttributes(G4Color::Brown());
		//AddSensitiveDetector(lv->GetName());

		auto tf = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, -0.5*10*m));
		new G4PVPlacement(tf, lv, "PV_Moon", lv_periodic_world, false, ++copyNb, true);
	}

}

