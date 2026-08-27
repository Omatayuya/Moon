#include "Geometry.hh"
#include "GeometryMessenger.hh"
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

/* virtual scoring detector */
void Geometry::ConstructScorer()
{
	/*********************************************************/
	constexpr G4double ej270X = 10*cm;
	constexpr G4double ej270Y = 10*cm;
	constexpr G4double ej270Z = 20*cm;

	constexpr G4double outerRingBottomZ = 1.4*m; //bottom face of LV_Lander_Outer
	constexpr G4double gap = 5*cm;
	constexpr G4double posX = 1.25*m;

	G4int copyNb = 3000;
	/*********************************************************/

	auto sl = new G4Box("SL_EJ270", 0.5*ej270X, 0.5*ej270Y, 0.5*ej270Z);
	auto lv = new G4LogicalVolume(sl, EJ270, "LV_EJ270");
	lv->SetVisAttributes(G4Color::Yellow());
	AddSensitiveDetector(lv->GetName());

	G4double posZ = outerRingBottomZ - gap - 0.5*ej270Z;
	auto tf = G4Transform3D(G4RotationMatrix(), G4ThreeVector(posX, 0, posZ));
	new G4PVPlacement(tf, lv, "PV_EJ270", lv_periodic_world, false, ++copyNb, true);
}

