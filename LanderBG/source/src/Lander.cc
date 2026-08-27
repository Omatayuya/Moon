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

/* Lunar lander */
void Geometry::ConstructLander()
{
	/*********************************************************/
	G4int copyNb = 2000;
	/*********************************************************/

	/* lunar lander */
	if(true){
		{ //outer
			auto sl = new G4Tubs("SL_Lander_Outer", 0.5*2*m, 0.5*3*m, 0.5*0.7*m, 0*deg, 360*deg);
			auto lv = new G4LogicalVolume(sl, Lander, "LV_Lander_Outer");
			lv->SetVisAttributes(G4Color::Cyan());
			//AddSensitiveDetector(lv->GetName());

			auto tf = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 1.4*m+0.5*0.7*m));
			new G4PVPlacement(tf, lv, "PV_Lander_Outer", lv_periodic_world, false, ++copyNb, true);
		}

		{ //inner
			auto sl = new G4Tubs("SL_Lander_Inner", 0, 0.5*2*m, 0.5*1.5*m, 0*deg, 360*deg);
			auto lv = new G4LogicalVolume(sl, Lander, "LV_Lander_Inner");
			lv->SetVisAttributes(G4Color::Cyan());
			//AddSensitiveDetector(lv->GetName());

			auto tf = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 0.6*m+0.5*1.5*m));
			new G4PVPlacement(tf, lv, "PV_Lander_Inner", lv_periodic_world, false, ++copyNb, true);
		
			{ //oxidizer, fuel
				auto sl_of = new G4Tubs("SL_Lander_OF", 0, 0.5*0.8*m, 0.5*1*m, 0*deg, 360*deg);

				//fuel
				auto lv_fuel = new G4LogicalVolume(sl_of, Fuel, "LV_Lander_Fuel");
				lv_fuel->SetVisAttributes(G4Color::Red());
				//AddSensitiveDetector(lv_fuel->GetName());

				auto tf_fuel1 = G4Transform3D(G4RotationMatrix(), G4ThreeVector(-0.5*0.8*m, 0.5*0.8*m, 0));
				new G4PVPlacement(tf_fuel1, lv_fuel, "PV_Lander_Fuel-1", lv, false, ++copyNb, true);

				auto tf_fuel2 = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0.5*0.8*m, -0.5*0.8*m, 0));
				new G4PVPlacement(tf_fuel2, lv_fuel, "PV_Lander_Fuel-2", lv, false, ++copyNb, true);

				//oxidizer
				auto lv_oxidizer = new G4LogicalVolume(sl_of, Oxidizer, "LV_Lander_Oxidizer");
				lv_oxidizer->SetVisAttributes(G4Color::Blue());
				//AddSensitiveDetector(lv_oxidizer->GetName());

				auto tf_oxidizer1 = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0.5*0.8*m, 0.5*0.8*m, 0));
				new G4PVPlacement(tf_oxidizer1, lv_oxidizer, "PV_Lander_Oxidizer-1", lv, false, ++copyNb, true);

				auto tf_oxidizer2 = G4Transform3D(G4RotationMatrix(), G4ThreeVector(-0.5*0.8*m, -0.5*0.8*m, 0));
				new G4PVPlacement(tf_oxidizer2, lv_oxidizer, "PV_Lander_Oxidizer-2", lv, false, ++copyNb, true);
			}
		}
	}
}

