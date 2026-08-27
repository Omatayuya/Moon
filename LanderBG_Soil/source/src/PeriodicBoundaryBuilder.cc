#include "PeriodicBoundaryBuilder.hh"

#include "LogicalVolumePeriodic.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

G4LogicalVolume* PeriodicBoundaryBuilder::Construct(G4LogicalVolume* logical_world)
{
	auto world = dynamic_cast<G4Box*>(logical_world->GetSolid());

	if(world==nullptr){
		G4ExceptionDescription ed;
		ed << " PeriodicBoundaryBuilder::Construct: Unsupported period boundary for this solid : " << logical_world->GetName() << G4endl;
		G4Exception("G4PeriodicBoundaryProcess::G4PeriodicBoundaryProcess", "Builder01", FatalException, ed, "Unsupported period boundary for this solid");
	}

	constexpr G4double buffer = 1.*um;

	G4double periodic_world_hx = world->GetXHalfLength();
	G4double periodic_world_hy = world->GetYHalfLength();
	G4double periodic_world_hz = world->GetZHalfLength();

	world->SetXHalfLength(world->GetXHalfLength()+0.5*buffer);
	world->SetYHalfLength(world->GetYHalfLength()+0.5*buffer);
	world->SetZHalfLength(world->GetZHalfLength()+0.5*buffer);

	auto periodic_world = new G4Box("PBC", periodic_world_hx, periodic_world_hy, periodic_world_hz);
	fLogicalPeriodic = new LogicalVolumePeriodic(periodic_world, logical_world->GetMaterial(), "PBC");
	fLogicalPeriodic->SetVisAttributes(G4Color::Magenta());
	new G4PVPlacement(nullptr, G4ThreeVector(), fLogicalPeriodic, "PBC", logical_world, false, 0, true);

	return fLogicalPeriodic;
}

