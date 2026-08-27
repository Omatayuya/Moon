#ifndef LogicalVolumePeriodic_hh
#define LogicalVolumePeriodic_hh 1

#include "G4LogicalVolume.hh"

class LogicalVolumePeriodic : public G4LogicalVolume
{
public:
	LogicalVolumePeriodic(G4VSolid* pSolid, G4Material* pMaterial, const G4String& name)
		: G4LogicalVolume(pSolid, pMaterial, name){};
	~LogicalVolumePeriodic() override = default;

	G4bool IsExtended() const override { return true; }
};
#endif

