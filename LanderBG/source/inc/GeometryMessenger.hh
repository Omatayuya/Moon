#ifndef GeometryMessenger_h
#define GeometryMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Geometry;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

class GeometryMessenger: public G4UImessenger
{
public:
	GeometryMessenger(Geometry* );
	~GeometryMessenger() override;

	void SetNewValue(G4UIcommand*, G4String) override;

private:
    Geometry* fGeometry = nullptr;
	G4UIdirectory* fDir = nullptr;

	G4UIcmdWithADouble* fLunarSoilHfracCmd = nullptr;
	G4UIcmdWithADoubleAndUnit* fLanderFuelMassCmd = nullptr;
};

#endif

