#ifndef Geometry_h
#define Geometry_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class GeometryMessenger;

class Geometry : public G4VUserDetectorConstruction
{
public:
	Geometry();
	~Geometry() override;

	G4VPhysicalVolume* Construct() override;
	void ConstructSDandField() override;

	void PrintParameters();
	void SetLunarSoilHfrac(G4double val){ fLunarSoilHfrac = val; };
	void SetLanderFuelMass(G4double val){ fLanderFuelMass = val; };

	void AddSensitiveDetector(G4String name){ fSensitiveLogVolName.emplace_back(name); };

private:  
	G4VPhysicalVolume* ConstructVolume();
	void ConstructMaterials();
	void ConstructMoon();
	void ConstructLander();
	void ConstructScorer();

	std::vector<G4String> fSensitiveLogVolName;

	//G4LogicalVolume* lv_world;
	G4LogicalVolume* lv_periodic_world;

	G4Material* Vacuum;
	//G4Material* Air;
	//G4Material* Helium3_10atm;
	//G4Material* Water;
	//G4Material* Mirrobor;
	G4Material* FHT;
	G4Material* Lander;
	G4Material* Fuel;
	G4Material* Oxidizer;
	G4Material* EJ270;

	//G4Material* TS_Al;
	//G4Material* TS_Fe;
	//G4Material* TS_Water;
	//G4Material* TS_A5052;

	GeometryMessenger* fMessenger = nullptr;
	G4double fLunarSoilHfrac = 0;
	G4double fLanderFuelMass = 0;
};

struct tubeParam{
	G4double rb; //radius begin
	G4double re; //radius end
	G4double lz; //length
	G4double pb; //phi begin
	G4double pe; //phi end
};

#endif

