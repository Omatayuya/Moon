#include "GeometryMessenger.hh"
#include "Geometry.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

GeometryMessenger::GeometryMessenger(Geometry* geom)
 : fGeometry(geom)
{
	fDir = new G4UIdirectory("/myGeom/");
	fDir->SetGuidance("UI commands specific to target station");

	fLunarSoilHfracCmd = new G4UIcmdWithADouble("/myGeom/setLunarSoilHfrac", this);
	fLunarSoilHfracCmd->SetGuidance("Define H fraction in lunar soil");
	fLunarSoilHfracCmd->SetParameterName("lunarSoilHfrac", true);
	fLunarSoilHfracCmd->SetDefaultValue(0);
	fLunarSoilHfracCmd->AvailableForStates(G4State_PreInit);

	fLanderFuelMassCmd = new G4UIcmdWithADoubleAndUnit("/myGeom/setLanderFuelMass", this);
	fLanderFuelMassCmd->SetGuidance("Define mass of fuel & oxidizer in lunar lander");
	fLanderFuelMassCmd->SetParameterName("landerFuelMass", true);
	fLanderFuelMassCmd->SetDefaultValue(0);
	fLanderFuelMassCmd->SetDefaultUnit("kg");
	fLanderFuelMassCmd->SetUnitCategory("Mass");
	fLanderFuelMassCmd->AvailableForStates(G4State_PreInit);
}

GeometryMessenger::~GeometryMessenger()
{
	delete fLunarSoilHfracCmd;
	delete fLanderFuelMassCmd;
	delete fDir;
}

void GeometryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command==fLunarSoilHfracCmd){
		fGeometry->SetLunarSoilHfrac(fLunarSoilHfracCmd->GetNewDoubleValue(newValue));
	}else if(command==fLanderFuelMassCmd){
		fGeometry->SetLanderFuelMass(fLanderFuelMassCmd->GetNewDoubleValue(newValue));
	}
}

