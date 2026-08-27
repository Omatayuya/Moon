#include "Geometry.hh"
#include "GeometryMessenger.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4VisAttributes.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
// #include "G4GDMLParser.hh"

#include "SensitiveVolume.hh"
#include "G4SDManager.hh"

#include "PeriodicBoundaryBuilder.hh"

Geometry::Geometry()
	: G4VUserDetectorConstruction()
{
	fSensitiveLogVolName.clear();

	fMessenger = new GeometryMessenger(this);
}

Geometry::~Geometry()
{
	delete fMessenger;
}

G4VPhysicalVolume *Geometry::Construct()
{
	ConstructMaterials();
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	PrintParameters();

	return ConstructVolume();
}

G4VPhysicalVolume *Geometry::ConstructVolume()
{
	/* world */
	constexpr G4double buffer = 1 * um;
	constexpr G4double worldSizeX = 5 * m;
	constexpr G4double worldSizeY = 5 * m;
	constexpr G4double worldSizeZ = 6 * m;
	auto sl_world = new G4Box("SL_World", 0.5 * worldSizeX + buffer, 0.5 * worldSizeY + buffer, 0.5 * worldSizeZ + buffer);
	auto lv_world = new G4LogicalVolume(sl_world, Vacuum, "LV_World");
	lv_world->SetVisAttributes(G4VisAttributes::GetInvisible());
	// lv_world->SetVisAttributes(G4Color::Gray());
	auto tf_world = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 0));
	auto pv_world = new G4PVPlacement(tf_world, lv_world, "PV_World", nullptr, false, 0);

	/* periodic boundary volume */
	PeriodicBoundaryBuilder *pbb = new PeriodicBoundaryBuilder();
	lv_periodic_world = pbb->Construct(lv_world);

	// ConstructMoon();
	ConstructLander();
	ConstructScorer();

	return pv_world;
}

void Geometry::ConstructSDandField()
{
	auto aSensVol = new SensitiveVolume("SensitiveVolume", "HitsCollection");

	G4cout << G4endl;
	G4cout << "----- Sensitive detectors -----" << G4endl;
	for (const auto &el : fSensitiveLogVolName)
	{
		G4cout << el << G4endl;
		SetSensitiveDetector(el, aSensVol);
	}
	G4cout << "-------------------------------" << G4endl;

	auto SDManager = G4SDManager::GetSDMpointer();
	SDManager->AddNewDetector(aSensVol);
}

void Geometry::ConstructMaterials()
{
	/*********************************************************/
	// constexpr G4double roomTemp = 293.15*kelvin;
	constexpr G4double moonTemp = 175. * kelvin;
	constexpr G4double landerTemp = 293.15 * kelvin;
	/*********************************************************/

	auto nistManager = G4NistManager::Instance();

	/* isotope */
	// auto iso3He = new G4Isotope("He3", 2, 3, 3.0160293201*g/mole);

	/* element */
	auto elH = nistManager->FindOrBuildElement("H");
	// auto elLi = nistManager->FindOrBuildElement("Li");
	// auto elB = nistManager->FindOrBuildElement("B");
	auto elC = nistManager->FindOrBuildElement("C");
	auto elN = nistManager->FindOrBuildElement("N");
	auto elO = nistManager->FindOrBuildElement("O");
	auto elNa = nistManager->FindOrBuildElement("Na");
	auto elMg = nistManager->FindOrBuildElement("Mg");
	auto elAl = nistManager->FindOrBuildElement("Al");
	auto elSi = nistManager->FindOrBuildElement("Si");
	auto elP = nistManager->FindOrBuildElement("P");
	// auto elS = nistManager->FindOrBuildElement("S");
	auto elK = nistManager->FindOrBuildElement("K");
	auto elCa = nistManager->FindOrBuildElement("Ca");
	auto elTi = nistManager->FindOrBuildElement("Ti");
	auto elCr = nistManager->FindOrBuildElement("Cr");
	auto elMn = nistManager->FindOrBuildElement("Mn");
	auto elFe = nistManager->FindOrBuildElement("Fe");
	// auto elCu = nistManager->FindOrBuildElement("Cu");
	// auto elZn = nistManager->FindOrBuildElement("Zn");
	auto elSm = nistManager->FindOrBuildElement("Sm");
	auto elEu = nistManager->FindOrBuildElement("Eu");
	auto elGd = nistManager->FindOrBuildElement("Gd");
	auto elTh = nistManager->FindOrBuildElement("Th");
	auto elU = nistManager->FindOrBuildElement("U");

	// element w/ enriched isotope
	// auto el3He = new G4Element("helium3", "3He", 1);
	// el3He->AddIsotope(iso3He, 1);

	auto iso6Li = new G4Isotope("Li6", 3, 6, 6.02 * g / mole);
	auto iso7Li = new G4Isotope("Li7", 3, 7, 7.02 * g / mole);
	auto elIsoLi = new G4Element("isoLi", "Li", 2);
	elIsoLi->AddIsotope(iso6Li, 95.5 * perCent);
	elIsoLi->AddIsotope(iso7Li, 4.5 * perCent);

	// element w/ thermal scattering
	// auto elTSH_H2O = new G4Element("TS_H_of_Water", "h_water", 1, 1.0079*g/mole);
	auto elTSAl = new G4Element("TS_Aluminium_Metal", "al_metal", 13, 26.982 * g / mole);
	// auto elTSFe = new G4Element("TS_Iron_Metal", "fe_metal", 26, 55.845*g/mole);

	/* material */
	Vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
	// Air = nistManager->FindOrBuildMaterial("G4_AIR");

	// Water w/ thermal scattering
	// TS_Water = new G4Material("TS_water", 1.*g/cm3, 2, kStateLiquid, roomTemp);
	// TS_Water->AddElement(elTSH_H2O, 2);
	// TS_Water->AddElement(elO, 1);

	// 3He gas (10atm)
	// Helium3_10atm = new G4Material("Helium3_10atm", 1.346*mg/cm3, 1, kStateGas, roomTemp);
	// Helium3_10atm->AddElement(el3He, 1);

	// Al w/ thermal scattering
	// TS_Al = new G4Material("TS_Al", 2.7*g/cm3, 1, kStateSolid, roomTemp);
	// TS_Al->AddElement(elTSAl, 1);

	// Fe w/ thermal scattering
	// TS_Fe = new G4Material("TS_Fe", 7.874*g/cm3, 1, kStateSolid, roomTemp);
	// TS_Fe->AddElement(elTSFe, 1);

	/*
		//Al A5052 w/ thermal scattering
		std::vector< std::pair<G4Element*, G4double> > vEl_TS_A5052{
			{elMg, 0.025},
			{elTSAl, 0.96775},
			{elSi, 0.00125},
			{elCr, 0.0025},
			{elMn, 0.0005},
			{elTSFe, 0.002},
			{elCu, 0.0005},
			{elZn, 0.0005}
		};
		TS_A5052 = new G4Material("TS_A5052", 2.68*g/cm3, vEl_TS_A5052.size(), kStateSolid, roomTemp);
		for(const auto& el : vEl_TS_A5052) TS_A5052->AddElement(el.first, el.second);
	*/

	/*
		//B4C sheet (Mirrobor)
		std::vector< std::pair<G4Element*, G4double> > vEl_Mirrobor{
			{elH, 0.01},
			{elB, 0.6321},
			{elC, 0.3061},
			{elO, 0.05},
			{elFe, 0.0018}
		};
		Mirrobor = new G4Material("Mirrobor", 1.36*g/cm3, vEl_Mirrobor.size(), kStateSolid, roomTemp);
		for(const auto& el : vEl_Mirrobor) Mirrobor->AddElement(el.first, el.second);
	*/

	// Feldspathic crust surface based on feldspathic lunar meteorites
	std::vector<std::pair<G4Element *, G4double>> vEl_FHT{
		{elH, fLunarSoilHfrac},
		{elO, 0.453142467 * (1. - fLunarSoilHfrac)},
		{elNa, 0.002596501 * (1. - fLunarSoilHfrac)},
		{elMg, 0.032563938 * (1. - fLunarSoilHfrac)},
		{elAl, 0.149248695 * (1. - fLunarSoilHfrac)},
		{elSi, 0.20894341 * (1. - fLunarSoilHfrac)},
		{elP, 0.000117834 * (1. - fLunarSoilHfrac)},
		{elK, 0.00022414 * (1. - fLunarSoilHfrac)},
		{elCa, 0.116494595 * (1. - fLunarSoilHfrac)},
		{elTi, 0.001318554 * (1. - fLunarSoilHfrac)},
		{elCr, 0.000656834 * (1. - fLunarSoilHfrac)},
		{elMn, 0.000487908 * (1. - fLunarSoilHfrac)},
		{elFe, 0.034201413 * (1. - fLunarSoilHfrac)},
		{elSm, 0.0000011 * (1. - fLunarSoilHfrac)},
		{elEu, 0.00000078 * (1. - fLunarSoilHfrac)},
		{elGd, 0.0000013 * (1. - fLunarSoilHfrac)},
		{elTh, 0.00000037 * (1. - fLunarSoilHfrac)},
		{elU, 0.00000016 * (1. - fLunarSoilHfrac)},
	};
	FHT = new G4Material("FHT", 1.6 * g / cm3, vEl_FHT.size(), kStateSolid, moonTemp);
	for (const auto &el : vEl_FHT)
		FHT->AddElement(el.first, el.second);

	// Lunar lander
	Lander = new G4Material("Lander", 2.018 * 1e-1 * g / cm3, 1, kStateSolid, landerTemp);
	Lander->AddElement(elTSAl, 1);

	// Lunar lander fuel
	constexpr G4double tankVolume = M_PI * 40 * 40 * 100 * 4; // volume for 4 tanks (cm3)

	G4double fuelMass = fLanderFuelMass / g;				// total fuel mass (g)
	G4double fuelDensity = fuelMass / tankVolume * g / cm3; // density for fuel & oxidizer (g/cm3)
	if (fuelMass < 0.1)
		fuelDensity = 1e-25 * g / cm3;

	// fuel: hydrazine N2H4, 1.103 g/cm3
	Fuel = new G4Material("Fuel", fuelDensity, 2, kStateLiquid, landerTemp);
	Fuel->AddElement(elH, 4);
	Fuel->AddElement(elN, 2);

	// oxidizer: N2O4, 1.443 g/cm3
	Oxidizer = new G4Material("Oxidizer", fuelDensity, 2, kStateLiquid, landerTemp);
	Oxidizer->AddElement(elN, 2);
	Oxidizer->AddElement(elO, 4);

	// EJ270 plastic scintillator (Eljen Technology, https://eljentechnology.com/products/plastic-scintillators/ej-270)
	constexpr G4double EJ270_H_frac = 5.20 * 1;
	constexpr G4double EJ270_C_frac = 4.62 * 12;
	constexpr G4double EJ270_Li_frac = 0.0575 * 6;
	constexpr G4double EJ270_O_frac = 0.479 * 16;
	constexpr G4double EJ270_N_frac = 0.0931 * 14;
	constexpr G4double EJ270_tot_frac = EJ270_H_frac + EJ270_C_frac + EJ270_Li_frac + EJ270_O_frac + EJ270_N_frac;
	EJ270 = new G4Material("EJ270", 1.14 * g / cm3, 5);
	EJ270->AddElement(elH, EJ270_H_frac / EJ270_tot_frac * 100 * perCent);
	EJ270->AddElement(elC, EJ270_C_frac / EJ270_tot_frac * 100 * perCent);
	EJ270->AddElement(elIsoLi, EJ270_Li_frac / EJ270_tot_frac * 100 * perCent);
	EJ270->AddElement(elO, EJ270_O_frac / EJ270_tot_frac * 100 * perCent);
	EJ270->AddElement(elN, EJ270_N_frac / EJ270_tot_frac * 100 * perCent);
}

void Geometry::PrintParameters()
{
	G4cout << std::boolalpha;
	G4cout << G4endl;
	G4cout << "----- Variable parameters -----" << G4endl;
	G4cout << "Lunar soil H fraction: " << fLunarSoilHfrac << G4endl;
	G4cout << "Lander fuel mass: " << G4BestUnit(fLanderFuelMass, "Mass") << G4endl;
	G4cout << "-------------------------------" << G4endl;
	G4cout << G4endl;
	G4cout << std::noboolalpha;
}
