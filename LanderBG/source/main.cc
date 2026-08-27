#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
// #include "G4GDMLParser.hh"
#include "Randomize.hh"
#include "time.h"

#include "G4PhysListFactory.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicParameters.hh"
// #include "G4EmParameters.hh"

#include "SensitiveVolume.hh"
#include "Geometry.hh"
#include "UserActionInitialization.hh"
#include "PeriodicBoundaryPhysics.hh"

int main(int argc, char **argv)
{
	// ============= [ Setting up the application environment ] ================
	G4String mainMacro = "GlobalSetup.mac";
	G4bool visFlag = true;

	G4bool messengerFlag = true;
	G4String messengerMacro = "preInitSetup.mac";
	if (!visFlag)
		messengerMacro.insert(0, "../../");
	// =========================================================================

	G4Random::setTheSeed(static_cast<G4long>(time(NULL)));

	// default run manager
	auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
#ifdef G4MULTITHREADED
	if (visFlag)
		runManager->SetNumberOfThreads(1);
	else
		runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
#endif

	// mandatory user initialization: Geometry
	runManager->SetUserInitialization(new Geometry);

	// mandatory user initialization: Physics-List
	auto physListFactory = new G4PhysListFactory();
	auto physList = physListFactory->GetReferencePhysList("Shielding_HPT");

	auto hpManager = G4ParticleHPManager::GetInstance();
	hpManager->SetDoNotAdjustFinalState(true);
	hpManager->SetSkipMissingIsotopes(true);
	hpManager->SetUseOnlyPhotoEvaporation(false);
	hpManager->SetNeglectDoppler(false);
	hpManager->SetProduceFissionFragments(false); // default fission model
	hpManager->SetUseWendtFissionModel(false);	  // alternative fission model
	hpManager->SetUseNRESP71Model(false);

	auto hadronicParam = G4HadronicParameters::Instance();
	hadronicParam->SetTimeThresholdForRadioactiveDecay(1e+12 * CLHEP::year);
	hadronicParam->SetMinEnergyTransitionFTF_Cascade(1.5 * CLHEP::GeV); // 追加：デフォルト3.0GeVから引き下げ
	hadronicParam->SetMaxEnergyTransitionFTF_Cascade(3.0 * CLHEP::GeV); // 追加：デフォルト6.0GeVから引き下げ

	// auto emParam = G4EmParameters::Instance();
	// emParam->SetFluo(true);
	// emParam->SetAuger(true);
	// emParam->SetPixe(true);
	// emParam->SetBeardenFluoDir(true);
	// emParam->SetANSTOFluoDir(true);

	auto pbc = new PeriodicBoundaryPhysics("Periodic", true, true, false);
	pbc->SetVerboseLevel(0);
	physList->RegisterPhysics(pbc);

	runManager->SetUserInitialization(physList);

	// user initialization: User Actions
	runManager->SetUserInitialization(new UserActionInitialization());

	// UI manager
	auto uiManager = G4UImanager::GetUIpointer();
	if (messengerFlag)
		uiManager->ApplyCommand("/control/execute " + messengerMacro);

	// Initialize G4 kernel
	runManager->Initialize();

	// visualization environment
	G4VisExecutive *visManager = nullptr;
	if (visFlag)
	{
		visManager = new G4VisExecutive;
		visManager->Initialize();
	}

	// interactive session or process macro
	G4UIExecutive *ui = nullptr;
	if (visFlag)
		ui = new G4UIExecutive(argc, argv, "qt");

	uiManager->ApplyCommand("/control/execute " + mainMacro);

	if (visFlag)
		ui->SessionStart();

	delete ui;
	delete visManager;
	delete runManager;

	return 0;
}
