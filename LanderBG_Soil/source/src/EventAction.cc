#include "EventAction.hh"
#include "RunAction.hh"
#include "SensitiveVolume.hh"
#include "TargetHit.hh"
#include "Analysis.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>
#include <map>
#include <vector>

EventAction::EventAction(RunAction *runAction)
	: G4UserEventAction(),
	  fRunAction(runAction),
	  fTargetHCID(-1)
{
}

EventAction::~EventAction()
{
}

HitsCollection *EventAction::GetHitsCollection(G4int hcID, const G4Event *event) const
{
	auto hitsCollection = static_cast<HitsCollection *>(event->GetHCofThisEvent()->GetHC(hcID));
	if (!hitsCollection)
	{
		G4ExceptionDescription msg;
		msg << "Cannot access hitsCollection ID " << hcID;
		G4Exception("EventAction::GetHitsCollection()", "MyCode0003", FatalException, msg);
	}

	return hitsCollection;
}

void EventAction::PrintEventStatistics() const
{
}

void EventAction::BeginOfEventAction(const G4Event *)
{
	fTrackIDToParticleNameMap.clear();
	fRunAction->ClearScorerData();
}

void EventAction::EndOfEventAction(const G4Event *event)
{
	// Get hits collections IDs (only once)
	if (fTargetHCID == -1)
	{
		fTargetHCID = G4SDManager::GetSDMpointer()->GetCollectionID("HitsCollection");
	}
	auto targetHC = GetHitsCollection(fTargetHCID, event);
	auto analysisManager = G4AnalysisManager::Instance();

	// Primary particle
	auto primaryVertex = event->GetPrimaryVertex();
	auto primaryParticle = primaryVertex->GetPrimary();
	auto primEnergy = primaryParticle->GetKineticEnergy() / MeV;

	analysisManager->FillH1(0, primEnergy);
	analysisManager->FillH1(1, primaryVertex->GetX0() / m);
	analysisManager->FillH1(2, primaryVertex->GetY0() / m);
	analysisManager->FillH1(3, primaryVertex->GetZ0() / m);
	analysisManager->FillH1(4, std::acos(primaryParticle->GetMomentumDirection().z()) / M_PI * 180);
	analysisManager->FillH1(5, std::atan2(primaryParticle->GetMomentumDirection().y(), primaryParticle->GetMomentumDirection().x()) / M_PI * 180);
	analysisManager->FillH2(0, primaryVertex->GetX0() / m, primaryVertex->GetY0() / m);
	analysisManager->FillH2(1, primaryVertex->GetX0() / m, primaryVertex->GetZ0() / m);

	// Hit container
	std::map<G4int, G4double> neutronEntranceEnergyMap;
	G4bool EnteringFlag = false;
	G4bool scorerHitFlag = false;
	for (size_t i = 0, nHit = targetHC->entries(); i < nHit; ++i)
	{
		auto hit = (*targetHC)[i];

		if (hit->isEnteringVolume() && hit->GetParticleName() == "neutron")
		{
			EnteringFlag = true;
			fRunAction->AddTrackID(hit->GetTrackID());
			fRunAction->AddNeutronEntranceEnergy(hit->GetPreKineticEnergy() / MeV);
			fRunAction->AddNeutronEntrancePos(hit->GetPreStepLocalPosition().x() / mm, hit->GetPreStepLocalPosition().y() / mm, hit->GetPreStepLocalPosition().z() / mm);

			neutronEntranceEnergyMap[hit->GetTrackID()] = hit->GetPreKineticEnergy() / MeV;
		}

		G4double edep = hit->GetEdep() / MeV;
		if (edep <= 0.0)
			continue;

		scorerHitFlag = true;
		G4String parentParticleName = GetTrackParticleName(hit->GetParentID());

		auto neIt = neutronEntranceEnergyMap.find(hit->GetParentID());
		G4double parentEntranceEnergy = (neIt != neutronEntranceEnergyMap.end()) ? neIt->second : -1.0; // -1 = 親が直接入射した中性子ではない(孫世代など)

		fRunAction->AddParentEntranceEnergy(parentEntranceEnergy);
		fRunAction->AddParentID(hit->GetParentID());
		fRunAction->AddParentParticleID(parentParticleName == "Primary" ? -1 : ParticleID(parentParticleName));
		fRunAction->AddParticleID(ParticleID(hit->GetParticleName()));
		fRunAction->AddStepProcessID(ProcessID(hit->GetStepProcessName()));
		fRunAction->AddCreatorProcessID(ProcessID(hit->GetCreatorProcessName()));
		fRunAction->AddPreEnergy(hit->GetPreKineticEnergy() / MeV);
		fRunAction->AddPostEnergy(hit->GetPostKineticEnergy() / MeV);
		fRunAction->AddEdep(edep);
		fRunAction->AddGTime(hit->GetPostGlobalTime() / ns);
		fRunAction->AddPrePos(hit->GetPreStepLocalPosition().x() / mm, hit->GetPreStepLocalPosition().y() / mm, hit->GetPreStepLocalPosition().z() / mm);
		fRunAction->AddPostPos(hit->GetPostStepLocalPosition().x() / mm, hit->GetPostStepLocalPosition().y() / mm, hit->GetPostStepLocalPosition().z() / mm);
		fRunAction->AddPreMom(hit->GetPreMomentumDirection().x(), hit->GetPreMomentumDirection().y(), hit->GetPreMomentumDirection().z());
	}

	if (EnteringFlag)
	{
		analysisManager->FillNtupleIColumn(0, 0, event->GetEventID());
		analysisManager->AddNtupleRow(0);
		if (scorerHitFlag)
		{
			analysisManager->FillNtupleIColumn(1, 0, event->GetEventID());
			analysisManager->AddNtupleRow(1);
		}
	}
}

G4int EventAction::ParticleID(G4String particle)
{
	G4int id = -1;

	if (particle.find("nu_") != std::string::npos)
		id = -2;
	else if (particle == "proton")
		id = 0;
	else if (particle == "neutron")
		id = 1;
	else if (particle == "gamma")
		id = 2;
	else if (particle == "deuteron")
		id = 3;
	else if (particle == "triton")
		id = 4;
	else if (particle == "He3")
		id = 5;
	else if (particle == "alpha")
		id = 6;
	else if (particle == "e-")
		id = 7;
	else if (particle == "e+")
		id = 8;
	else if (particle == "mu-")
		id = 9;
	else if (particle == "mu+")
		id = 10;
	else if (particle == "pi-")
		id = 11;
	else if (particle == "pi+")
		id = 12;
	else if (particle == "kaon-")
		id = 13;
	else if (particle == "kaon+")
		id = 14;
	else if (particle == "kaon0S")
		id = 15;
	else if (particle == "kaon0L")
		id = 16;
	else if (particle == "lambda")
		id = 17;
	else if (particle == "sigma-")
		id = 18;
	else if (particle == "sigma+")
		id = 19;
	else
	{
		id = 20;
		// G4cout << "other particle: " << particle << G4endl;
	}

	return id;
}

G4int EventAction::ProcessID(G4String process)
{
	G4int id = -1;

	if (process == "Primary")
		id = 0;
	else if (process == "protonInelastic")
		id = 1;
	else if (process == "dInelastic")
		id = 2;
	else if (process == "tInelastic")
		id = 3;
	else if (process == "He3Inelastic")
		id = 4;
	else if (process == "alphaInelastic")
		id = 5;
	else if (process == "ionInelastic")
		id = 6;
	else if (process == "pi+Inelastic")
		id = 7;
	else if (process == "pi-Inelastic")
		id = 8;
	else if (process == "kaon+Inelastic")
		id = 9;
	else if (process == "kaon-Inelastic")
		id = 10;
	else if (process == "kaon0SInelastic")
		id = 11;
	else if (process == "kaon0LInelastic")
		id = 12;
	else if (process == "lambdaInelastic")
		id = 13;
	else if (process == "sigma+Inelastic")
		id = 14;
	else if (process == "sigma-Inelastic")
		id = 15;
	else if (process == "neutronInelastic")
		id = 16;
	else if (process == "nCapture")
		id = 17;
	else if (process == "nFission")
		id = 18;
	else if (process == "Decay")
		id = 19;
	else if (process == "RadioactiveDecay")
		id = 20;
	else if (process == "hBrems")
		id = 21;
	else if (process == "eBrem")
		id = 22;
	else if (process == "phot")
		id = 23;
	else if (process == "compt")
		id = 24;
	else if (process == "annihil")
		id = 25;
	else if (process == "photonNuclear")
		id = 26;
	else if (process == "electronNuclear")
		id = 27;
	else if (process == "positronNuclear")
		id = 28;
	else if (process == "hBertiniCaptureAtRest")
		id = 29;
	else if (process == "muMinusCaptureAtRest")
		id = 30;
	else if (process == "hadElastic")
		id = 31;
	else if (process == "hIoni")
		id = 32;
	else if (process == "ionIoni")
		id = 33;
	else if (process == "eIoni")
		id = 34;
	else if (process == "muIoni")
		id = 35;
	else if (process == "conv")
		id = 36;
	else
	{
		id = 37;
		// G4cout << "other process: " << process << G4endl;
	}

	std::cout << "Process: " << process << ", ID: " << id << std::endl;

	return id;
}
