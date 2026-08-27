#include "ParticleChangeForPeriodic.hh"

#include "G4DynamicParticle.hh"
#include "G4Step.hh"
#include "G4Track.hh"

ParticleChangeForPeriodic::ParticleChangeForPeriodic() : G4VParticleChange()
{}

G4Step* ParticleChangeForPeriodic::UpdateStepForPostStep(G4Step* pStep)
{
	auto pPostStepPoint = pStep->GetPostStepPoint();
	pPostStepPoint->SetMomentumDirection(fProposedMomentumDirection);
	pPostStepPoint->SetPolarization(fProposedPolarization);
	pPostStepPoint->SetPosition(fProposedPosition);

	if(isParentWeightProposed) pPostStepPoint->SetWeight(theParentWeight);

	pStep->AddTotalEnergyDeposit(theLocalEnergyDeposit);
	pStep->AddNonIonizingEnergyDeposit(theNonIonizingEnergyDeposit);

	return pStep;
}

void ParticleChangeForPeriodic::AddSecondary(G4DynamicParticle* aParticle)
{
	auto aTrack = new G4Track(aParticle, fTrack->GetGlobalTime(), fTrack->GetPosition());
	aTrack->SetTouchableHandle(fTrack->GetTouchableHandle());
	G4VParticleChange::AddSecondary(aTrack);
}

void ParticleChangeForPeriodic::DumpInfo() const
{
	G4VParticleChange::DumpInfo();
	G4int oldprc = G4cout.precision(3);

	G4cout << "        Momentum Direction: " << std::setw(20) << fProposedMomentumDirection << G4endl;
	G4cout << "        Polarization: " << std::setw(20) << fProposedPolarization << G4endl;
	G4cout << "        Position: " << std::setw(20) << fProposedPosition << G4endl;
	G4cout.precision(oldprc);
}

