#ifndef ParticleChangeForPeriodic_hh
#define ParticleChangeForPeriodic_hh 1

#include "G4VParticleChange.hh"
#include "globals.hh"

class G4DynamicParticle;

class ParticleChangeForPeriodic : public G4VParticleChange
{
public:
	ParticleChangeForPeriodic();
	~ParticleChangeForPeriodic() override = default;

	G4Step* UpdateStepForPostStep(G4Step* Step) override;

	void InitializeForPostStep(const G4Track&);

	void AddSecondary(G4DynamicParticle* aParticle);

	const G4ThreeVector& GetProposedMomentumDirection() const;
	void ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz);
	void ProposeMomentumDirection(const G4ThreeVector& Pfinal);

	const G4ThreeVector& GetProposedPolarization() const;
	void ProposePolarization(const G4ThreeVector& dir);
	void ProposePolarization(G4double Px, G4double Py, G4double Pz);

	const G4ThreeVector& GetProposedPosition() const;
	void ProposePosition(const G4ThreeVector& pos);
	void ProposePosition(G4double x, G4double y, G4double z);

	const G4Track* GetCurrentTrack() const;

	void DumpInfo() const override;

	ParticleChangeForPeriodic(const ParticleChangeForPeriodic& right) = delete;
	ParticleChangeForPeriodic& operator=(const ParticleChangeForPeriodic& right) = delete;

private:
	const G4Track* fTrack{};
	G4ThreeVector fProposedMomentumDirection;
	G4ThreeVector fProposedPolarization;
	G4ThreeVector fProposedPosition;

};

inline const G4ThreeVector& ParticleChangeForPeriodic::GetProposedMomentumDirection() const
{
	return fProposedMomentumDirection;
}

inline void ParticleChangeForPeriodic::ProposeMomentumDirection(const G4ThreeVector& dir)
{
	fProposedMomentumDirection = dir;
}

inline void ParticleChangeForPeriodic::ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz)
{
	fProposedMomentumDirection.setX(Px);
	fProposedMomentumDirection.setY(Py);
	fProposedMomentumDirection.setZ(Pz);
}

inline const G4ThreeVector& ParticleChangeForPeriodic::GetProposedPolarization() const
{
	return fProposedPolarization;
}

inline void ParticleChangeForPeriodic::ProposePolarization(const G4ThreeVector& dir)
{
	fProposedPolarization = dir;
}

inline void ParticleChangeForPeriodic::ProposePolarization(G4double Px, G4double Py, G4double Pz)
{
	fProposedPolarization.setX(Px);
	fProposedPolarization.setY(Py);
	fProposedPolarization.setZ(Pz);
}

inline const G4ThreeVector& ParticleChangeForPeriodic::GetProposedPosition() const
{
	return fProposedPosition;
}

inline void ParticleChangeForPeriodic::ProposePosition(const G4ThreeVector& dir)
{
	fProposedPosition = dir;
}

inline void ParticleChangeForPeriodic::ProposePosition(G4double Px, G4double Py, G4double Pz)
{
	fProposedPosition.setX(Px);
	fProposedPosition.setY(Py);
	fProposedPosition.setZ(Pz);
}

inline void ParticleChangeForPeriodic::InitializeForPostStep(const G4Track& track)
{
	theStatusChange = track.GetTrackStatus();
	theLocalEnergyDeposit = 0.0;
	theNonIonizingEnergyDeposit = 0.0;
	InitializeSecondaries();
	theParentWeight = track.GetWeight();
	isParentWeightProposed = false;
	fProposedMomentumDirection = track.GetMomentumDirection();
	fProposedPolarization = track.GetPolarization();
	fProposedPosition = track.GetPosition();
	fTrack = &track;
}

inline const G4Track* ParticleChangeForPeriodic::GetCurrentTrack() const
{
	return fTrack;
}
#endif

