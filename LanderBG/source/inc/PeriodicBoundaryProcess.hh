#ifndef PeriodicBoundaryProcess_h
#define PeriodicBoundaryProcess_h 1

#include "ParticleChangeForPeriodic.hh"

#include "G4DynamicParticle.hh"
#include "G4OpticalPhoton.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoTau.hh"
#include "G4TransportationManager.hh"
#include "G4VDiscreteProcess.hh"
#include "globals.hh"

class G4Step;

enum ProcessStatus
{
	Undefined,
	Reflection,
	Cycling,
	StepTooSmall,
	NotAtBoundary
};

class PeriodicBoundaryProcess : public G4VDiscreteProcess
{
public:
	explicit PeriodicBoundaryProcess(const G4String& processName = "PBC", G4ProcessType type = fNotDefined, bool per_x = true, bool per_y = true, bool per_z = false);
	~PeriodicBoundaryProcess() override = default;

	PeriodicBoundaryProcess(const PeriodicBoundaryProcess& right) = delete;
	PeriodicBoundaryProcess& operator=(const PeriodicBoundaryProcess& right) = delete;

	G4bool IsApplicable(const G4ParticleDefinition&) override;
	G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) override;
	ProcessStatus GetStatus() const;
	G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

protected:
	ParticleChangeForPeriodic fParticleChange;

private:
	void BoundaryProcessVerbose();
	std::map<ProcessStatus, G4String> fStatusMessages = {
		{Undefined, " *** Undefined *** "},
		{NotAtBoundary, " *** NotAtBoundary *** "},
		{Reflection, " *** Reflection *** "},
		{Cycling, " *** Periodic *** "},
		{StepTooSmall, " *** StepTooSmall *** "}
	};

	ProcessStatus fTheStatus = Undefined;
	G4ThreeVector fOldPosition;
	G4ThreeVector fNewPosition;
	G4ThreeVector fOldMomentum;
	G4ThreeVector fNewMomentum;
	G4ThreeVector fOldPolarization;
	G4ThreeVector fNewPolarization;
	G4ThreeVector fTheGlobalNormal;
	G4double fkCarTolerance;
	G4bool fPeriodicX = true;
	G4bool fPeriodicY = true;
	G4bool fPeriodicZ = false;
};

inline G4bool PeriodicBoundaryProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
	G4bool applicable = true;

	if(&aParticleType==G4AntiNeutrinoE::AntiNeutrinoE()) applicable = false;
	else if(&aParticleType==G4NeutrinoE::NeutrinoE()) applicable = false;
	else if(&aParticleType==G4AntiNeutrinoMu::AntiNeutrinoMu()) applicable = false;
	else if(&aParticleType==G4NeutrinoMu::NeutrinoMu()) applicable = false;
	else if(&aParticleType==G4AntiNeutrinoTau::AntiNeutrinoTau()) applicable = false;
	else if(&aParticleType==G4NeutrinoTau::NeutrinoTau()) applicable = false;
	else if(&aParticleType==G4OpticalPhoton::OpticalPhoton()) applicable = false;

	return applicable;
}

inline ProcessStatus PeriodicBoundaryProcess::GetStatus() const
{
	return fTheStatus;
}

#endif

