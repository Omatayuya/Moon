#ifndef PeriodicBoundaryPhysics_hh
#define PeriodicBoundaryPhysics_hh 1

#include "G4VPhysicsConstructor.hh"

class PeriodicBoundaryProcess;

class PeriodicBoundaryPhysics : public G4VPhysicsConstructor
{
public:
	explicit PeriodicBoundaryPhysics(const G4String& name = "Periodic", G4bool per_x = true, G4bool per_y = true, G4bool per_z = false);
	~PeriodicBoundaryPhysics() override = default;

	void ConstructParticle() override;
	void ConstructProcess() override;

private:
	G4bool fPeriodicX = true, fPeriodicY = true, fPeriodicZ = false;

	static void ThrowException(const G4String& particleName);
	static void AddDiscreteProcess(PeriodicBoundaryProcess* periodicBoundaryProcess, G4ParticleDefinition& particle, G4ProcessManager* processManager);

};
#endif

