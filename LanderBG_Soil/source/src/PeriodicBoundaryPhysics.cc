#include "PeriodicBoundaryPhysics.hh"
#include "PeriodicBoundaryProcess.hh"

#include "G4PhysicsConstructorFactory.hh"
#include "G4ProcessManager.hh"
#include "globals.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(PeriodicBoundaryPhysics);

PeriodicBoundaryPhysics::PeriodicBoundaryPhysics(const G4String& name, G4bool per_x, G4bool per_y, G4bool per_z)
  : G4VPhysicsConstructor(name), fPeriodicX(per_x), fPeriodicY(per_y), fPeriodicZ(per_z)
{
	verboseLevel = 0;
}

void PeriodicBoundaryPhysics::ConstructParticle()
{}

void PeriodicBoundaryPhysics::ConstructProcess()
{
	if(verboseLevel>0) G4cout << "Constructing cyclic boundary physics process" << G4endl;

	auto* pbc = new PeriodicBoundaryProcess("Cyclic", fNotDefined, fPeriodicX, fPeriodicY, fPeriodicZ);
	if(verboseLevel>0) pbc->SetVerboseLevel(verboseLevel);

	auto aParticleIterator = GetParticleIterator();
	aParticleIterator->reset();

	G4ProcessManager* processManager = nullptr;

	while((*aParticleIterator)()){
		G4ParticleDefinition* particle = aParticleIterator->value();
		G4String particleName = particle->GetParticleName();

		processManager = particle->GetProcessManager();
		if(!processManager){ ThrowException(particleName); return; }

		AddDiscreteProcess(pbc, *particle, processManager);
	}
}

void PeriodicBoundaryPhysics::ThrowException(const G4String& particleName)
{
	std::ostringstream o;
	o << "Particle " << particleName << "without a Process Manager";
	G4Exception("G4PeriodicBoundaryPhysics::ConstructProcess()", "", FatalException, o.str().c_str());
}

void PeriodicBoundaryPhysics::AddDiscreteProcess(PeriodicBoundaryProcess* periodicBoundaryProcess, G4ParticleDefinition& particle, G4ProcessManager* processManager)
{
	if(periodicBoundaryProcess->IsApplicable(particle)) processManager->AddDiscreteProcess(periodicBoundaryProcess);
}

