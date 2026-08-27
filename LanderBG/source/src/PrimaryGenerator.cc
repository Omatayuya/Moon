#include "PrimaryGenerator.hh"
#include "G4GeneralParticleSource.hh"

PrimaryGenerator::PrimaryGenerator()
: fParticleGPS(0)
{
	fParticleGPS = new G4GeneralParticleSource();
}

PrimaryGenerator::~PrimaryGenerator()
{
	delete fParticleGPS;
}

void PrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
	fParticleGPS->GeneratePrimaryVertex(anEvent);
}

