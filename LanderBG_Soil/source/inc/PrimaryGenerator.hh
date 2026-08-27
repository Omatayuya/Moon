#ifndef PrimaryGenerator_h
#define PrimaryGenerator_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
class G4GeneralParticleSource;

class PrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
	PrimaryGenerator();
	~PrimaryGenerator();

	void GeneratePrimaries(G4Event*);

private:
	G4GeneralParticleSource* fParticleGPS;
};

#endif

