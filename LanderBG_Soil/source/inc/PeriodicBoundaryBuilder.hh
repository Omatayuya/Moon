#ifndef PeriodicBoundaryBuilder_h
#define PeriodicBoundaryBuilder_h 1

#include "LogicalVolumePeriodic.hh"

class PeriodicBoundaryBuilder
{
public:
	PeriodicBoundaryBuilder() = default;
	~PeriodicBoundaryBuilder() = default;

	G4LogicalVolume* Construct(G4LogicalVolume*);

private:
	LogicalVolumePeriodic* fLogicalPeriodic = nullptr;

};
#endif

