#ifndef B2SteppingAction_h
#define B2SteppingAction_h 1
#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>

class G4LogicalVolume;

class B2SteppingAction : public G4UserSteppingAction {
public:
  B2SteppingAction() = default;
  void UserSteppingAction(const G4Step*) override;
private:
  void Init();
  std::vector<G4double> fPlaneZ;                // 各仮想Cd面のワールドZ
  G4double fZFront = 0.;                        // 前面(深さ0)のワールドZ
  const G4LogicalVolume* fCubeLV = nullptr;     // EJ-270 cube
};
#endif
