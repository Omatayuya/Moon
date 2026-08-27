#include "B2SteppingAction.hh"
#include "B2aDetectorConstruction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Neutron.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

void B2SteppingAction::Init()
{
  const auto* det = dynamic_cast<const B2aDetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4ThreeVector c(0, 0, 500.*mm), h(50.*mm, 50.*mm, 100.*mm);
  if (det) { c = det->GetSensorCenter(); h = det->GetSensorHalfSize(); }

  fZFront = c.z() - h.z();                        // 深さ0のワールドZ
  const G4double pitch    = 5.*mm;
  const G4double depthMax = 2.*h.z();             // 200 mm（背面）
  // 深さ 5,10,...,195 mm。前面(0)・背面(200)は含めない
  for (G4double d = pitch; d < depthMax - 0.5*pitch; d += pitch)
    fPlaneZ.push_back(fZFront + d);

  fCubeLV = G4LogicalVolumeStore::GetInstance()->GetVolume("LV_EJ270cube");
  if (!fCubeLV)
    G4Exception("B2SteppingAction::Init", "NoCubeLV", FatalException,
                "LV_EJ270cube not found");
}

void B2SteppingAction::UserSteppingAction(const G4Step* step)
{
  const G4Track* trk = step->GetTrack();
  if (trk->GetParticleDefinition() != G4Neutron::Definition()) return;

  if (fPlaneZ.empty()) Init();

  const G4VPhysicalVolume* pv = step->GetPreStepPoint()->GetPhysicalVolume();
  if (!pv || pv->GetLogicalVolume() != fCubeLV) return;

  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  const G4double zPre = pre->GetPosition().z();
  const G4double dz   = post->GetPosition().z() - zPre;
  if (dz == 0.) return;

  const G4double ekin = pre->GetKineticEnergy();
  const G4double uz   = pre->GetMomentumDirection().z();
  const G4double tPre = pre->GetGlobalTime();
  const G4double tPost= post->GetGlobalTime();

  auto* ana = G4AnalysisManager::Instance();
  const G4int evtID =
      G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  for (std::size_t i = 0; i < fPlaneZ.size(); ++i) {
    const G4double f = (fPlaneZ[i] - zPre) / dz;
    if (f < 0. || f >= 1.) continue;

    const G4double tC = tPre + f * (tPost - tPre);
    const G4int depthMm = static_cast<G4int>(std::lround(fPlaneZ[i] - fZFront));

    G4int col = 0;
    ana->FillNtupleIColumn(1, col++, evtID);
    ana->FillNtupleIColumn(1, col++, trk->GetTrackID());
    ana->FillNtupleIColumn(1, col++, trk->GetParentID());
    ana->FillNtupleIColumn(1, col++, depthMm);              // 5,10,...,195
    ana->FillNtupleDColumn(1, col++, ekin / MeV);
    ana->FillNtupleDColumn(1, col++, uz);
    ana->FillNtupleDColumn(1, col++, tC / ns);
    ana->AddNtupleRow(1);
  }
}
