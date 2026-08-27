//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2TrackerHit.hh
/// \brief Definition of the B2TrackerHit class

#ifndef B2TrackerHit_h
#define B2TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"


/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fGTime, fLTime, fPTime, fDTime, fPostPos

class B2TrackerHit : public G4VHit
{
  public:
    B2TrackerHit();
    B2TrackerHit(const B2TrackerHit&);
    virtual ~B2TrackerHit();

    // operators
    const B2TrackerHit& operator=(const B2TrackerHit&);
    G4bool operator==(const B2TrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetEventID  (G4int event)             { fEventID = event; };//yk
	void SetTrackID  (G4int track)             { fTrackID = track; };
    void SetChamberNb(G4int chamb)             { fChamberNb = chamb; };
    void SetEdep     (G4double de)             { fEdep = de; };
    void SetStepLength(G4double stepLength)    { fStepLength = stepLength; };
    void SetGlobalTime (G4double gtime)        { fGTime = gtime; };//yk
    void SetLocalTime  (G4double ltime)        { fLTime = ltime; };//yk
    void SetProperTime (G4double ptime)        { fPTime = ptime; };//yk
    void SetDeltaTime  (G4double dtime)        { fDTime = dtime; };//yk
	void SetPrePos   (G4ThreeVector prexyz)    { fPrePos = prexyz; };//yk
    void SetPostPos  (G4ThreeVector postxyz)   { fPostPos = postxyz; };
	void SetParticleName(G4String pname)       { fPname = pname; };//yk
	void SetParticleID(G4int pid)              { fPID = pid; };//yk
	void SetParentParticleID(G4int ppid)       { fPPID = ppid; };//yk
	void SetPreProcessName(G4String preproc)   { fPreProc = preproc; };//yk
	void SetPostProcessName(G4String postproc) { fPostProc = postproc; };//yk
	void SetCreatorProcessName(G4String cproc) { fCProc = cproc; };//yk
	void SetCreatorModelName(G4String cmodel)  { fCModel = cmodel; };//yk test
	void SetPreKineticEnergy(G4double preke)   { fPreKinE = preke; };//yk
	void SetPostKineticEnergy(G4double postke) { fPostKinE = postke; };//yk
	void SetLogicalVolumeAtVertexName(G4String llvn) { fLLVname = llvn; };//yk
    void SetIsEntering(G4bool b){fIsEntering = b;}

    // Get methods
	G4int GetEventID() const               { return fEventID; };//yk
    G4int GetTrackID() const               { return fTrackID; };
    G4int GetChamberNb() const             { return fChamberNb; };
    G4double GetEdep() const               { return fEdep; };
    G4double GetStepLength() const         { return fStepLength; };
    G4double GetGlobalTime() const         { return fGTime; };//yk
    G4double GetLocalTime() const          { return fLTime; };//yk
    G4double GetProperTime() const         { return fPTime; };//yk
    G4double GetDeltaTime() const          { return fDTime; };//yk
	G4ThreeVector GetPrePos() const        { return fPrePos; };//yk
    G4ThreeVector GetPostPos() const       { return fPostPos; };
	G4String GetParticleName() const       { return fPname; };//yk
	G4int GetParticleID() const            { return fPID; };//yk
	G4int GetParentParticleID() const      { return fPPID; };//yk
	G4String GetPreProcessName() const     { return fPreProc; };//yk
	G4String GetPostProcessName() const    { return fPostProc; };//yk
	G4String GetCreatorProcessName() const { return fCProc; };//yk
	G4String GetCreatorModelName() const   { return fCModel; };//yk test
	G4double GetPreKineticEnergy() const   { return fPreKinE; };//yk
	G4double GetPostKineticEnergy() const  { return fPostKinE; };//yk
	G4String GetLogicalVolumeAtVertex() const { return fLLVname; };//yk
    G4bool GetIsEntering() const { return fIsEntering; }

  private:

      G4int         fEventID;//yk
	  G4int         fTrackID;
      G4int         fChamberNb;
      G4double      fEdep;
      G4double      fStepLength;
      G4double      fGTime;//yk
      G4double      fLTime;//yk
      G4double      fPTime;//yk
      G4double      fDTime;//yk
	  G4ThreeVector fPrePos;//yk
      G4ThreeVector fPostPos;
	  G4String      fPname;//yk
	  G4int         fPID;//yk
	  G4int         fPPID;//yk
	  G4String      fPreProc;//yk
	  G4String      fPostProc;//yk
	  G4String      fCProc;//yk
	  G4String      fCModel;//yk
	  G4double      fPreKinE;//yk
	  G4double      fPostKinE;//yk
	  G4String      fLLVname;//yk
      G4bool        fIsEntering;
    G4int fRunID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<B2TrackerHit> B2TrackerHitsCollection;

extern G4ThreadLocal G4Allocator<B2TrackerHit>* B2TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* B2TrackerHit::operator new(size_t)
{
  if(!B2TrackerHitAllocator)
      B2TrackerHitAllocator = new G4Allocator<B2TrackerHit>;
  return (void *) B2TrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void B2TrackerHit::operator delete(void *hit)
{
  B2TrackerHitAllocator->FreeSingle((B2TrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
