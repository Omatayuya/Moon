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
/// \file B2aDetectorConstruction.hh
/// \brief Definition of the B2aDetectorConstruction class

#ifndef B2aDetectorConstruction_h
#define B2aDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"

#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class B2aDetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class B2aDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    B2aDetectorConstruction();
    virtual ~B2aDetectorConstruction();

    void AddSensitiveDetector(G4String name) { fSensitiveLogVolName.emplace_back(name); };

public:
    virtual G4VPhysicalVolume *Construct();
    virtual void ConstructSDandField();

    // Set methods
    void SetTargetMaterial(G4String);
    void SetTargetMaterial_0(G4String);
    void SetTargetMaterial_1(G4String);
    void SetFilmMaterial(G4String);
    void SetChamberMaterial(G4String);
    void SetMaxStep(G4double);
    void SetCheckOverlaps(G4bool);
    void SetAbsorberMaterial(G4String); // byHN

    const G4ThreeVector &GetSensorCenter() const { return fSensorCenter; }
    const G4ThreeVector &GetSensorHalfSize() const { return fSensorHalfSize; }
    G4double GetGdsheetT() const { return fGdsheetT; }

private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume *DefineVolumes();
    void ConstructSensor();
    void ConstructLab();

    // data members
    G4int fNbOfChambers;

    G4VPhysicalVolume *worldPV;
    std::vector<G4String> fSensitiveLogVolName;

    G4LogicalVolume *fLogicTarget; // pointer to the logical Target
    G4LogicalVolume *fLogicTarget_0;
    G4LogicalVolume *fLogicTarget_1;
    G4LogicalVolume *fLogicTarget_2;
    G4LogicalVolume *fLogicTarget_3;
    G4LogicalVolume *fLogicTarget_4;
    G4LogicalVolume *fLogicTarget_5;
    G4LogicalVolume *fLogicTarget_6;
    G4LogicalVolume *fLogicTarget_7;
    // pointer to the logical Target
    G4LogicalVolume *fLogicFilm;     // pointer to the logical Film
    G4LogicalVolume **fLogicChamber; // pointer to the logical Chamber
    G4LogicalVolume *fLogicSphere;   // pointer to the logical Sphere
    G4LogicalVolume *fLogicAbsorber; // pointer to the logical Absorber(byHN)
    G4LogicalVolume *fLogicAbsorber_4;
    G4LogicalVolume *fLogicAbsorber_5;
    G4LogicalVolume *fLogicAbsorber_6;
    G4LogicalVolume *fLogicAbsorber_7;
    G4Material *fTargetMaterial;      // pointer to the target  material
    G4Material *fTargetMaterial_GAGG; // pointer to the target  material (GAGG)
    G4Material *fFilmMaterial;        // pointer to the film  material
    G4Material *fChamberMaterial;     // pointer to the chamber material
    G4Material *fSphereMaterial;      // pointer to the hemisphere material
    G4Material *fAbsorberMaterial;    // pointer to the absorber material(byHN)

    G4UserLimits *fStepLimit; // pointer to user step limits

    B2aDetectorMessenger *fMessenger; // messenger

    static G4ThreadLocal G4GlobalMagFieldMessenger *fMagFieldMessenger;
    // magnetic field messenger

    G4bool fCheckOverlaps; // option to activate checking of volumes overlaps

    G4Material *PbBlock;
    G4Material *G10;
    G4Material *PTFE;
    G4Material *Concrete;
    G4Material *Cu;
    G4Material *Air;

    G4Material *EJ270;
    G4Material *Mirrobor;
    G4Material *Polyethylene;
    G4Material *Gdsheet;

    G4ThreeVector fSensorCenter;
    G4ThreeVector fSensorHalfSize;
    G4double fGdsheetT = 0.;
};

struct CadParam
{
    G4String filename;
    G4Material *material;
    G4VisAttributes vis = G4VisAttributes::GetInvisible();
    G4bool sensVolFlag = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
