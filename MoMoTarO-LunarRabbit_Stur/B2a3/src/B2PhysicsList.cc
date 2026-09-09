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
/// \file B2PhysicsList.cc
/// \brief Implementation of the B2PhysicsList class

#include "B2PhysicsList.hh"

#include "G4SystemOfUnits.hh"

#include "G4EmLivermorePhysics.hh"//yk
#include "G4EmPenelopePhysics.hh"//yk

// Decays
#include <G4DecayPhysics.hh>
#include <G4RadioactiveDecayPhysics.hh>
// Hadron Elastic scattering
#include <G4HadronElasticPhysics.hh>
#include <G4HadronElasticPhysicsHP.hh>
// Hadron Physics
#include <G4HadronPhysicsFTFP_BERT_HP.hh>
// Absorption
#include <G4StoppingPhysics.hh>
// Ions
#include <G4IonPhysics.hh>
// Step Limiter
#include <G4StepLimiterPhysics.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PhysicsList::B2PhysicsList() 
: G4VModularPhysicsList(){
  SetVerboseLevel(1);
  defaultCutValue = 1.0 * nm;//yk

  // ver = 1;

  // RegisterPhysics(new G4EmLivermorePhysics());//yk
  // RegisterPhysics(new G4EmPenelopePhysics());//yk


  // // Decay
  // RegisterPhysics(new G4DecayPhysics());
  // RegisterPhysics(new G4RadioactiveDecayPhysics());

  // // Hadron process
  // RegisterPhysics(new G4HadronElasticPhysicsHP());
  // RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP());

  // // Absorption
  // RegisterPhysics(new G4StoppingPhysics());

  // // Ions
  // RegisterPhysics(new G4IonPhysics());

  // // Step Limiter
  // RegisterPhysics(new G4StepLimiterPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PhysicsList::~B2PhysicsList()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2PhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
}  
