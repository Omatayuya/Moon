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
/// \file B2TrackerHit.cc
/// \brief Implementation of the B2TrackerHit class

#include "B2TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh" // GetCurrentRun() を使うため
#include "G4Run.hh"		   // GetRunID() を使うため

#include <iomanip>
#include <fstream>

G4ThreadLocal G4Allocator<B2TrackerHit> *B2TrackerHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerHit::B2TrackerHit()
	: G4VHit(),
	  fEventID(-1),
	  fTrackID(-1),
	  fChamberNb(-1),
	  fEdep(0.),
	  fStepLength(0.),
	  fGTime(0.),
	  fLTime(0.),
	  fPTime(0.),
	  fDTime(0.),
	  fPrePos(G4ThreeVector()),
	  fPostPos(G4ThreeVector()),
	  fPID(-1),
	  fPPID(-1),
	  fRunID(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerHit::~B2TrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerHit::B2TrackerHit(const B2TrackerHit &right)
	: G4VHit()
{
	fTrackID = right.fTrackID;
	fChamberNb = right.fChamberNb;
	fEdep = right.fEdep;
	fStepLength = right.fStepLength;
	fGTime = right.fGTime;
	fLTime = right.fLTime;
	fPTime = right.fPTime;
	fDTime = right.fDTime;
	fPrePos = right.fPrePos;
	fPostPos = right.fPostPos;
	fEventID = right.fEventID;
	fPID = right.fPID;
	fPPID = right.fPPID;
	fPname = right.fPname;
	fPreProc = right.fPreProc;
	fPostProc = right.fPostProc;
	fCProc = right.fCProc;
	fCModel = right.fCModel;
	fPreKinE = right.fPreKinE;
	fPostKinE = right.fPostKinE;
	fLLVname = right.fLLVname;
	fRunID = right.fRunID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const B2TrackerHit &B2TrackerHit::operator=(const B2TrackerHit &right)
{
	fTrackID = right.fTrackID;
	fChamberNb = right.fChamberNb;
	fEdep = right.fEdep;
	fStepLength = right.fStepLength;
	fGTime = right.fGTime;
	fLTime = right.fLTime;
	fPTime = right.fPTime;
	fDTime = right.fDTime;
	fPrePos = right.fPrePos;
	fPostPos = right.fPostPos;
	fEventID = right.fEventID;
	fPID = right.fPID;
	fPPID = right.fPPID;
	fPname = right.fPname;
	fPreProc = right.fPreProc;
	fPostProc = right.fPostProc;
	fCProc = right.fCProc;
	fCModel = right.fCModel;
	fPreKinE = right.fPreKinE;
	fPostKinE = right.fPostKinE;
	fLLVname = right.fLLVname;
	fRunID = right.fRunID;

	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B2TrackerHit::operator==(const B2TrackerHit &right) const
{
	return (this == &right) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerHit::Draw()
{
	G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance();
	if (pVVisManager)
	{
		G4Circle circle(fPostPos);
		circle.SetScreenSize(4.);
		circle.SetFillStyle(G4Circle::filled);
		G4Colour colour(1., 0., 0.);
		G4VisAttributes attribs(colour);
		circle.SetVisAttributes(attribs);
		pVVisManager->Draw(circle);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerHit::Print()
{
	G4String datapath = "../log/output.csv";
	std::ofstream outFile("output.csv", std::ios::app);

	// Space should not exist.
	if (fCModel == "Binary Cascade")
	{
		fCModel = "BinaryCascade";
	}
	if (fCModel == "Binary Light Ion Reaction")
	{
		fCModel = "BinaryLightIonReaction";
	}

	// if (fCProc == "phot" || fCProc == "compt") {

	fRunID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();

	outFile
		// << fRunID
		// << "," << fEventID // yk
		//<< "eventID: "
		<< fEventID
		//<< "\ttrackID: " << fTrackID
		//<< "\tchamberNb: "
		<< ",target" << fChamberNb
		//<< " PID: "
		<< "," << fPID // yk
		//<< " PPID: "     //ParentParticleID
		<< "," << fPPID // yk
		//<< " Particle: "
		<< "," << fPname // yk
		//<< " PreProc: "
		<< "," << fPreProc // yk
		//<< " PostProc: "
		<< "," << fPostProc // yk
		//<< " CProc: "
		<< "," << fCProc // yk
		//<< " CModel: "
		<< "," << fCModel // yk
		//<< "\tVerVol: " << fLLVname //yk
		//<< " PreKinE: " << std::setw(4) << G4BestUnit(fPreKinE, "Energy") //yk
		//<< " PreKinE: "
		<< "," << std::scientific << fPreKinE * 1000 * 1000 // yk [eV]
		//<< " PostKinE: "
		<< "," << std::scientific << fPostKinE * 1000 * 1000 // yk [eV]
		//<< "\tEdep: " << std::setw(4) << G4BestUnit(fEdep, "Energy")
		//<< " Edep: "
		<< "," << std::scientific << fEdep * 1000 * 1000 // yk [eV]
		//<< " GlobalT: "
		<< "," << std::scientific << fGTime // yk [ns]
		//<< " LocalT: "
		<< "," << std::scientific << fLTime // yk [ns]
		//<< " ProperT: "
		<< "," << std::scientific << fPTime // yk [ns]
		//<< " DeltaT: "
		<< "," << std::scientific << fDTime // yk [ns]
											//<< "\tDiff: " << std::setw(4) << G4BestUnit(fPreKinE - fKinE - fEdep, "Energy")
											//<< " PrePos: " << std::setw(4) << G4BestUnit(fPrePos, "Length") //yk
		//<< " PrePos: "
		<< "," << std::scientific << fPrePos[0] / 1000 // yk [m]
		<< "," << fPrePos[1] / 1000					   // yk [m]
		<< "," << fPrePos[2] / 1000					   // yk [m]
													   //<< " PostPos: " << std::setw(4) << G4BestUnit(fPostPos, "Length")
		//<< " PostPos: "
		<< "," << std::scientific << fPostPos[0] / 1000 // yk [m]
		<< "," << fPostPos[1] / 1000					// yk [m]
		<< "," << fPostPos[2] / 1000					// yk [m]
		<< "," << fStepLength / 1000                    // yk [m]
		<< G4endl;
	//}
	//}
	outFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
