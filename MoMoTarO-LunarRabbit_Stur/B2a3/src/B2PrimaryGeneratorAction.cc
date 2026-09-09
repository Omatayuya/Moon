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
/// \file B2PrimaryGeneratorAction.cc
/// \brief Implementation of the B2PrimaryGeneratorAction class

#include "B2PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh" //20211222 HN
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh" // G4RunManager::GetCurrentRun() を使うため
#include "G4Event.hh"	   // G4Event::GetEventID() を使うため
#include "G4Run.hh"		   // G4Run::GetRunID() を使うため
#include "G4PhysicalConstants.hh"
#include "G4AnalysisManager.hh"
#include <cstdlib>

#include "Randomize.hh"
#include <fstream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::B2PrimaryGeneratorAction()
	: fParticleGun(nullptr), fpParticleGPS(nullptr)
{
	fParticleGun = new G4ParticleGun();
	fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("neutron"));

	fpParticleGPS = new G4GeneralParticleSource(); // 20211222 HN

	const char *envSpectrumFile = std::getenv("MOONNEUTRON_SPECTRUM_FILE");
	G4String spectrumFile = envSpectrumFile ? G4String(envSpectrumFile) : "/home/yomata/work/materials/Kusano/Code/Response/leakage/Neutron_380MV_H_0ppm_angle.root";
	LoadSpectrum(spectrumFile);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2PrimaryGeneratorAction::~B2PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fpParticleGPS; // 20211222 HN
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2PrimaryGeneratorAction::LoadSpectrum(const G4String &filename)
{
	TFile fin(filename.c_str());
	if (fin.IsZombie())
	{
		G4Exception("B2PrimaryGeneratorAction::LoadSpectrum", "MyCode0001", FatalException, ("Cannot open " + filename).c_str());
	}

	fEnergyEdges.clear();
	fEnergyCDF.clear();
	fTotalFlux = 0.0;
	std::vector<G4double> angleIntegral;

	G4int j = 0;
	while (true)
	{
		std::stringstream ssName;
		ssName << "neutron_E*Current_" << j;
		auto h = dynamic_cast<TH1D *>(fin.Get(ssName.str().c_str()));
		if (!h)
			break;

		const G4int nBins = h->GetNbinsX();
		if (fEnergyEdges.empty())
		{
			fEnergyEdges.resize(nBins + 1);
			for (G4int i = 1; i <= nBins + 1; ++i)
				fEnergyEdges.at(i - 1) = h->GetXaxis()->GetBinLowEdge(i);
		}

		std::vector<G4double> cdf(nBins + 1);
		cdf.at(0) = 0.;
		for (G4int i = 1; i <= nBins; ++i)
		{
			cdf.at(i) = cdf.at(i - 1) + h->GetBinContent(i);
		}
		angleIntegral.emplace_back(cdf.back());
		if (cdf.back() > 0)
		{
			for (auto &el : cdf)
				el /= cdf.back();
		}
		fEnergyCDF.emplace_back(std::move(cdf));

		++j;
	}
	fin.Close();

	fNTheta = j;
	fAngleCumWeight.clear();

	for (const auto &el : angleIntegral)
		fTotalFlux += el;

	fAngleCumWeight.resize(fNTheta);
	G4double running = 0;
	for (G4int k = 0; k < fNTheta; ++k)
	{
		running += angleIntegral.at(k);
		fAngleCumWeight.at(k) = running / fTotalFlux;
	}
}

void B2PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
	G4double worldZHalfLength = 0;
	G4LogicalVolume *worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
	G4Box *worldBox = NULL;
	if (worldLV)
		worldBox = dynamic_cast<G4Box *>(worldLV->GetSolid());
	if (worldBox)
		worldZHalfLength = worldBox->GetZHalfLength();
	else
	{
		G4cerr << "World volume of box not found." << G4endl;
		G4cerr << "Perhaps you have changed geometry." << G4endl;
		G4cerr << "The gun will be place in the center." << G4endl;
	}

	// 1. 角度ビンを重み付きで選ぶ
	G4double u = G4UniformRand();
	G4int j = 0;
	while (j < fNTheta - 1 && u > fAngleCumWeight.at(j))
		++j;

	// 2. 選ばれたビンのエネルギーCDFから逆引き(ビン内線形補間)
	const auto &cdf = fEnergyCDF.at(j);
	G4double v = G4UniformRand();
	auto it = std::upper_bound(cdf.begin(), cdf.end(), v);
	size_t i = std::distance(cdf.begin(), it);
	if (i == 0)
		i = 1;
	if (i >= cdf.size())
		i = cdf.size() - 1;
	--i;

	G4double eLow = fEnergyEdges.at(i), eHigh = fEnergyEdges.at(i + 1);
	G4double frac = (cdf.at(i + 1) > cdf.at(i)) ? (v - cdf.at(i)) / (cdf.at(i + 1) - cdf.at(i)) : 0.5;
	G4double energy = eLow + frac * (eHigh - eLow);

	// 3. ビン内でthetaを一様サンプル
	G4double thetaBinWidth = 90. / fNTheta;
	G4double thetaDeg = j * thetaBinWidth + G4UniformRand() * thetaBinWidth;
	G4double theta = thetaDeg * deg;

	// 4. phiは一様
	G4double phi = twopi * G4UniformRand();

	G4ThreeVector direction(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));

	// 5. 位置
	G4double posX = (2 * G4UniformRand() - 1) * 3 * m;
	G4double posY = (2 * G4UniformRand() - 1) * 3 * m;
	G4double posZ = 0 * mm;

	fParticleGun->SetParticleEnergy(energy * MeV);
	fParticleGun->SetParticleMomentumDirection(direction);
	fParticleGun->SetParticlePosition(G4ThreeVector(posX, posY, posZ));
	fParticleGun->GeneratePrimaryVertex(anEvent);

	auto analysisManager = G4AnalysisManager::Instance();

	if (!fHistIdsReady)
	{
		const G4int nThetaBins = 18; // B2RunAction.cc の nThetaBins と必ず一致させる
		fHistIdEnergyByTheta.clear();
		for (G4int k = 0; k < nThetaBins; ++k)
		{
			G4String name = "hPrimEnergyByTheta_" + std::to_string(k);
			fHistIdEnergyByTheta.push_back(analysisManager->GetH1Id(name));
		}
		fHistIdEnergy = analysisManager->GetH1Id("hPrimEnergy");
		fHistIdTheta = analysisManager->GetH1Id("hPrimTheta");
		fHistIdPhi = analysisManager->GetH1Id("hPrimPhi");
		fHistIdPosZ = analysisManager->GetH1Id("hPrimPosZ");
		fHistIdPosXY = analysisManager->GetH2Id("hPrimPosXY");
		fHistIdsReady = true;
	}

	const G4int nThetaBinsOut = static_cast<G4int>(fHistIdEnergyByTheta.size());
	G4int thetaBinOut = static_cast<G4int>(thetaDeg / (90.0 / nThetaBinsOut));
	if (thetaBinOut >= nThetaBinsOut)
		thetaBinOut = nThetaBinsOut - 1;
	if (thetaBinOut < 0)
		thetaBinOut = 0;

	analysisManager->FillH1(fHistIdEnergyByTheta.at(thetaBinOut), energy);
	analysisManager->FillH1(fHistIdEnergy, energy);
	analysisManager->FillH1(fHistIdTheta, thetaDeg);
	analysisManager->FillH1(fHistIdPhi, phi / deg);
	analysisManager->FillH1(fHistIdPosZ, posZ / mm);
	analysisManager->FillH2(fHistIdPosXY, posX / mm, posY / mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
