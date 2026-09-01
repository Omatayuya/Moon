#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <array>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

struct ParticleData{
	TString name; //particle name for ROOT file
	TString title; //particle name for histogram etc.
	Double_t minBin; //histogram min energy (MeV)
	Double_t maxBin; //histogram max energy (MeV)
	Int_t nBins; //number of histogram bins
};

int toRootNeutron_angle(TString strShape = "box_10x10x1m", TString strTarget = "H_0ppm", TString strDepth = "depth_0.5m"){
	gROOT->Reset();
	gROOT->SetBatch();

	/*-----------------------------------*/
	//TString strTarget("H_0ppm");
	TString strPhi("380MV");
	vector<TString> vStrPrimary{"Proton", "Helium"};
	const UInt_t nGCR = vStrPrimary.size(); //number of sources (max=2: Proton, Helium)

	constexpr Double_t sensRadius = 2.5; //sensitive area radius (m)
	constexpr Double_t sensArea = M_PI*sensRadius*100*sensRadius*100; //sensitive surface area (cm^2)
	constexpr Double_t irrArea = 2000*2000; //irradiation surface area (cm^2) for Proton, Helium

	constexpr Int_t nTheta = 18; //number of bins for theta (0-90 deg)

	vector<ParticleData> vParticle{
	//	{"incident", "Incident", 10,   1e+7, 1000},
	//	{"proton",   "Proton",   1,    1e+4, 1000},
		{"neutron",  "Neutron",  1e-9, 1e+4, 1000},
	//	{"gamma",    "Gamma",    1e-2, 1e+4, 1000},
	//	{"deuteron", "Deuteron", 1,    1e+4, 1000},
	//	{"triton",   "Triton",   1,    1e+4, 1000},
	//	{"alpha",    "Alpha",    1,    1e+4, 1000},
	//	{"electron", "e-",       1e-2, 1e+4, 1000},
	//	{"positron", "e+",       1e-2, 1e+4, 1000},
	//	{"posiMuon", "mu+",      1,    1e+4, 1000},
	//	{"negaMuon", "mu-",      1,    1e+4, 1000},
	//	{"posiPion", "pi+",      1,    1e+4, 1000},
	//	{"negaPion", "pi-",      1,    1e+4, 1000},
	//	{"others",   "Others",   1e-9, 1e+4, 1000},
	};

	/*-----------------------------------*/

	//GCR flux
	vector<Double_t> vGcrFlux; //GCR 2pi flux (cm^-2 s^-1) @ 10MeV/n-100GeV/n
	constexpr Double_t gcrFactor = 4.; //GCR normalization factor
	if(strPhi=="380MV"){
		vGcrFlux.emplace_back(5.37992/gcrFactor); //1-H
		vGcrFlux.emplace_back(0.48827/gcrFactor); //2-He
	}else if(strPhi=="1000MV"){
		vGcrFlux.emplace_back(1.76038/gcrFactor); //1-H
		vGcrFlux.emplace_back(0.20587/gcrFactor); //2-He
	}else{ cerr << "--- Phi error" << endl; return 1; }
	if(vGcrFlux.size()!=nGCR){ cerr << "--- GCR nuclide number error" << endl; return 1; }

	constexpr Double_t binMin = 1e-9; //min. energy (MeV)
	constexpr Double_t binMax = 1e+4; //max. energy (MeV)
	constexpr Int_t nBins = 130;
	const Double_t binWidthLog = (log10(binMax)-log10(binMin))/nBins;
	vector<Double_t> vEdge;
	for(Int_t i=0; i<(nBins+1); i++) vEdge.emplace_back(pow(10,log10(binMin)+binWidthLog*i));

	vector<TH1D*> vHECurrent;
	for(Int_t i=0; i<nTheta; i++){
		stringstream ssName;
		ssName << "neutron_E*Current_" << i;

		stringstream ssTitle;
		ssTitle << "Neutron energy (" << i*90/nTheta << "-" << (i+1)*90/nTheta << "deg); Energy (MeV); Energy #times Current (cm^{-2} s^{-1})";
		
		vHECurrent.emplace_back(new TH1D(ssName.str().c_str(), ssTitle.str().c_str(), nBins, &vEdge.front()));
	}


	//read data
	for(size_t i=0; i<vStrPrimary.size(); i++){
		stringstream ssInput;
		ssInput << "/extra1/kusano/G4work/SELENE-R/EmissionDistribution/H_block/" << strShape << "/" << strTarget << "/" << strDepth << "/" << vStrPrimary.at(i) << "/" << strPhi << "/results.root";
		TString filename = ssInput.str();
		ifstream ifs(filename);
		if(ifs.fail()){ cerr << "--- File " << filename << " not found" << endl; return 1; }
		auto fin = new TFile(filename);

		//particle at a depth
		TString strBranch("surface");
		TTreeReader reader(strBranch, fin);
		const Int_t nEvents = reader.GetEntries(true);

		Double_t cntToCurrent = 0;
		if(i<nGCR) cntToCurrent = 1./sensArea/(static_cast<Double_t>(nEvents)/(vGcrFlux.at(i)*irrArea)); //counts -> counts/cm2/s
		//else cntToCurrent = 1./sensArea/(nEvents/(log(2)/vNatIsot.at(i-nGCR).halfLife*natIsotVolume*targetDensity*vNatIsot.at(i-nGCR).wFrac*avogadroConst/vNatIsot.at(i-nGCR).massNb*vNatIsot.at(i-nGCR).natAbundance)); //counts -> counts/cm2/s

		vector< TTreeReaderArray<double> > vaReaderEnergy, vaReaderMomentumZ, vaReaderPosX, vaReaderPosY;
		for(const auto& el : vParticle){
			TTreeReaderArray<double> aEnergy(reader, el.name+"Energy");
			vaReaderEnergy.emplace_back(std::move(aEnergy));
			
			TTreeReaderArray<double> aMomentumZ(reader, el.name+"EmissionAngle");
			vaReaderMomentumZ.emplace_back(std::move(aMomentumZ));
			
			TTreeReaderArray<double> aPosX(reader, el.name+"EmissionPositionX");
			vaReaderPosX.emplace_back(std::move(aPosX));
			
			TTreeReaderArray<double> aPosY(reader, el.name+"EmissionPositionY");
			vaReaderPosY.emplace_back(std::move(aPosY));
		}
		
		while(reader.Next()){
			for(size_t i=0; i<vParticle.size(); i++){
				for(const auto& el : vaReaderEnergy[i]){
					if(el==0) break;

					size_t idx = &el-&vaReaderEnergy[i][0];
					Double_t momZ = vaReaderMomentumZ[i][idx];
					
					Double_t posX = vaReaderPosX[i][idx];
					Double_t posY = vaReaderPosY[i][idx];
					Double_t radius = TMath::Sqrt(posX*posX+posY*posY);
					if(radius>sensRadius) continue;

					for(Int_t j=0; j<nTheta; j++){
						Double_t theta = TMath::ACos(momZ)*TMath::RadToDeg();
						if(theta>=1.*j*90/nTheta && theta<1.*(j+1)*90/nTheta) vHECurrent.at(j)->Fill(el, cntToCurrent);
					}
				}
			}
		}
	}

	//Double_t flux = 0;
	//for(const auto& el : vHECurrent) flux += el->Integral();
	//cout << flux << endl;

	//root file
	if(true){
		stringstream ssRoot;
		ssRoot << "Neutron_" << strShape << "_" << strDepth << "_" << strPhi << "_" << strTarget << "_angle.root";
		TString fRootOut = ssRoot.str();
		auto fOut = new TFile(fRootOut, "recreate");

		for(const auto& el : vHECurrent) el->Write();

		fOut->Close();
		cout << "--- " << fRootOut << " was created." << endl;
	}

	return 0;
}

