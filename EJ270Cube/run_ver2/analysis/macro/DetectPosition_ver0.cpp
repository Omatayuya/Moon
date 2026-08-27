#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm> // std::max
#include <set>
#include <map>
#include <tuple>
#include <sstream>
#include <limits>

using namespace std;

#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TBox.h>
#include <TFile.h>
#include <cctype> // std::isspace
#include <TPaveText.h>
#include <TMultiGraph.h>
#include <TLegend.h>

void DetectPosition_ver0()
{
    vector<TString> folder = {"0ppm", "10ppm", "20ppm", "50ppm", "100ppm", "200ppm", "500ppm", "1000ppm", "2000ppm", "5000ppm", "10000ppm"};
    // vector<TString> folder = {"0ppm"};

    for (int folderID = 0; folderID < folder.size(); folderID++)
    {

        // ==================================================================

        const string TargetParticle = "neutron";
        // const Double_t irrArea = pow(450, 2);                    // irradiation surface area (cm^2) for Proton, Helium
        const Double_t DetectorOffsetZ = 400; // Detector offset in Z (mm)
        const Double_t irrArea = 600 * 600;   // irradiation surface area (cm^2) for Proton, Helium

        // energy window
        constexpr double scatterEdepLow = 1.0; // MeV
        // constexpr double scatterEdepHigh = 3.0; // MeV
        constexpr double captureEdepLow = 4.5; // MeV

        // Thermal neutron cut (109Cd)
        constexpr double TNEnergyCut = 5e-7; // MeV

        for (auto axis : {"X", "Y", "Z"})
        {
            gStyle->SetLabelFont(62, axis);
            gStyle->SetTitleFont(62, axis);
            gStyle->SetTitleOffset(1.2, axis); // 軸タイトルのオフセット
            gStyle->SetLabelSize(0.04, axis);  // 目盛り数字のサイズ
            gStyle->SetTitleSize(0.04, axis);  // 軸タイトルのサイズ
        }
        gStyle->SetTextFont(62);
        gStyle->SetTitleFont(62, "");

        gStyle->SetPadGridX(true);
        gStyle->SetPadGridY(true);
        // gStyle->SetPalette(kRainBow);
        gStyle->SetOptStat(0);

        gStyle->SetHistLineWidth(2);

        // ==================================================================

        TFile *fin = TFile::Open("../../" + folder[folderID] + "/results.root");
        if (!fin || fin->IsZombie())
        {
            cout << "Failed to open results.root" << endl;
            return;
        }

        TTree *HitTree = (TTree *)fin->Get("Hit");

        /*Input ParticleData*/
        TTree *RunInfoTree = (TTree *)fin->Get("RunInfo");
        double moonNeutronFlux = 0.0;
        RunInfoTree->SetBranchAddress("TotalFlux", &moonNeutronFlux);
        RunInfoTree->GetEntry(0);
        const Double_t MoonNeutronFlux = moonNeutronFlux;
        cout << "Moon neutron flux: " << MoonNeutronFlux << " cm^-2 s^-1" << endl;

        constexpr int nThetaBins = 18;        // 0-90度を5度刻みで18分割
        constexpr double thetaBinWidth = 5.0; // deg

        vector<TH1F *> vH_ip_theta(nThetaBins);
        for (int k = 0; k < nThetaBins; ++k)
        {
            TString hname = Form("hPrimEnergyByTheta_%d", k);
            vH_ip_theta[k] = (TH1F *)fin->Get(hname);
            if (!vH_ip_theta[k])
            {
                cout << "Histogram " << hname << " not found!" << endl;
                return;
            }
            vH_ip_theta[k]->SetDirectory(nullptr);
            vH_ip_theta[k]->SetStats(0);
            vH_ip_theta[k]->SetLineWidth(2);
            vH_ip_theta[k]->SetLineColor(TColor::GetColorPalette(
                static_cast<int>(1.0 * k * (TColor::GetNumberOfColors() - 1) / (nThetaBins - 1))));
        }

        TH1F *hPrimEnergy = (TH1F *)fin->Get("hPrimEnergy");
        int nEvents = hPrimEnergy->GetEntries();
        const Double_t eqTime = nEvents / (irrArea * MoonNeutronFlux);

        cout << "Total number of events in InputDataTree: " << nEvents << endl;
        cout << "Equivalent time: " << eqTime << " s" << endl;

        // InputDataTree->BuildIndex("eventID");
        Double_t cntToCurrent = MoonNeutronFlux / nEvents; // Convert to s^-1 cm^-2
        for (int k = 0; k < nThetaBins; ++k)
        {
            vH_ip_theta[k]->Scale(cntToCurrent);
        }

        /* Output Data*/
        struct ChamberEventData
        {
            double edepSum = 0.0;
            double cpEdepSum = 0.0;
            double scEdepSum = 0.0;
            double cpPosX = 0.0;
            double cpPosY = 0.0;
            double cpPosZ = 0.0;
            double scPosX = 0.0;
            double scPosY = 0.0;
            double scPosZ = 0.0;
            double cpTriggerTime = numeric_limits<double>::infinity();
            double scTriggerTime = numeric_limits<double>::infinity();
            bool captureflag = false;
            bool scatterflag = false;
        } chamberEventData;

        struct EventChamberID
        {
            int eventID = 0;
            double primEnergy = 0.0;
            string chamberNb;

            bool operator<(const EventChamberID &other) const
            {
                if (eventID != other.eventID)
                    return eventID < other.eventID;
                return chamberNb < other.chamberNb;
            }
        } eventChamberID;

        map<EventChamberID, ChamberEventData> DetectorchamberMap;

        TTree *DetectorchamberTree = new TTree("DetectorchamberTree", "event-chamber-level tree from csvfile");
        DetectorchamberTree->Branch("eventID", &eventChamberID.eventID, "eventID/I");
        DetectorchamberTree->Branch("primEnergy", &eventChamberID.primEnergy, "primEnergy/D");
        DetectorchamberTree->Branch("chamberNb", &eventChamberID.chamberNb);
        DetectorchamberTree->Branch("edepSum", &chamberEventData.edepSum, "edepSum/D");
        DetectorchamberTree->Branch("captureEdepSum", &chamberEventData.cpEdepSum, "captureEdepSum/D");
        DetectorchamberTree->Branch("scatterEdepSum", &chamberEventData.scEdepSum, "scatterEdepSum/D");
        DetectorchamberTree->Branch("capturePosX", &chamberEventData.cpPosX, "capturePosX/D");
        DetectorchamberTree->Branch("capturePosY", &chamberEventData.cpPosY, "capturePosY/D");
        DetectorchamberTree->Branch("capturePosZ", &chamberEventData.cpPosZ, "capturePosZ/D");
        DetectorchamberTree->Branch("scatterPosX", &chamberEventData.scPosX, "scatterPosX/D");
        DetectorchamberTree->Branch("scatterPosY", &chamberEventData.scPosY, "scatterPosY/D");
        DetectorchamberTree->Branch("scatterPosZ", &chamberEventData.scPosZ, "scatterPosZ/D");
        DetectorchamberTree->Branch("captureTriggerTime", &chamberEventData.cpTriggerTime, "captureTriggerTime/D");
        DetectorchamberTree->Branch("scatterTriggerTime", &chamberEventData.scTriggerTime, "scatterTriggerTime/D");
        DetectorchamberTree->Branch("captureflag", &chamberEventData.captureflag, "captureflag/O");
        DetectorchamberTree->Branch("scatterflag", &chamberEventData.scatterflag, "scatterflag/O");
        DetectorchamberTree->SetDirectory(nullptr);

        double primEnergy;
        int eventID, fPID, fPPID;
        char chamberNb[256], fpname[256], fCProc[256];
        char collection[256], fPreProc[256], fPostProc[256], fCModel[256];
        double fPreKinE, fPostKinE, fEdep, fGTime, fLTime, fPTime, fDTime;
        double fPrePosX, fPrePosY, fPrePosZ, fPostPosX, fPostPosY, fPostPosZ, fStepLength;

        HitTree->SetBranchAddress("eventID", &eventID);
        HitTree->SetBranchAddress("primEnergy", &primEnergy);
        HitTree->SetBranchAddress("chamberNb", chamberNb);
        HitTree->SetBranchAddress("collection", collection);
        HitTree->SetBranchAddress("PID", &fPID);
        HitTree->SetBranchAddress("PPID", &fPPID);
        HitTree->SetBranchAddress("Pname", fpname);
        HitTree->SetBranchAddress("PreProc", fPreProc);
        HitTree->SetBranchAddress("PostProc", fPostProc);
        HitTree->SetBranchAddress("CProc", fCProc);
        HitTree->SetBranchAddress("CModel", fCModel);
        HitTree->SetBranchAddress("PreKinE", &fPreKinE);
        HitTree->SetBranchAddress("PostKinE", &fPostKinE);
        HitTree->SetBranchAddress("Edep", &fEdep);
        HitTree->SetBranchAddress("GTime", &fGTime);
        HitTree->SetBranchAddress("LTime", &fLTime);
        HitTree->SetBranchAddress("PTime", &fPTime);
        HitTree->SetBranchAddress("DTime", &fDTime);
        HitTree->SetBranchAddress("PrePosX", &fPrePosX);
        HitTree->SetBranchAddress("PrePosY", &fPrePosY);
        HitTree->SetBranchAddress("PrePosZ", &fPrePosZ);
        HitTree->SetBranchAddress("PostPosX", &fPostPosX);
        HitTree->SetBranchAddress("PostPosY", &fPostPosY);
        HitTree->SetBranchAddress("PostPosZ", &fPostPosZ);
        HitTree->SetBranchAddress("StepLength", &fStepLength);

        const vector<string> captureStepParticles = {"alpha", "triton"};
        // const vector<string> scatterStepParticles = {"proton", "C12", "O16", "N14"};
        const vector<string> scatterStepParticles = {"proton"};

        int nHits = HitTree->GetEntries();
        for (int i = 0; i < nHits; ++i)
        {
            HitTree->GetEntry(i);

            // double fEdepQ = 0.0;     // Initialize quenched energy deposit
            // const double kB = 0.012; // Birks' constant (mm/MeV)
            // const double S = 1.0;    // Scintillation efficiency

            // if (fStepLength > 0. && fEdep > 0.)
            // {
            //     fStepLength = fStepLength * 1000;  // Convert m to mm
            //     double dedx = fEdep / fStepLength; // [MeV/mm]
            //     double quenchingFactor = 1. / (1. + kB * dedx);
            //     fEdepQ = S * fEdep * quenchingFactor;
            // }

            // const double preKinEMeV = fPreKinE / 1e6; // Convert eV to MeV
            // const double edepMeV = fEdep / 1e6;       // Convert eV to MeV

            if (string(chamberNb) == "target1")
            {
                EventChamberID eventChamberID{eventID, primEnergy, string(chamberNb)};
                auto &acc = DetectorchamberMap[eventChamberID];
                acc.edepSum += fEdep;

                /* screening conditions for Capture & Scatter event*/
                const bool isCaptureParticle =
                    find(captureStepParticles.begin(), captureStepParticles.end(), fpname) != captureStepParticles.end();
                const bool isScatterParticle =
                    find(scatterStepParticles.begin(), scatterStepParticles.end(), fpname) != scatterStepParticles.end();

                const bool isCaptureDepositStep =
                    (fEdep > 0.0) &&
                    isCaptureParticle &&
                    (string(fCProc) == "neutronInelastic");
                const bool isScatterDepositStep =
                    (fEdep > 0.0) &&
                    isScatterParticle &&
                    // isPrimaryNeutronSecondary &&
                    (string(fCProc) == "hadElastic");

                if (isCaptureDepositStep)
                {
                    acc.cpEdepSum += fEdep;

                    if (acc.captureflag == false && acc.cpEdepSum > captureEdepLow)
                    {
                        acc.captureflag = true;
                        acc.cpTriggerTime = min(acc.cpTriggerTime, fGTime);
                        acc.cpPosX = fPrePosX;
                        acc.cpPosY = fPrePosY;
                        acc.cpPosZ = fPrePosZ;
                    }
                }

                if (isScatterDepositStep)
                {
                    acc.scEdepSum += fEdep;
                    if (acc.scatterflag == false && acc.scEdepSum > scatterEdepLow)
                    {
                        acc.scatterflag = true;
                        acc.scTriggerTime = min(acc.scTriggerTime, fGTime);
                        acc.scPosX = fPrePosX;
                        acc.scPosY = fPrePosY;
                        acc.scPosZ = fPrePosZ;
                    }
                }
            }
        }

        for (const auto &entry : DetectorchamberMap)
        {
            eventChamberID = entry.first;
            chamberEventData = entry.second;
            DetectorchamberTree->Fill();
        }
        DetectorchamberMap.clear();
        const int entries = DetectorchamberTree->GetEntries();

        int BinWidthZ = 5; // mm
        int minZ = 0;      // mm
        int maxZ = 200;    // mm
        int nBinsZ = (maxZ - minZ) / BinWidthZ;
        int BinWidthXY = 5; // mm
        int minXY = -50;    // mm
        int maxXY = 50;     // mm
        int nBinsXY = (maxXY - minXY) / BinWidthXY;

        // Match vH_ip_theta's log-uniform energy binning (see B2RunAction.cc)
        constexpr int nEnergyBins = 500;
        constexpr double energyMin = 1e-10; // MeV
        constexpr double energyMax = 1e4;   // MeV
        double energyBins[nEnergyBins + 1];
        for (int i = 0; i <= nEnergyBins; ++i)
            energyBins[i] = energyMin * pow(energyMax / energyMin, static_cast<double>(i) / nEnergyBins);

        TH1F *h1_cpposZ = new TH1F("h1_cpposZ", Form("Capture Position Z ;Z (mm);Count rate (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        TH1F *h1_cpposZ_TNcut = new TH1F("h1_cpposZ_TNcut", Form("Capture Position Z  (TN cut / %.1f eV <);Z (mm);Count rate (s^{-1} %d mm^{-1})", TNEnergyCut * 1e6, BinWidthZ), nBinsZ, minZ, maxZ);
        TH2D *h2_cpposZ = new TH2D("h2_cpposZ", Form("Capture Position Z; Energy (MeV); Z (mm); Count rate"), nEnergyBins, energyBins, nBinsZ, minZ, maxZ);
        TH2D *h2_cpposXY = new TH2D("h2_cpposXY", Form("Capture Position XY Distribution;X (mm);Y (mm);Count rate (s^{-1} %d #times %d mm^{-2})", BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        TH2D *h2_cpposXY_TNcut = new TH2D("h2_cpposXY_TNcut", Form("Capture Position XY Distribution (TN cut / %.1f eV <);X (mm);Y (mm);Count rate (s^{-1} %d#times%d mm^{-2})", TNEnergyCut * 1e6, BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        TH1F *h1_scposZ = new TH1F("h1_scposZ", Form("Scatter Position Z Distribution;Z (mm);Count rate (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        TH2D *h2_scposXY = new TH2D("h2_scposXY", Form("Scatter Position XY Distribution;X (mm);Y (mm);Count rate (s^{-1} %d#times%d mm^{-2})", BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);

        constexpr int nEBins = 4;
        const double primEnergyEdges[nEBins + 1] = {0.0, 5e-7, 1e-3, 1.0, numeric_limits<double>::infinity()};
        const TString primEnergyLabels[nEBins] = {"E < 0.5 eV", "0.5 eV #leq E < 1 keV", "1 keV #leq E < 1 MeV", "E #geq 1 MeV"};
        const int primEnergyColors[nEBins] = {kOrange + 8, kGreen - 7, kGreen + 2, kBlue};
        vector<TH1F *> vH_cpposZ_byE(nEBins);
        for (int e = 0; e < nEBins; ++e)
        {
            TString hname = Form("h1_cpposZ_E%d", e);
            vH_cpposZ_byE[e] = new TH1F(hname,
                                        Form("Capture Position Z (%s);Z (mm);Count rate (s^{-1} %d mm^{-1})", primEnergyLabels[e].Data(), BinWidthZ),
                                        nBinsZ, minZ, maxZ);
        }

        // Z region bins for capture XY histograms
        vector<double> zRegionEdges = {0, 40, 160, 200};
        const int nZRegions = zRegionEdges.size() - 1;
        vector<TString> zRegionLabels;
        for (int z = 0; z < nZRegions; ++z)
        {
            zRegionLabels.push_back(Form("%.0f <= Z < %.0f ", zRegionEdges[z], zRegionEdges[z + 1]));
        }
        vector<TH2D *> vH2_cpposXY_byZ(nZRegions);
        for (int z = 0; z < nZRegions; ++z)
        {
            TString hname = Form("h2_cpposXY_Z%d", z);
            vH2_cpposXY_byZ[z] = new TH2D(hname,
                                          Form("Capture Position XY Distribution (%s);X (mm);Y (mm);Count rate (s^{-1} %d#times%d mm^{-2})",
                                               zRegionLabels[z].Data(), BinWidthXY, BinWidthXY),
                                          nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        }
        constexpr double xFidCutLow = -30.0; // mm, フィデューシャルカット下限
        constexpr double xFidCutHigh = 30.0; // mm, フィデューシャルカット上限
        constexpr int nEnergyBins_cppos = 200;
        for (int i = 0; i <= nEnergyBins_cppos; ++i)
            energyBins[i] = energyMin * pow(energyMax / energyMin, static_cast<double>(i) / nEnergyBins_cppos);
        vector<TH1D *> vH1_cpposY_byZ(nZRegions);
        for (int z = 0; z < nZRegions; ++z)
        {
            TString hname = Form("h1_cpposY_Z%d", z);
            vH1_cpposY_byZ[z] = new TH1D(hname,
                                         Form("Capture Position Y (%s, %.0f < X < %.0f );Y (mm);Count rate (s^{-1} %d mm^{-1})",
                                              zRegionLabels[z].Data(), xFidCutLow, xFidCutHigh, BinWidthXY),
                                         nBinsXY, minXY, maxXY);
        }
        vector<TH2D *> vH2_cpposEY_byZ(nZRegions);
        for (int z = 0; z < nZRegions; ++z)
        {
            TString hname = Form("h2_cpposEY_Z%d", z);
            vH2_cpposEY_byZ[z] = new TH2D(hname,
                                          Form("Capture Position Y (%s, %.0f < X < %.0f );Energy (MeV);Y (mm);Count rate",
                                               zRegionLabels[z].Data(), xFidCutLow, xFidCutHigh),
                                          nEnergyBins_cppos, energyBins, nBinsXY, minXY, maxXY);
        }
        vector<vector<TH1D *>> vH1_cpposY_byZ_byE(nZRegions, vector<TH1D *>(nEBins));
        for (int z = 0; z < nZRegions; ++z)
        {
            for (int e = 0; e < nEBins; ++e)
            {
                TString hname = Form("h1_cpposY_Z%d_E%d", z, e);
                vH1_cpposY_byZ_byE[z][e] = new TH1D(hname,
                                                    Form("Capture Position Y (%s, %s);Y (mm);Count rate (s^{-1} %d mm^{-1})",
                                                         zRegionLabels[z].Data(), primEnergyLabels[e].Data(), BinWidthXY),
                                                    nBinsXY, minXY, maxXY);
                vH1_cpposY_byZ_byE[z][e]->SetLineColor(primEnergyColors[e]);
            }
        }

        using HistLegPair = pair<TH1 *, TString>;
        vector<vector<HistLegPair>> vvHist;
        vvHist.push_back({{h2_cpposZ, "Capture Position Z"}});
        vvHist.push_back({{h2_cpposXY, "Capture Position XY"}});
        for (int z = 0; z < nZRegions; ++z)
        {
            // vvHist.push_back({{vH2_cpposXY_byZ[z], Form("Capture Position XY (%s)", zRegionLabels[z].Data())}});
            vvHist.push_back({{vH2_cpposEY_byZ[z], Form("Capture Position EY (%s)", zRegionLabels[z].Data())}});

            vH1_cpposY_byZ[z]->SetLineColor(kBlack);
            vector<HistLegPair> group;
            group.push_back({vH1_cpposY_byZ[z], "Total"});
            for (int e = 0; e < nEBins; ++e)
                group.push_back({vH1_cpposY_byZ_byE[z][e], primEnergyLabels[e]});
            vvHist.push_back(group);
        }
        vvHist.push_back({{h1_scposZ, "Scatter Position Z"}});
        vvHist.push_back({{h2_scposXY, "Scatter Position XY"}});

        for (int i = 0; i < entries; ++i)
        {
            DetectorchamberTree->GetEntry(i);
            double capturePosZ = chamberEventData.cpPosZ - DetectorOffsetZ; // Convert m to mm and subtract detector offset
            double capturePosX = chamberEventData.cpPosX;                   // Convert m to mm
            double capturePosY = chamberEventData.cpPosY;                   // Convert m to mm
            if (chamberEventData.captureflag)
            {
                h1_cpposZ->Fill(capturePosZ);
                h2_cpposZ->Fill(eventChamberID.primEnergy, capturePosZ);
                h2_cpposXY->Fill(capturePosX, capturePosY);
                for (int z = 0; z < nZRegions; ++z)
                {
                    if (capturePosZ >= zRegionEdges[z] && capturePosZ < zRegionEdges[z + 1])
                    {
                        vH2_cpposXY_byZ[z]->Fill(capturePosX, capturePosY);
                        if (capturePosX > xFidCutLow && capturePosX < xFidCutHigh)
                        {
                            vH1_cpposY_byZ[z]->Fill(capturePosY);
                            vH2_cpposEY_byZ[z]->Fill(eventChamberID.primEnergy, capturePosY);
                            for (int e = 0; e < nEBins; ++e)
                            {
                                if (eventChamberID.primEnergy >= primEnergyEdges[e] && eventChamberID.primEnergy < primEnergyEdges[e + 1])
                                {
                                    vH1_cpposY_byZ_byE[z][e]->Fill(capturePosY);
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }

                if (eventChamberID.primEnergy > TNEnergyCut) // TN cut
                {
                    h1_cpposZ_TNcut->Fill(capturePosZ);
                    h2_cpposXY_TNcut->Fill(capturePosX, capturePosY);
                }

                // primEnergyで分岐してFill
                for (int e = 0; e < nEBins; ++e)
                {
                    if (eventChamberID.primEnergy >= primEnergyEdges[e] && eventChamberID.primEnergy < primEnergyEdges[e + 1])
                    {
                        vH_cpposZ_byE[e]->Fill(capturePosZ);
                        break;
                    }
                }
            }

            double scatterPosZ = chamberEventData.scPosZ - DetectorOffsetZ; // Convert m to mm and subtract detector offset
            double scatterPosX = chamberEventData.scPosX;                   // Convert m to mm
            double scatterPosY = chamberEventData.scPosY;                   // Convert
            if (chamberEventData.scatterflag)
            {
                h1_scposZ->Fill(scatterPosZ);
                h2_scposXY->Fill(scatterPosX, scatterPosY);
            }
        }

        // Draw
        vector<TCanvas *> vCan;
        {
            TCanvas *c1 = new TCanvas("c1", "Initial Particle Energy Distribution", 800, 600);
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetGridx();
            gPad->SetGridy();

            double yMax = 0;
            for (int k = 0; k < nThetaBins; ++k)
                yMax = max(yMax, vH_ip_theta[k]->GetMaximum());

            TLegend *leg = new TLegend(0.05, 0.7, 0.45, 0.9);
            leg->SetNColumns(4);
            for (int k = 0; k < nThetaBins; ++k)
            {
                if (k == 0)
                {
                    vH_ip_theta[k]->SetTitle(Form("Moon Neutron (Total current: %.4f cm^{-2} s^{-1} / eqTime: %.2f s);Energy (MeV); Energy #times Current (cm^{-2} s^{-1} bin^{-1})", MoonNeutronFlux, eqTime));
                    vH_ip_theta[k]->SetMaximum(yMax * 1.5);
                    vH_ip_theta[k]->SetMinimum(1e-5);
                    vH_ip_theta[k]->Draw("HIST");
                }
                else
                {
                    vH_ip_theta[k]->Draw("HIST SAME");
                }
                leg->AddEntry(vH_ip_theta[k],
                              Form("%d-%d deg", static_cast<int>(k * thetaBinWidth), static_cast<int>((k + 1) * thetaBinWidth)), "l");
            }
            leg->Draw();
            vCan.push_back(c1);
        }

        TCanvas *cCpposZByE = new TCanvas("cCpposZByE", "Capture Position Z by Primary Energy", 800, 600);
        gPad->SetGridx();
        gPad->SetGridy();
        double yMaxE = 0;
        for (int e = 0; e < nEBins; ++e)
        {
            vH_cpposZ_byE[e]->Scale(1.0 / eqTime); // 他のヒストグラムと同様にcps換算
        }
        TLegend *legE = new TLegend(0.55, 0.7, 0.88, 0.88);
        h1_cpposZ->Scale(1.0 / eqTime); // 他のヒストグラムと同様にcps換算
        h1_cpposZ->SetLineColor(kBlack);
        h1_cpposZ->SetMinimum(0);
        h1_cpposZ->Draw("HIST E");
        legE->AddEntry(h1_cpposZ, "Total", "l");
        for (int e = 0; e < nEBins; ++e)
        {
            vH_cpposZ_byE[e]->SetLineColor(primEnergyColors[e]);
            // vH_cpposZ_byE[e]->SetMaximum(yMaxE * 1.3);
            vH_cpposZ_byE[e]->Draw(e == 0 ? "HIST E SAME" : "HIST E SAME");
            legE->AddEntry(vH_cpposZ_byE[e], primEnergyLabels[e], "l");
        }
        legE->Draw();
        vCan.push_back(cCpposZByE);

        for (int i = 0; i < vvHist.size(); ++i)
        {
            auto &vHist = vvHist[i];
            int nvHist = vHist.size();
            vCan.push_back(new TCanvas(Form("c%d", i + 2), vHist[0].first->GetTitle(), 800, 600));
            TLegend *legend = (nvHist > 1) ? new TLegend(0.55, 0.7, 0.88, 0.88) : nullptr;

            double yMax = 0;
            for (auto &entry : vHist)
            {
                auto *h = entry.first;
                h->Scale(1.0 / eqTime); // Convert to cps
                yMax = max(yMax, entry.first->GetMaximum());
            }

            for (int j = 0; j < vHist.size(); ++j)
            {
                auto &h = vHist[j].first;
                auto &title = vHist[j].second;

                if (auto h2 = dynamic_cast<TH2 *>(h))
                {
                    gPad->SetRightMargin(0.2);
                    if (h2 == h2_cpposZ)
                    {
                        gPad->SetLogx();
                        gPad->SetLogy();
                        gPad->SetLogz();
                        h2->SetMinimum(0);
                    }
                    else if (std::find(vH2_cpposEY_byZ.begin(), vH2_cpposEY_byZ.end(), h2) != vH2_cpposEY_byZ.end())
                    {
                        gPad->SetLogx();
                        gPad->SetLogz();
                        h2->SetMinimum(1e-5);
                    }
                    else if (h2 == h2_cpposXY || h2 == h2_scposXY)
                    {
                        h2->SetMinimum(0);
                    }
                    h2->Draw("COLZ");
                }
                else if (auto h1 = dynamic_cast<TH1 *>(h))
                {
                    gPad->SetGridx();
                    gPad->SetGridy();
                    h1->SetMaximum(yMax * 1.3);
                    h1->SetMinimum(0);
                    h1->Draw(j == 0 ? "HIST E" : "HIST E SAME");
                }

                if (legend)
                    legend->AddEntry(h, title, "l");
            }
            if (legend)
                legend->Draw();
        }

        // save PDF
        if (true)
        {
            TString fPdfOut = "../fig/" + folder[folderID] + "_DetectPosition.pdf";
            if (vCan.size() == 1)
                vCan.at(0)->Print(fPdfOut);
            else
            {
                for (size_t i = 0; i < vCan.size(); ++i)
                {
                    if (i == 0)
                        vCan.at(i)->Print(fPdfOut + "(");
                    else if (i == vCan.size() - 1)
                        vCan.at(i)->Print(fPdfOut + ")");
                    else
                        vCan.at(i)->Print(fPdfOut);
                }
            }
        }

        for (auto *c : vCan)
            delete c;
        vCan.clear();
    }

    return;
}
