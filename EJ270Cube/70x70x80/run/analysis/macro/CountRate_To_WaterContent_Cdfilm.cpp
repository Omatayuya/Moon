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

void CountRate_To_WaterContent_Cdfilm()
{
    // ==================================================================
    // const Double_t irrArea = pow(450, 2);                    // irradiation surface area (cm^2) for Proton, Helium
    const Double_t DetectorOffsetZ = 460; // Detector offset in Z (mm)
    const Double_t irrArea = 600 * 600;   // irradiation surface area (cm^2) for Proton, Helium
    const Double_t SDArea = 10 * 10;      // Sensitive detector area (cm^2)

    // energy window
    constexpr double scatterEdepLow = 1.0; // MeV
    // constexpr double scatterEdepHigh = 3.0; // MeV
    constexpr double captureEdepLow = 4.5;  // MeV
    constexpr double captureEdepHigh = 5.0; // MeV

    // Thermal neutron cut (109Cd)
    constexpr double TNEnergyCut = 5e-7; // MeV

    // Cd plane position (mm)
    constexpr double sideCut = 5;               // mm
    const double fidHalfWidth = 35.0 - sideCut; // mm
    vector<double> cdPlanePosXY = {
        20 + DetectorOffsetZ,
    };
    vector<double> cdPlanePosZX = {-35, -35 + sideCut, 35 - sideCut, 35};
    vector<double> cdPlanePosZY = {-35, -35 + sideCut, 35 - sideCut, 35};

    // Primary energy bins
    constexpr int nEBins = 4;
    const double primEnergyEdges[nEBins + 1] = {0.0, 5e-7, 1e-3, 1.0, numeric_limits<double>::infinity()};
    const TString primEnergyLabels[nEBins] = {"E < 0.5 eV", "0.5 eV #leq E < 1 keV", "1 keV #leq E < 1 MeV", "E #geq 1 MeV"};
    const int primEnergyColors[nEBins] = {kOrange + 8, kGreen - 7, kGreen + 2, kBlue};

    // Sensitivity (0 ppm vs X ppm)
    const double nSigma = 5.0; // 要求する分離有意度 (σ)

    for (const TString &axis : {"X", "Y", "Z"})
    {
        gStyle->SetLabelFont(62, axis);
        gStyle->SetTitleFont(62, axis);
        if (axis == "Y")
            gStyle->SetTitleOffset(1.4, axis); // 軸タイトルのオフセット
        else
            gStyle->SetTitleOffset(1.2, axis); // 軸タイトルのオフセット
        gStyle->SetLabelSize(0.04, axis);      // 目盛り数字のサイズ
        gStyle->SetTitleSize(0.04, axis);      // 軸タイトルのサイズ
    }
    gStyle->SetTextFont(62);
    gStyle->SetTitleFont(62, "");

    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    // gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    gStyle->SetPadLeftMargin(0.15);

    // ==================================================================

    vector<TString> folder = {"0ppm", "10ppm", "20ppm", "50ppm", "100ppm", "200ppm", "500ppm", "1000ppm", "2000ppm", "5000ppm", "10000ppm"};
    // vector<TString> folder = {"0ppm"};

    vector<TH1F *> vHistcpPosZ;
    vector<TH1F *> vHistcpPosZ_TNcut;
    vector<TH1F *> vHistscPosZ;
    vector<double> vEqTime;

    for (int folderID = 0; folderID < folder.size(); folderID++)
    {

        TFile *fin = TFile::Open("../../" + folder[folderID] + "/results.root");
        if (fin)
        {
            cout << "Opened results.root in folder: " << folder[folderID] << endl;
        }
        else if (!fin || fin->IsZombie())
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
        vEqTime.push_back(eqTime);

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
        char fpname[256], fCProc[256];
        char collection[256], fPreProc[256], fPostProc[256];
        double fPreKinE, fPostKinE, fEdep, fGTime;
        double fPrePosX, fPrePosY, fPrePosZ, fPostPosX, fPostPosY, fPostPosZ;

        HitTree->SetBranchAddress("eventID", &eventID);
        HitTree->SetBranchAddress("primEnergy", &primEnergy);
        HitTree->SetBranchAddress("collection", collection);
        HitTree->SetBranchAddress("PID", &fPID);
        HitTree->SetBranchAddress("PPID", &fPPID);
        HitTree->SetBranchAddress("Pname", fpname);
        HitTree->SetBranchAddress("PreProc", fPreProc);
        HitTree->SetBranchAddress("PostProc", fPostProc);
        HitTree->SetBranchAddress("CProc", fCProc);
        HitTree->SetBranchAddress("PreKinE", &fPreKinE);
        HitTree->SetBranchAddress("PostKinE", &fPostKinE);
        HitTree->SetBranchAddress("Edep", &fEdep);
        HitTree->SetBranchAddress("GTime", &fGTime);
        HitTree->SetBranchAddress("PrePosX", &fPrePosX);
        HitTree->SetBranchAddress("PrePosY", &fPrePosY);
        HitTree->SetBranchAddress("PrePosZ", &fPrePosZ);
        HitTree->SetBranchAddress("PostPosX", &fPostPosX);
        HitTree->SetBranchAddress("PostPosY", &fPostPosY);
        HitTree->SetBranchAddress("PostPosZ", &fPostPosZ);
        
        const vector<string> captureStepParticles = {"alpha", "triton"};
        // const vector<string> scatterStepParticles = {"proton", "C12", "O16", "N14"};
        const vector<string> scatterStepParticles = {"proton"};

        set<int> CdCaptureEvents;
        int nHits = HitTree->GetEntries();
        for (int i = 0; i < nHits; ++i)
        {
            HitTree->GetEntry(i);

            if (CdCaptureEvents.count(eventID) > 0)
            {
                continue; // Skip if this event has already been processed for Cd capture
            }

            if (string(fpname) == "neutron")
            {
                bool cdAbsorbed = false;
                for (double zPlane : cdPlanePosXY) // XY平面（法線Z）
                    if ((fPrePosZ - zPlane) * (fPostPosZ - zPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                        cdAbsorbed = true;
                for (double yPlane : cdPlanePosZX) // ZX平面（法線Y）
                    if ((fPrePosY - yPlane) * (fPostPosY - yPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                        cdAbsorbed = true;
                for (double xPlane : cdPlanePosZY) // ZY平面（法線X）
                    if ((fPrePosX - xPlane) * (fPostPosX - xPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                        cdAbsorbed = true;
                if (cdAbsorbed)
                {
                    CdCaptureEvents.insert(eventID);
                    continue;
                }
            }

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

            if (string(collection) == "TrackerHitsCollection")
            {
                if (fabs(fPrePosX) > fidHalfWidth || fabs(fPrePosY) > fidHalfWidth)
                    continue;

                EventChamberID eventChamberID{eventID, primEnergy, string(collection)};
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

                    if (acc.captureflag == false && acc.cpEdepSum > captureEdepLow && acc.cpEdepSum < captureEdepHigh)
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
        cout << "Virtual Cd layer absorbed: " << CdCaptureEvents.size() << " / " << nEvents << " events" << endl;

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
        int maxZ = 80;     // mm
        int nBinsZ = (maxZ - minZ) / BinWidthZ;
        int BinWidthXY = 5; // mm
        int minXY = -35;    // mm
        int maxXY = 35;     // mm
        int nBinsXY = (maxXY - minXY) / BinWidthXY;

        TH1F *h1_cpposZ = new TH1F("h1_cpposZ", Form("Capture Position Z Distribution;Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        TH1F *h1_cpposZ_TNcut = new TH1F("h1_cpposZ_TNcut", Form("Capture Position Z Distribution (TN cut / %.1f eV <);Z (mm);Counts (s^{-1} %d mm^{-1})", TNEnergyCut * 1e6, BinWidthZ), nBinsZ, minZ, maxZ);
        TH2D *h2_cpposXY = new TH2D("h2_cpposXY", Form("Capture Position XY Distribution;X (mm);Y (mm);Counts (s^{-1} %d #times %d mm^{-2})", BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        TH2D *h2_cpposXY_TNcut = new TH2D("h2_cpposXY_TNcut", Form("Capture Position XY Distribution (TN cut / %.1f eV <);X (mm);Y (mm);Counts (s^{-1} %d #times %d mm^{-2})", TNEnergyCut * 1e6, BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        TH1F *h1_scposZ = new TH1F("h1_scposZ", Form("Scatter Position Z Distribution;Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        TH2D *h2_scposXY = new TH2D("h2_scposXY", Form("Scatter Position XY Distribution;X (mm);Y (mm);Counts (s^{-1} %d #times %d mm^{-2})", BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);

        // vector<TH1F *> vH_cpposZ_byE(nEBins);
        // for (size_t e = 0; e < nEBins; ++e)
        // {
        //     TString hname = Form("h1_cpposZ_E%d", e);
        //     vH_cpposZ_byE[e] = new TH1F(hname,
        //                                 Form("Capture Position Z (%s);Z (mm);Counts (s^{-1} %d mm^{-1})", primEnergyLabels[e].Data(), BinWidthZ),
        //                                 nBinsZ, minZ, maxZ);
        // }

        for (int i = 0; i < entries; ++i)
        {
            DetectorchamberTree->GetEntry(i);
            double capturePosZ = chamberEventData.cpPosZ - DetectorOffsetZ; // Convert m to mm and subtract detector offset
            double capturePosX = chamberEventData.cpPosX;                   // Convert m to mm
            double capturePosY = chamberEventData.cpPosY;                   // Convert m to mm
            if (chamberEventData.captureflag)
            {
                h1_cpposZ->Fill(capturePosZ);
                h2_cpposXY->Fill(capturePosX, capturePosY);
                if (eventChamberID.primEnergy > TNEnergyCut) // TN cut
                {
                    h1_cpposZ_TNcut->Fill(capturePosZ);
                    h2_cpposXY_TNcut->Fill(capturePosX, capturePosY);
                }

                // // primEnergyで分岐してFill
                // for (int e = 0; e < nEBins; ++e)
                // {
                //     if (eventChamberID.primEnergy >= primEnergyEdges[e] && eventChamberID.primEnergy < primEnergyEdges[e + 1])
                //     {
                //         vH_cpposZ_byE[e]->Fill(capturePosZ);
                //         break;
                //     }
                // }
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

        h1_cpposZ->Scale(1.0 / eqTime);
        h1_cpposZ_TNcut->Scale(1.0 / eqTime);
        h2_cpposXY->Scale(1.0 / eqTime);
        h2_cpposXY_TNcut->Scale(1.0 / eqTime);
        h1_scposZ->Scale(1.0 / eqTime);
        h2_scposXY->Scale(1.0 / eqTime);

        vHistcpPosZ.push_back(h1_cpposZ);
        vHistcpPosZ_TNcut.push_back(h1_cpposZ_TNcut);
        vHistscPosZ.push_back(h1_scposZ);
    }

    // === 各層 capture_TNcut / 全 scatter vs 含水率 ===
    vector<double> layerEdges = {0, 5, 10, 20, 40, 60, 70, 75, 80}; // mm（前面から）。DetectPosition の zRegionEdges と同じ流儀
    const int nLayers = layerEdges.size() - 1;
    vector<TString> layerLabels;
    for (int k = 0; k < nLayers; ++k)
        layerLabels.push_back(Form("%.0f #leq Z < %.0f mm", layerEdges[k], layerEdges[k + 1]));

    vector<double> vPpm;
    vector<vector<double>> vLayerRatio(nLayers), vLayerRatioErr(nLayers); // [layer][ppmIndex]
    for (size_t i = 0; i < vHistcpPosZ.size(); ++i)
    {
        double scatterErr = 0;
        double scatterSum = vHistscPosZ[i]->IntegralAndError(1, vHistscPosZ[i]->GetNbinsX(), scatterErr);

        for (int k = 0; k < nLayers; ++k)
        {
            int binLo = vHistcpPosZ[i]->GetXaxis()->FindBin(layerEdges[k]);
            int binHi = vHistcpPosZ[i]->GetXaxis()->FindBin(layerEdges[k + 1]) - 1;
            double capErr = 0;
            double capSum = vHistcpPosZ[i]->IntegralAndError(binLo, binHi, capErr);

            double ratio = (capSum > 0 && scatterSum > 0) ? capSum / scatterSum : 0.0;
            double ratioErr = (ratio > 0)
                                  ? ratio * sqrt(pow(capErr / capSum, 2) + pow(scatterErr / scatterSum, 2))
                                  : 0.0;
            vLayerRatio[k].push_back(ratio);
            vLayerRatioErr[k].push_back(ratioErr);
        }

        TString folderName = folder[i];
        TString ppmStr = folderName(0, folderName.Length() - 3);
        vPpm.push_back(ppmStr.Atof());
    }

    // --- 各層 比率の 0ppm 規格化 ---
    size_t iZero = std::distance(vPpm.begin(), std::find(vPpm.begin(), vPpm.end(), 0.0));

    gStyle->SetPalette(kBird);
    int nColors = gStyle->GetNumberOfColors();

    TMultiGraph *mgLayerRatio_0ppm = new TMultiGraph();
    mgLayerRatio_0ppm->SetTitle("Count rate ratio (capture_{TNcut} / scatter, normalized to 0 ppm) vs H content;H content (ppm);Relative count rate ratio");
    TLegend *legLayerRatio_0ppm = new TLegend(0.18, 0.60, 0.45, 0.88);

    for (int k = 0; k < nLayers; ++k)
    {
        vector<double> &y = vLayerRatio[k];
        vector<double> &yErr = vLayerRatioErr[k];
        double y0 = y.at(iZero), y0Err = yErr.at(iZero);
        vector<double> yNorm(y.size()), yNormErr(y.size());
        for (size_t j = 0; j < y.size(); ++j)
        {
            yNorm.at(j) = (y0 > 0) ? y.at(j) / y0 : 0.0;
            yNormErr.at(j) = (yNorm.at(j) > 0 && y.at(j) > 0)
                                 ? yNorm.at(j) * sqrt(pow(yErr.at(j) / y.at(j), 2) + pow(y0Err / y0, 2))
                                 : 0.0;
        }
        TGraphErrors *gr = new TGraphErrors(vPpm.size(), vPpm.data(), yNorm.data(), 0, yNormErr.data());
        Color_t col = gStyle->GetColorPalette(nLayers > 1 ? k * (nColors - 1) / (nLayers - 1) : 0);
        gr->SetMarkerColor(col);
        gr->SetLineColor(col);
        gr->SetMarkerStyle(20 + (k % 10));
        mgLayerRatio_0ppm->Add(gr, "PL");
        legLayerRatio_0ppm->AddEntry(gr, layerLabels[k], "lp");
    }

    constexpr Double_t xmin = 5;    // H (ppm)
    constexpr Double_t xmax = 2e+4; // H (ppm)

    // Draw
    vector<TCanvas *> vCan;
    {
        TCanvas *cLayerRatio_0ppm = new TCanvas("cLayerRatio_0ppm", "Layer count rate ratio (normalized to 0 ppm)", 800, 600);
        cLayerRatio_0ppm->SetLogx();
        mgLayerRatio_0ppm->Draw("A");
        legLayerRatio_0ppm->Draw();
        vCan.push_back(cLayerRatio_0ppm);
    }

    // save PDF
    if (true)
    {
        TString fPdfOut = "../fig/CountRate_To_WaterContent.pdf";
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

    return;
}
