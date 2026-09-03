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

void CountRate_To_WaterContent()
{
    // ==================================================================
    // const Double_t irrArea = pow(450, 2);                    // irradiation surface area (cm^2) for Proton, Helium
    const Double_t DetectorOffsetZ = 460; // Detector offset in Z (mm)
    const Double_t irrArea = 600 * 600;   // irradiation surface area (cm^2) for Proton, Helium

    // energy window
    constexpr double scatterEdepLow = 1.0;  // MeV
    constexpr double scatterEdepHigh = 3.0; // MeV
    constexpr double captureEdepLow = 4.5;  // MeV
    constexpr double captureEdepHigh = 5.0; // MeV

    // Thermal neutron cut (109Cd)
    constexpr double TNEnergyCut = 5e-7; // MeV

    // Fiducial cut & virtual Cd layer (DetectPosition.cpp と同じ)
    constexpr Double_t EJ270HalfWidth = 35;                                                           // mm
    constexpr double sideCut = 5;                                                                     // mm
    const double fidHalfWidth = EJ270HalfWidth - sideCut;                                             // mm
    vector<double> cdPlanePosXY = {10 + DetectorOffsetZ, 20 + DetectorOffsetZ, 70 + DetectorOffsetZ}; // XY平面 (法線: Z軸) [mm]
    vector<double> cdPlanePosZX = {-EJ270HalfWidth + sideCut, EJ270HalfWidth - sideCut};              // ZX平面 (法線: Y軸) [mm]
    vector<double> cdPlanePosZY = {-EJ270HalfWidth + sideCut, EJ270HalfWidth - sideCut};              // ZY平面 (法線: X軸) [mm]

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

        set<int> CdCutEventSets; // Set to store event IDs that pass the Cd cut
        int nHits = HitTree->GetEntries();
        for (int i = 0; i < nHits; ++i)
        {
            HitTree->GetEntry(i);

            if (CdCutEventSets.count(eventID) > 0)
            {
                continue;
            }

            if (string(fpname) == "neutron")
            {
                bool cdAbsorbed = false;
                for (double zPlane : cdPlanePosXY) // XY平面 (法線: Z軸)
                {
                    if ((fPrePosZ - zPlane) * (fPostPosZ - zPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                        cdAbsorbed = true;
                }
                for (double yPlane : cdPlanePosZX) // ZX平面 (法線: Y軸)
                {
                    if ((fPrePosY - yPlane) * (fPostPosY - yPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                        cdAbsorbed = true;
                }
                for (double xPlane : cdPlanePosZY) // ZY平面 (法線: X軸)
                {
                    if ((fPrePosX - xPlane) * (fPostPosX - xPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                        cdAbsorbed = true;
                }
                if (cdAbsorbed)
                {
                    CdCutEventSets.insert(eventID);
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
                    if (acc.scatterflag == false && acc.scEdepSum > scatterEdepLow && acc.scEdepSum < scatterEdepHigh)
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

        cout << "Virtual Cd layer absorbed: " << CdCutEventSets.size() << " / " << nEvents << " events" << endl;

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

    // === 各検出層ごとの capture レートと、全 scatter に対する比率 vs 含水率 ===
    vector<double> vSensThick = {0, 5, 10, 20, 40}; // mm、DetectPosition.cpp と同じ流儀
    const int nSensThick = vSensThick.size() - 1;
    vector<TString> zRegionLabels;
    for (int z = 0; z < nSensThick; ++z)
        zRegionLabels.push_back(Form("%.0f #leq Z < %.0f mm", vSensThick[z], vSensThick[z + 1]));

    vector<double> vPpm;
    vector<double> vSc, vScErr;                                                 // 全 scatter レート (s^-1)
    vector<vector<double>> vCpLayer(nSensThick), vCpLayerErr(nSensThick);       // [layer][ppm] capture レート (s^-1)
    vector<vector<double>> vRatioLayer(nSensThick), vRatioLayerErr(nSensThick); // [layer][ppm] capture_layer / scatter_total
    for (size_t i = 0; i < vHistcpPosZ.size(); ++i)
    {
        double scatterErr = 0;
        double scatterSum = vHistscPosZ[i]->IntegralAndError(1, vHistscPosZ[i]->GetNbinsX(), scatterErr);
        vSc.push_back(scatterSum);
        vScErr.push_back(scatterErr);

        for (int z = 0; z < nSensThick; ++z)
        {
            int binLo = vHistcpPosZ[i]->GetXaxis()->FindBin(vSensThick[z]);
            int binHi = vHistcpPosZ[i]->GetXaxis()->FindBin(vSensThick[z + 1]) - 1;
            double capErr = 0;
            double capSum = vHistcpPosZ[i]->IntegralAndError(binLo, binHi, capErr);
            vCpLayer[z].push_back(capSum);
            vCpLayerErr[z].push_back(capErr);

            double ratio = (capSum > 0 && scatterSum > 0) ? capSum / scatterSum : 0.0;
            double ratioErr = (ratio > 0)
                                  ? ratio * sqrt(pow(capErr / capSum, 2) + pow(scatterErr / scatterSum, 2))
                                  : 0.0;
            vRatioLayer[z].push_back(ratio);
            vRatioLayerErr[z].push_back(ratioErr);
        }

        TString folderName = folder[i];
        TString ppmStr = folderName(0, folderName.Length() - 3);
        vPpm.push_back(ppmStr.Atof());
    }

    size_t iZero = std::distance(vPpm.begin(), std::find(vPpm.begin(), vPpm.end(), 0.0));

    gStyle->SetPalette(kRainBow);
    int nColors = gStyle->GetNumberOfColors();

    TMultiGraph *mgCpLayer = new TMultiGraph();
    mgCpLayer->SetTitle("Capture count rate per layer vs H content;H content (ppm);Count rate (s^{-1})");
    TLegend *legCpLayer = new TLegend(0.18, 0.12, 0.45, 0.25);

    TMultiGraph *mgRatioLayer = new TMultiGraph();
    mgRatioLayer->SetTitle("Count rate ratio (capture layer / scatter total) vs H content;H content (ppm);Count rate ratio");
    mgRatioLayer->SetMinimum(0);
    TLegend *legRatioLayer = new TLegend(0.18, 0.12, 0.45, 0.25);

    TMultiGraph *mgRatioLayer_0ppm = new TMultiGraph();
    mgRatioLayer_0ppm->SetTitle("Count rate ratio (normalized to 0 ppm) vs H content;H content (ppm);Relative count rate ratio");
    TLegend *legRatioLayer_0ppm = new TLegend(0.18, 0.75, 0.45, 0.88);

    TMultiGraph *mgSigTimeLayer = new TMultiGraph();
    mgSigTimeLayer->SetTitle(Form("Observation time for %.0f#sigma separation from 0 ppm;H content (ppm);Observation time (s)", nSigma));
    TLegend *legSigTimeLayer = new TLegend(0.55, 0.75, 0.85, 0.88);

    for (int z = 0; z < nSensThick; ++z)
    {
        Color_t col = gStyle->GetColorPalette(nSensThick > 1 ? static_cast<int>((0.1 + 0.8 * z / (nSensThick - 1)) * (nColors - 1)) : 0);
        int mstyle = 20 + (z % 10);

        TGraphErrors *grCp = new TGraphErrors(vPpm.size(), vPpm.data(), vCpLayer[z].data(), 0, vCpLayerErr[z].data());
        grCp->SetMarkerColor(col);
        grCp->SetLineColor(col);
        grCp->SetMarkerStyle(mstyle);
        mgCpLayer->Add(grCp, "PL");
        legCpLayer->AddEntry(grCp, zRegionLabels[z], "lp");

        TGraphErrors *grRatio = new TGraphErrors(vPpm.size(), vPpm.data(), vRatioLayer[z].data(), 0, vRatioLayerErr[z].data());
        grRatio->SetMarkerColor(col);
        grRatio->SetLineColor(col);
        grRatio->SetMarkerStyle(mstyle);
        mgRatioLayer->Add(grRatio, "PL");
        legRatioLayer->AddEntry(grRatio, zRegionLabels[z], "lp");

        vector<double> &y = vRatioLayer[z];
        vector<double> &yErr = vRatioLayerErr[z];
        double y0 = y.at(iZero), y0Err = yErr.at(iZero);
        vector<double> yNorm(y.size()), yNormErr(y.size());
        for (size_t j = 0; j < y.size(); ++j)
        {
            yNorm.at(j) = (y0 > 0) ? y.at(j) / y0 : 0.0;
            yNormErr.at(j) = (yNorm.at(j) > 0 && y.at(j) > 0)
                                 ? yNorm.at(j) * sqrt(pow(yErr.at(j) / y.at(j), 2) + pow(y0Err / y0, 2))
                                 : 0.0;
        }
        TGraphErrors *grNorm = new TGraphErrors(vPpm.size(), vPpm.data(), yNorm.data(), 0, yNormErr.data());
        grNorm->SetMarkerColor(col);
        grNorm->SetLineColor(col);
        grNorm->SetMarkerStyle(mstyle);
        mgRatioLayer_0ppm->Add(grNorm, "PL");
        legRatioLayer_0ppm->AddEntry(grNorm, zRegionLabels[z], "lp");

        double ratio0 = vRatioLayer[z].at(iZero);
        double rateCap0 = vCpLayer[z].at(iZero);
        double rateSc0 = vSc.at(iZero);
        if (ratio0 <= 0 || rateCap0 <= 0 || rateSc0 <= 0)
            continue;

        vector<double> vTimePpm, vTimeReq;
        for (size_t j = 0; j < vPpm.size(); ++j)
        {
            if (j == iZero)
                continue;

            double ratioX = vRatioLayer[z].at(j);
            double rateCapX = vCpLayer[z].at(j);
            double rateScX = vSc.at(j);
            double delta = fabs(ratioX - ratio0);
            if (delta <= 0 || ratioX <= 0 || rateCapX <= 0)
                continue;

            // ratioErr(T)^2 = ratio^2 (1/rateCap + 1/rateSc) / T   (ポアソン統計、T は観測時間 s)
            double C0 = pow(ratio0, 2) * (1.0 / rateSc0 + 1.0 / rateCap0);
            double CX = pow(ratioX, 2) * (1.0 / rateScX + 1.0 / rateCapX);
            double reqTime = pow(nSigma, 2) * (C0 + CX) / pow(delta, 2); // s

            vTimePpm.push_back(vPpm.at(j));
            vTimeReq.push_back(reqTime);
        }
        if (vTimePpm.empty())
            continue;

        TGraph *grTime = new TGraph(vTimePpm.size(), vTimePpm.data(), vTimeReq.data());
        grTime->SetMarkerColor(col);
        grTime->SetLineColor(col);
        grTime->SetMarkerStyle(mstyle);
        mgSigTimeLayer->Add(grTime, "PL");
        legSigTimeLayer->AddEntry(grTime, zRegionLabels[z], "lp");
    }

    TGraphErrors *grSc = new TGraphErrors(vPpm.size(), vPpm.data(), vSc.data(), 0, vScErr.data());
    grSc->SetMarkerColor(kBlack);
    grSc->SetLineColor(kBlack);
    grSc->SetMarkerStyle(24);
    mgCpLayer->Add(grSc, "PL");
    legCpLayer->AddEntry(grSc, "Scatter (total)", "lp");

    constexpr Double_t xmin = 5;    // H (ppm)
    constexpr Double_t xmax = 2e+4; // H (ppm)

    // Draw
    vector<TCanvas *> vCan;
    {
        TCanvas *cCpLayer = new TCanvas("cCpLayer", "Capture count rate per layer vs H content", 800, 600);
        cCpLayer->SetLogx();
        cCpLayer->SetLogy();
        mgCpLayer->Draw("A");
        mgCpLayer->GetXaxis()->SetLimits(xmin, xmax);
        legCpLayer->Draw();
        vCan.push_back(cCpLayer);

        TCanvas *cRatioLayer = new TCanvas("cRatioLayer", "Count rate ratio per layer vs H content", 800, 600);
        cRatioLayer->SetLogx();
        mgRatioLayer->SetMinimum(0);
        mgRatioLayer->Draw("A");
        mgRatioLayer->GetXaxis()->SetLimits(xmin, xmax);
        legRatioLayer->Draw();
        vCan.push_back(cRatioLayer);

        TCanvas *cRatioLayer_0ppm = new TCanvas("cRatioLayer_0ppm", "Count rate ratio per layer (normalized to 0 ppm)", 800, 600);
        cRatioLayer_0ppm->SetLogx();
        mgRatioLayer_0ppm->Draw("A");
        mgRatioLayer_0ppm->GetXaxis()->SetLimits(xmin, xmax);
        legRatioLayer_0ppm->Draw();
        vCan.push_back(cRatioLayer_0ppm);

        TCanvas *cSigTimeLayer = new TCanvas("cSigTimeLayer", "Observation time for Nsigma separation vs H content", 800, 600);
        cSigTimeLayer->SetLogx();
        cSigTimeLayer->SetLogy();
        cSigTimeLayer->SetGridx(0);
        cSigTimeLayer->SetGridy(0);
        mgSigTimeLayer->Draw("A");
        mgSigTimeLayer->GetXaxis()->SetLimits(xmin, xmax);
        legSigTimeLayer->Draw();
        vCan.push_back(cSigTimeLayer);

        vector<double> vTimeLine{60, 3600, 3600 * 24, 3600 * 24 * 7};
        vector<TString> vTimeText{"1m", "1h", "1d", "1w"};
        for (size_t k = 0; k < vTimeLine.size(); ++k)
        {
            TLine *l = new TLine(xmin, vTimeLine[k], xmax, vTimeLine[k]);
            l->SetLineColor(kGray + 1);
            l->SetLineStyle(kDashed);
            l->Draw();
            TText *t = new TText(xmin * 1.2, vTimeLine[k] * 0.5, vTimeText[k]);
            t->Draw();
        }
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
