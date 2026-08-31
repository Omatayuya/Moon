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
    const Double_t DetectorOffsetZ = 400; // Detector offset in Z (mm)
    const Double_t irrArea = 600 * 600;   // irradiation surface area (cm^2) for Proton, Helium
    const Double_t SDArea = 10 * 10;      // Sensitive detector area (cm^2)

    // energy window
    constexpr double scatterEdepLow = 1.0; // MeV
    // constexpr double scatterEdepHigh = 3.0; // MeV
    constexpr double captureEdepLow = 4.5; // MeV

    // Thermal neutron cut (109Cd)
    constexpr double TNEnergyCut = 5e-7; // MeV

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
        int minXY = -75;    // mm
        int maxXY = 75;     // mm
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

    // === scatter/capture, scatter/capture_TNcut, capture/capture_TNcut vs 含水率 ===
    double TNLayerRange[] = {0.0, 20.0};
    double ENLayerRange[] = {40.0, 160.0};
    double TNlayerThickness = (TNLayerRange[1] - TNLayerRange[0]) * 0.1; // cm
    double ENlayerThickness = (ENLayerRange[1] - ENLayerRange[0]) * 0.1; // cm
    double scatterLayerThickness = TNlayerThickness + ENlayerThickness;  // cm
    double captureToRate = 1 / (TNlayerThickness * SDArea);              // cm^-1 cm^-2
    double captureTNcutToRate = 1 / (ENlayerThickness * SDArea);         // cm^-1 cm^-2
    double scatterToRate = 1 / (scatterLayerThickness * SDArea);         // cm^-1 cm^-2

    vector<double> vPpm;
    vector<double> vCp, vCpErr, vCpTNcut, vCpTNcutErr, vSc, vScErr;
    vector<double> vRatio_scCp, vRatio_scCpTNcut, vRatio_cpCpTNcut;
    vector<double> vRatio_scCpErr, vRatio_scCpTNcutErr, vRatio_cpCpTNcutErr;
    for (size_t i = 0; i < vHistcpPosZ.size(); ++i)
    {
        int bin_TN_low = vHistcpPosZ[i]->GetXaxis()->FindBin(TNLayerRange[0]);
        int bin_TN_high = vHistcpPosZ[i]->GetXaxis()->FindBin(TNLayerRange[1]) - 1;
        int bin_EN_low = vHistcpPosZ_TNcut[i]->GetXaxis()->FindBin(ENLayerRange[0]);
        int bin_EN_high = vHistcpPosZ_TNcut[i]->GetXaxis()->FindBin(ENLayerRange[1]) - 1;

        double scatterErr = 0;
        double scatterSum = vHistscPosZ[i]->IntegralAndError(1, vHistscPosZ[i]->GetNbinsX(), scatterErr);
        double captureErr = 0;
        double captureSum = vHistcpPosZ[i]->IntegralAndError(bin_TN_low, bin_TN_high, captureErr); // Z = 0-40 mm
        double captureTNcutErr = 0;
        double captureTNcutSum = vHistcpPosZ_TNcut[i]->IntegralAndError(bin_EN_low, bin_EN_high, captureTNcutErr); // Z = 40-120 mm

        double captureRate = captureSum * captureToRate;
        double captureErrRate = captureErr * captureToRate;
        double captureTNcutRate = captureTNcutSum * captureTNcutToRate;
        double captureTNcutErrRate = captureTNcutErr * captureTNcutToRate;
        double scatterRate = scatterSum * scatterToRate;
        double scatterErrRate = scatterErr * scatterToRate;

        vCp.push_back(captureRate);
        vCpErr.push_back(captureErrRate);
        vCpTNcut.push_back(captureTNcutRate);
        vCpTNcutErr.push_back(captureTNcutErrRate);
        vSc.push_back(scatterRate);
        vScErr.push_back(scatterErrRate);

        double ratio_scCp = captureRate / scatterRate;
        double ratio_scCpErr = ratio_scCp * sqrt(pow(scatterErrRate / scatterRate, 2) + pow(captureErrRate / captureRate, 2));
        double ratio_scCpTNcut = captureTNcutRate / scatterRate;
        double ratio_scCpTNcutErr = ratio_scCpTNcut * sqrt(pow(scatterErrRate / scatterRate, 2) + pow(captureTNcutErrRate / captureTNcutRate, 2));
        double ratio_cpCpTNcut = captureTNcutRate / captureRate;
        double ratio_cpCpTNcutErr = ratio_cpCpTNcut * sqrt(pow(captureErrRate / captureRate, 2) + pow(captureTNcutErrRate / captureTNcutRate, 2));

        vRatio_scCp.push_back(ratio_scCp);
        vRatio_scCpErr.push_back(ratio_scCpErr);
        vRatio_scCpTNcut.push_back(ratio_scCpTNcut);
        vRatio_scCpTNcutErr.push_back(ratio_scCpTNcutErr);
        vRatio_cpCpTNcut.push_back(ratio_cpCpTNcut);
        vRatio_cpCpTNcutErr.push_back(ratio_cpCpTNcutErr);

        // Extract ppm value from folder name
        TString folderName = folder[i];
        TString ppmStr = folderName(0, folderName.Length() - 3); // Remove "ppm" suffix
        double ppmValue = ppmStr.Atof();
        vPpm.push_back(ppmValue);
    }

    // --- count rate vs H content (capture / capture_TNcut / scatter) ---
    vector<vector<double> *> countRateY = {&vCp, &vCpTNcut, &vSc};
    vector<vector<double> *> countRateYErr = {&vCpErr, &vCpTNcutErr, &vScErr};
    vector<Color_t> countRateColor = {kRed + 1, kGreen + 2, kBlue + 1};
    vector<TString> countRateLabel = {"Capture", "Capture_TNcut", "Scatter"};

    TMultiGraph *mgCountRate = new TMultiGraph();
    mgCountRate->SetTitle("Count rate vs H content;H content (ppm);Count rate (s^{-1} cm^{-3})");
    mgCountRate->SetMinimum(0);
    TLegend *legCountRate = new TLegend(0.20, 0.20, 0.4, 0.38);

    size_t iZero = std::distance(vPpm.begin(), std::find(vPpm.begin(), vPpm.end(), 0.0));
    TMultiGraph *mgCountRate_0ppmRatio = new TMultiGraph();
    mgCountRate_0ppmRatio->SetTitle("Count rate (normalized to 0 ppm) vs H content;H content (ppm);Relative count rate ");
    TLegend *legCountRate_0ppmRatio = new TLegend(0.20, 0.20, 0.4, 0.38);

    for (size_t i = 0; i < countRateY.size(); ++i)
    {
        TGraphErrors *gr = new TGraphErrors(vPpm.size(), vPpm.data(),
                                            countRateY.at(i)->data(), 0, countRateYErr.at(i)->data());
        gr->SetMarkerColor(countRateColor.at(i));
        gr->SetLineColor(countRateColor.at(i));
        gr->SetMarkerStyle(20);
        mgCountRate->Add(gr, "PL");
        legCountRate->AddEntry(gr, countRateLabel.at(i), "lp");

        // --- 0 ppm で規格化 ---
        vector<double> &y = *countRateY.at(i);
        vector<double> &yErr = *countRateYErr.at(i);
        double y0 = y.at(iZero), y0Err = yErr.at(iZero);
        vector<double> yNorm(y.size()), yNormErr(y.size());
        for (size_t j = 0; j < y.size(); ++j)
        {
            yNorm.at(j) = y.at(j) / y0;
            yNormErr.at(j) = yNorm.at(j) * sqrt(pow(yErr.at(j) / y.at(j), 2) + pow(y0Err / y0, 2));
        }
        TGraphErrors *grNorm = new TGraphErrors(vPpm.size(), vPpm.data(), yNorm.data(), 0, yNormErr.data());
        grNorm->SetMarkerColor(countRateColor.at(i));
        grNorm->SetLineColor(countRateColor.at(i));
        grNorm->SetMarkerStyle(20);
        mgCountRate_0ppmRatio->Add(grNorm, "PL");
        legCountRate_0ppmRatio->AddEntry(grNorm, countRateLabel.at(i), "lp");
    }

    vector<TString> vCountRateRatioLabel = {"Capture / Scatter ", "Capture_TNcut / Scatter ", "Capture_TNcut  / Capture "};
    vector<Color_t> vCountRateRatioColor = {kPink - 1, kSpring - 1, kOrange + 1};
    vector<vector<double> *> vCountRateRatioY = {&vRatio_scCp, &vRatio_scCpTNcut, &vRatio_cpCpTNcut};
    vector<vector<double> *> vCountRateRatioYErr = {&vRatio_scCpErr, &vRatio_scCpTNcutErr, &vRatio_cpCpTNcutErr};

    TMultiGraph *mgCountRateRatio = new TMultiGraph();
    mgCountRateRatio->SetTitle("Count rate ratio vs H content;H content (ppm);Count rate ratio");
    mgCountRateRatio->SetMinimum(0);
    TLegend *legCountRateRatio = new TLegend(0.20, 0.50, 0.5, 0.68);

    TMultiGraph *mgCountRateRatio_0ppmRatio = new TMultiGraph();
    mgCountRateRatio_0ppmRatio->SetTitle("Count rate ratio (normalized to 0 ppm) vs H content;H content (ppm);Relative count rate ratio ");
    TLegend *legCountRateRatio_0ppmRatio = new TLegend(0.20, 0.50, 0.5, 0.68);

    for (size_t i = 0; i < vCountRateRatioY.size(); ++i)
    {
        TGraphErrors *gr = new TGraphErrors(vPpm.size(), vPpm.data(),
                                            vCountRateRatioY.at(i)->data(), 0, vCountRateRatioYErr.at(i)->data());
        gr->SetMarkerColor(vCountRateRatioColor.at(i));
        gr->SetLineColor(vCountRateRatioColor.at(i));
        gr->SetMarkerStyle(20);
        mgCountRateRatio->Add(gr, "PL");
        legCountRateRatio->AddEntry(gr, vCountRateRatioLabel.at(i), "lp");
        // --- 0 ppm で規格化 ---
        vector<double> &y = *vCountRateRatioY.at(i);
        vector<double> &yErr = *vCountRateRatioYErr.at(i);
        double y0 = y.at(iZero), y0Err = yErr.at(iZero);
        vector<double> yNorm(y.size()), yNormErr(y.size());
        for (size_t j = 0; j < y.size(); ++j)
        {
            yNorm.at(j) = y.at(j) / y0;
            yNormErr.at(j) = yNorm.at(j) * sqrt(pow(yErr.at(j) / y.at(j), 2) + pow(y0Err / y0, 2));
        }
        TGraphErrors *grNorm = new TGraphErrors(vPpm.size(), vPpm.data(), yNorm.data(), 0, yNormErr.data());
        grNorm->SetMarkerColor(vCountRateRatioColor.at(i));
        grNorm->SetLineColor(vCountRateRatioColor.at(i));
        grNorm->SetMarkerStyle(20);
        mgCountRateRatio_0ppmRatio->Add(grNorm, "PL");
        legCountRateRatio_0ppmRatio->AddEntry(grNorm, vCountRateRatioLabel.at(i), "lp");
    }

    // === 0 ppm と各 ppm を Nσ で分離するために必要な観測時間 vs 体積 ===
    const double volumeStart = 100.0; // cm^3
    const double volumeStep = 50.0;   // cm^3
    const double volumeEnd = 1000.0;  // cm^3
    vector<double> vVolume;
    for (double volume = volumeStart; volume <= volumeEnd; volume += volumeStep)
    {
        vVolume.push_back(volume);
    }
    const int nVolumePoints = vVolume.size();

    vector<Color_t> vSigTimePpmColor;
    for (size_t j = 0; j < vPpm.size(); ++j)
    {
        vSigTimePpmColor.push_back(gStyle->GetColorPalette(j * gStyle->GetNumberOfColors() / vPpm.size()));
    }

    vector<TString> vSigTimeLabel = {"Capture / Scatter", "Capture_TNcut / Scatter", "Capture_TNcut / Capture"};
    vector<vector<double> *> vSigTimeRatioY = {&vRatio_scCp, &vRatio_scCpTNcut, &vRatio_cpCpTNcut};
    vector<vector<double> *> vSigTimeRateNumer = {&vCp, &vCpTNcut, &vCpTNcut}; // 各ratioの分子レート
    vector<vector<double> *> vSigTimeRateDenom = {&vSc, &vSc, &vCp};           // 各ratioの分母レート

    vector<TMultiGraph *> vMgSigTime;
    vector<TLegend *> vLegSigTime;

    for (size_t i = 0; i < vSigTimeLabel.size(); ++i)
    {
        TMultiGraph *mgSigTime = new TMultiGraph();
        mgSigTime->SetTitle(Form("Observation time for %.0f#sigma separation vs 0 ppm (%s);Volume (cm^{3});Observation time (s)", nSigma, vSigTimeLabel.at(i).Data()));
        TLegend *legSigTime = new TLegend(0.55, 0.70, 0.85, 0.88);
        legSigTime->SetNColumns(2);

        for (size_t j = 0; j < vPpm.size(); ++j)
        {
            if (j == iZero)
                continue;

            double ratio0 = vSigTimeRatioY.at(i)->at(iZero);
            double rateNumer0 = vSigTimeRateNumer.at(i)->at(iZero);
            double rateDenom0 = vSigTimeRateDenom.at(i)->at(iZero);
            double ratioX = vSigTimeRatioY.at(i)->at(j);
            double rateNumerX = vSigTimeRateNumer.at(i)->at(j);
            double rateDenomX = vSigTimeRateDenom.at(i)->at(j);

            // captureRate/scatterRate等を真の値と仮定し、体積V・観測時間Tでのポアソン統計を
            // その都度作り直す（過去のMC統計量eqTimeは引き継がない）
            double delta = fabs(ratioX - ratio0);
            double C0 = pow(ratio0, 2) * (1.0 / rateDenom0 + 1.0 / rateNumer0); // ratioErr(V,T)^2 = C0 / (V*T)
            double CX = pow(ratioX, 2) * (1.0 / rateDenomX + 1.0 / rateNumerX);
            double time1cm3 = pow(nSigma, 2) * (C0 + CX) / pow(delta, 2); // 1 cm^3 での必要観測時間 (s)

            vector<double> vTime(nVolumePoints);
            for (int v = 0; v < nVolumePoints; ++v)
            {
                vTime.at(v) = time1cm3 / vVolume.at(v); // 観測時間は体積に反比例
            }

            TGraph *gr = new TGraph(nVolumePoints, vVolume.data(), vTime.data());
            gr->SetMarkerColor(vSigTimePpmColor.at(j));
            gr->SetLineColor(vSigTimePpmColor.at(j));
            gr->SetMarkerStyle(20 + j);
            mgSigTime->Add(gr, "PL");
            legSigTime->AddEntry(gr, folder.at(j), "lp");
        }

        vMgSigTime.push_back(mgSigTime);
        vLegSigTime.push_back(legSigTime);
    }

    // Draw
    /*Color Palette*/
    gStyle->SetPalette(kBird);
    int nColors = gStyle->GetNumberOfColors();
    vector<Color_t> colorPalette;
    for (size_t i = 0; i < vHistcpPosZ.size(); ++i)
    {
        colorPalette.push_back(gStyle->GetColorPalette(i * nColors / vHistcpPosZ.size()));
    }
    vector<TCanvas *> vCan;
    {
        // --- count rate vs H content ---
        TCanvas *cCountRate = new TCanvas("cCountRate", "Count rate vs H content", 800, 600);
        cCountRate->SetLogx();
        mgCountRate->SetMinimum(0);
        mgCountRate->Draw("A");
        legCountRate->Draw();
        vCan.push_back(cCountRate);

        // --- count rate (normalized to 0 ppm) vs H content ---
        TCanvas *cCountRate_0ppmRatio = new TCanvas("cCountRate_0ppmRatio", "Count rate (normalized to 0 ppm) vs H content", 800, 600);
        cCountRate_0ppmRatio->SetLogx();
        mgCountRate_0ppmRatio->SetMinimum(0);
        mgCountRate_0ppmRatio->Draw("A");
        legCountRate_0ppmRatio->Draw();
        vCan.push_back(cCountRate_0ppmRatio);

        // --- count rate ratio vs H content ---
        TCanvas *cCountRateRatio = new TCanvas("cCountRateRatio", "Count rate ratio vs H content", 800, 600);
        cCountRateRatio->SetLogx();
        mgCountRateRatio->SetMinimum(0);
        mgCountRateRatio->Draw("A");
        legCountRateRatio->Draw();
        vCan.push_back(cCountRateRatio);

        // --- count rate ratio (normalized to 0 ppm) vs H content ---
        TCanvas *cCountRateRatio_0ppmRatio = new TCanvas("cCountRateRatio_0ppmRatio", "Count rate ratio (normalized to 0 ppm) vs H content", 800, 600);
        cCountRateRatio_0ppmRatio->SetLogx();
        mgCountRateRatio_0ppmRatio->SetMinimum(0);
        mgCountRateRatio_0ppmRatio->Draw("A");
        legCountRateRatio_0ppmRatio->Draw();
        vCan.push_back(cCountRateRatio_0ppmRatio);

        // --- observation time vs volume (0 ppm vs each ppm, Nσ separation), 比率ごとに別キャンバス ---
        for (size_t i = 0; i < vMgSigTime.size(); ++i)
        {
            TCanvas *cSigTime = new TCanvas(Form("cSigTime_%zu", i), vSigTimeLabel.at(i), 800, 600);
            cSigTime->SetGridx(0);
            cSigTime->SetGridy(0);
            gPad->SetLogy();
            vMgSigTime.at(i)->Draw("A");
            vLegSigTime.at(i)->Draw();
            vCan.push_back(cSigTime);

            vector<double> vTime{60, 3600, 3600 * 24, 3600 * 24 * 7};
            vector<TString> vText{"1m", "1h", "1d", "1w"};
            for (size_t k = 0; k < vTime.size(); ++k)
            {
                TLine *l = new TLine(volumeStart, vTime[k], volumeEnd, vTime[k]);
                l->SetLineColor(kGray + 1);
                l->SetLineStyle(kDashed);
                l->Draw();

                TText *t = new TText(volumeStart + 10, vTime[k] * 0.5, vText[k]);
                t->Draw();
            }
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
