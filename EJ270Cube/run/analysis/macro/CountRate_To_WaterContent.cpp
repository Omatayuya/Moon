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
    vector<TString> folder = {"0ppm", "10ppm", "20ppm", "50ppm", "100ppm", "200ppm", "500ppm", "1000ppm", "2000ppm", "5000ppm", "10000ppm"};
    // vector<TString> folder = {"0ppm"};

    vector<TH1F *> vHistcpPosZ;
    vector<TH1F *> vHistcpPosZ_TNcut;
    vector<TH1F *> vHistscPosZ;

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

        // Primary energy bins
        constexpr int nEBins = 4;
        const double primEnergyEdges[nEBins + 1] = {0.0, 5e-7, 1e-3, 1.0, numeric_limits<double>::infinity()};
        const TString primEnergyLabels[nEBins] = {"E < 0.5 eV", "0.5 eV #leq E < 1 keV", "1 keV #leq E < 1 MeV", "E #geq 1 MeV"};
        const int primEnergyColors[nEBins] = {kOrange + 8, kGreen - 7, kGreen + 2, kBlue};

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
        // for (int e = 0; e < nEBins; ++e)
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

    // Draw
    /*Color Palette*/
    gStyle->SetPalette(kBird);
    int nColors = gStyle->GetNumberOfColors();
    vector<Color_t> colorPalette;
    for (int i = 0; i < vHistcpPosZ.size(); ++i)
    {
        colorPalette.push_back(gStyle->GetColorPalette(i * nColors / vHistcpPosZ.size()));
    }

    vector<TCanvas *> vCan;
    // === scatter/capture, scatter/capture_TNcut, capture/capture_TNcut vs 含水率 ===
    double TNLayerRange[] = {0.0, 40.0};
    double ENLayerRange[] = {40.0, 120.0};
    vector<double> vPpm;
    vector<double> vRatio_scCp, vRatio_scCpTNcut, vRatio_cpCpTNcut;
    vector<double> vRatio_scCpErr, vRatio_scCpTNcutErr, vRatio_cpCpTNcutErr;
    for (int i = 0; i < vHistcpPosZ.size(); ++i)
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

        double ratio_scCp = captureSum / scatterSum;
        double ratio_scCpErr = ratio_scCp * sqrt(pow(scatterErr / scatterSum, 2) + pow(captureErr / captureSum, 2));
        double ratio_scCpTNcut = captureTNcutSum / scatterSum;
        double ratio_scCpTNcutErr = ratio_scCpTNcut * sqrt(pow(scatterErr / scatterSum, 2) + pow(captureTNcutErr / captureTNcutSum, 2));
        double ratio_cpCpTNcut = captureTNcutSum / captureSum;
        double ratio_cpCpTNcutErr = ratio_cpCpTNcut * sqrt(pow(captureErr / captureSum, 2) + pow(captureTNcutErr / captureTNcutSum, 2));

        vRatio_scCp.push_back(captureSum / scatterSum);
        vRatio_scCpErr.push_back(ratio_scCpErr);
        vRatio_scCpTNcut.push_back(captureTNcutSum / scatterSum);
        vRatio_scCpTNcutErr.push_back(ratio_scCpTNcutErr);
        vRatio_cpCpTNcut.push_back(captureTNcutSum / captureSum);
        vRatio_cpCpTNcutErr.push_back(ratio_cpCpTNcutErr);

        // Extract ppm value from folder name
        TString folderName = folder[i];
        TString ppmStr = folderName(0, folderName.Length() - 3); // Remove "ppm" suffix
        double ppmValue = ppmStr.Atof();
        vPpm.push_back(ppmValue);
    }

    {
        // scatter/capture vs H content
        TCanvas *cRatio_scCp = new TCanvas("cRatio_scCp", "Capture/Scatter vs H content", 800, 600);
        cRatio_scCp->SetLogx();
        TGraphErrors *grRatio_scCp = new TGraphErrors(vPpm.size(), vPpm.data(), vRatio_scCp.data(), 0, vRatio_scCpErr.data());
        grRatio_scCp->SetTitle(Form("Capture / Scatter (Z %g-%g mm) vs H content;H content (ppm);Capture / Scatter", TNLayerRange[0], TNLayerRange[1]));
        grRatio_scCp->SetMarkerStyle(20);
        grRatio_scCp->SetMarkerColor(kBlue + 1);
        grRatio_scCp->SetLineColor(kBlue + 1);
        grRatio_scCp->Draw("APL");
        vCan.push_back(cRatio_scCp);

        // scatter/capture_TNcut vs H content
        TCanvas *cRatio_scCpTNcut = new TCanvas("cRatio_scCpTNcut", "Capture_TNcut/Scatter vs H content", 800, 600);
        cRatio_scCpTNcut->SetLogx();
        TGraphErrors *grRatio_scCpTNcut = new TGraphErrors(vPpm.size(), vPpm.data(), vRatio_scCpTNcut.data(), 0, vRatio_scCpTNcutErr.data());
        grRatio_scCpTNcut->SetTitle(Form("Capture_TNcut / Scatter (Z %g-%g mm) vs H content;H content (ppm);Capture_TNcut / Scatter", ENLayerRange[0], ENLayerRange[1]));
        grRatio_scCpTNcut->SetMarkerStyle(20);
        grRatio_scCpTNcut->SetMarkerColor(kBlue + 1);
        grRatio_scCpTNcut->SetLineColor(kBlue + 1);
        grRatio_scCpTNcut->Draw("APL");
        vCan.push_back(cRatio_scCpTNcut);

        // capture/capture_TNcut vs H content
        TCanvas *cRatio_cpCpTNcut = new TCanvas("cRatio_cpCpTNcut", "Capture_TNcut/Capture vs H content", 800, 600);
        cRatio_cpCpTNcut->SetLogx();
        TGraphErrors *grRatio_cpCpTNcut = new TGraphErrors(vPpm.size(), vPpm.data(), vRatio_cpCpTNcut.data(), 0, vRatio_cpCpTNcutErr.data());
        grRatio_cpCpTNcut->SetTitle(Form("Capture_TNcut (Z %g-%g mm) / Capture (Z %g-%g mm) vs H content;H content (ppm);Capture_TNcut / Capture", ENLayerRange[0], ENLayerRange[1], TNLayerRange[0], TNLayerRange[1]));
        grRatio_cpCpTNcut->SetMarkerStyle(20);
        grRatio_cpCpTNcut->SetMarkerColor(kBlue + 1);
        grRatio_cpCpTNcut->SetLineColor(kBlue + 1);
        grRatio_cpCpTNcut->Draw("APL");
        vCan.push_back(cRatio_cpCpTNcut);
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
