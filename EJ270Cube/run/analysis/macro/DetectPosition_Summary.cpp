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

void DetectPosition_Summary()
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
    TCanvas *ccpPosZ = new TCanvas("cCpposZ", "Capture Position Z by H ppm", 800, 600);
    TLegend *legcpPosZ = new TLegend(0.55, 0.7, 0.88, 0.88);
    legcpPosZ->SetNColumns(2);
    for (int i = 0; i < vHistcpPosZ.size(); ++i)
    {
        vHistcpPosZ[i]->SetLineColor(colorPalette[i]);
        vHistcpPosZ[i]->SetLineWidth(2);
        vHistcpPosZ[i]->SetMaximum(25);
        vHistcpPosZ[i]->SetMinimum(0);
        vHistcpPosZ[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpPosZ->AddEntry(vHistcpPosZ[i], Form("%s", folder[i].Data()), "l");
    }
    legcpPosZ->Draw();
    vCan.push_back(ccpPosZ);

    // === フィットあり版 (別キャンバス。クローンを使い元のヒストグラムは変更しない) ===
    vector<TH1F *> vHistcpPosZ_fit;
    vector<double> vExpParam, vExpParamErr, vPpm;
    vector<TF1 *> vfExpFits;
    TCanvas *ccpPosZ_fit = new TCanvas("cCpposZ_fit", "Capture Position Z by H ppm (with fit)", 800, 600);
    TLegend *legcpPosZ_fit = new TLegend(0.55, 0.7, 0.88, 0.88);
    legcpPosZ_fit->SetNColumns(2);
    for (int i = 0; i < vHistcpPosZ.size(); ++i)
    {
        TH1F *h_fit = (TH1F *)vHistcpPosZ[i]->Clone(Form("h1_cpposZ_fit_%d", i));
        h_fit->SetLineColor(colorPalette[i]);
        h_fit->SetLineWidth(2);
        // h_fit->SetMaximum(25);
        // h_fit->SetMinimum(0);
        h_fit->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpPosZ_fit->AddEntry(h_fit, Form("%s", folder[i].Data()), "l");

        // --- exp + const フィット (Z < 40 mm) ---
        TF1 *fExp = new TF1(Form("fExp_%d", i), "[0]*exp([1]*x)+[2]", 0, 40);
        fExp->SetParameters(h_fit->GetMaximum(), -0.05, 1.0);
        fExp->SetLineColor(kRed);
        fExp->SetLineWidth(2);
        fExp->SetNpx(1000);
        h_fit->Fit(fExp, "RQ+");
        fExp->Draw(i == 0 ? "SAME" : "SAME");

        vfExpFits.push_back(fExp);
        double p1 = fExp->GetParameter(1);
        double p1_err = fExp->GetParError(1);
        vExpParam.push_back(p1);
        vExpParamErr.push_back(p1_err);

        TString ppmStr = folder[i];
        ppmStr.ReplaceAll("ppm", "");
        vPpm.push_back(ppmStr.Atof());

        vHistcpPosZ_fit.push_back(h_fit);
    }
    legcpPosZ_fit->Draw();
    vCan.push_back(ccpPosZ_fit);

    // === 指数パラメータ vs 含水率 ===
    TCanvas *cExpParam = new TCanvas("cExpParam", "Exponential Parameter vs H content", 800, 600);
    cExpParam->SetLogx();
    TGraphErrors *grExpParam = new TGraphErrors(vPpm.size(), vPpm.data(), vExpParam.data(), nullptr, vExpParamErr.data());
    grExpParam->SetTitle("Exponential Parameter vs H content;H content (ppm);Exp decay parameter");
    grExpParam->SetMarkerStyle(20);
    grExpParam->SetMarkerColor(kBlue + 1);
    grExpParam->SetLineColor(kBlue + 1);
    grExpParam->Draw("APL");
    vCan.push_back(cExpParam);

    TCanvas *ccpPosZ_ratio = new TCanvas("cCpposZ_ratio", "Capture Position Z by H ppm (Ratio)", 800, 600);
    TLegend *legcpPosZ_ratio = new TLegend(0.55, 0.12, 0.88, 0.3);
    legcpPosZ_ratio->SetNColumns(2);
    TH1F *h1_cpposZ_0ppm = (TH1F *)vHistcpPosZ[0]->Clone("h1_cpposZ_ratio");
    vector<TH1F *> vHistcpPosZ_ratio;
    for (int i = 0; i < vHistcpPosZ.size(); ++i)
    {
        TH1F *h1_cpposZ_ratio = (TH1F *)vHistcpPosZ[i]->Clone(Form("h1_cpposZ_ratio_%d", i));
        h1_cpposZ_ratio->Divide(h1_cpposZ_0ppm);
        vHistcpPosZ_ratio.push_back(h1_cpposZ_ratio);
        vHistcpPosZ_ratio[i]->SetTitle("Capture Position Z Distribution Ratio;Z (mm);Ratio to 0 ppm");
        vHistcpPosZ_ratio[i]->SetLineColor(colorPalette[i]);
        vHistcpPosZ_ratio[i]->SetLineWidth(2);
        vHistcpPosZ_ratio[i]->SetMaximum(1.2);
        vHistcpPosZ_ratio[i]->SetMinimum(0);
        vHistcpPosZ_ratio[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpPosZ_ratio->AddEntry(vHistcpPosZ_ratio[i], Form("%s", folder[i].Data()), "l");
    }
    legcpPosZ_ratio->Draw();
    vCan.push_back(ccpPosZ_ratio);

    TCanvas *ccpscRatioZ = new TCanvas("cCpScRatioZ", "Capture/Scatter Ratio by Z (per H ppm)", 800, 600);

    TLegend *legcpscRatioZ = new TLegend(0.55, 0.7, 0.88, 0.88);
    legcpscRatioZ->SetNColumns(2);
    vector<TH1F *> vHistcpscRatioZ;
    for (int i = 0; i < vHistcpPosZ.size(); ++i)
    {
        TH1F *h1_cpScRatioZ = (TH1F *)vHistcpPosZ[i]->Clone(Form("h1_cpScRatioZ_%d", i));
        h1_cpScRatioZ->Divide(vHistscPosZ[i]);
        h1_cpScRatioZ->SetTitle("Capture/Scatter Ratio;Z (mm);Capture / Scatter");
        h1_cpScRatioZ->SetLineColor(colorPalette[i]);
        h1_cpScRatioZ->SetLineWidth(2);
        h1_cpScRatioZ->SetMaximum(4);
        h1_cpScRatioZ->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpscRatioZ->AddEntry(h1_cpScRatioZ, Form("%s", folder[i].Data()), "l");
        vHistcpscRatioZ.push_back(h1_cpScRatioZ);
    }
    legcpscRatioZ->Draw();
    vCan.push_back(ccpscRatioZ);

    // === フィットあり版 (Capture/Scatter Ratio, 別キャンバス) ===
    vector<TH1F *> vHistcpscRatioZ_fit;
    vector<double> vExpParam_cpscRatio, vExpParamErr_cpscRatio, vPpm_cpscRatio;
    vector<TF1 *> vfExpFits_cpscRatio;
    TCanvas *ccpscRatioZ_fit = new TCanvas("cCpScRatioZ_fit", "Capture/Scatter Ratio by Z (per H ppm, with fit)", 800, 600);
    TLegend *legcpscRatioZ_fit = new TLegend(0.55, 0.7, 0.88, 0.88);
    legcpscRatioZ_fit->SetNColumns(2);
    for (int i = 0; i < vHistcpscRatioZ.size(); ++i)
    {
        TH1F *h_fit = (TH1F *)vHistcpscRatioZ[i]->Clone(Form("h1_cpScRatioZ_fit_%d", i));
        h_fit->SetLineColor(colorPalette[i]);
        h_fit->SetLineWidth(2);
        h_fit->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpscRatioZ_fit->AddEntry(h_fit, Form("%s", folder[i].Data()), "l");

        TF1 *fExp = new TF1(Form("fExp_cpscRatio_%d", i), "[0]*exp([1]*x)+[2]", 0, 40);
        fExp->SetParameters(h_fit->GetMaximum(), -0.05, 1.0);
        fExp->SetLineColor(kRed);
        fExp->SetLineWidth(2);
        fExp->SetNpx(1000);
        h_fit->Fit(fExp, "RQ+");
        fExp->Draw("same");

        vfExpFits_cpscRatio.push_back(fExp);
        double p1 = fExp->GetParameter(1);
        double p1_err = fExp->GetParError(1);
        vExpParam_cpscRatio.push_back(p1);
        vExpParamErr_cpscRatio.push_back(p1_err);

        TString ppmStr = folder[i];
        ppmStr.ReplaceAll("ppm", "");
        vPpm_cpscRatio.push_back(ppmStr.Atof());

        vHistcpscRatioZ_fit.push_back(h_fit);
    }
    legcpscRatioZ_fit->Draw();
    vCan.push_back(ccpscRatioZ_fit);

    // === 指数パラメータ vs 含水率 (Capture/Scatter Ratio) ===
    TCanvas *cExpParam_cpscRatio = new TCanvas("cExpParam_cpscRatio", "Exponential Parameter vs H content (Capture/Scatter Ratio)", 800, 600);
    cExpParam_cpscRatio->SetLogx();
    // cExpParam_cpscRatio->SetLogy();
    TGraphErrors *grExpParam_cpscRatio = new TGraphErrors(vPpm_cpscRatio.size(), vPpm_cpscRatio.data(), vExpParam_cpscRatio.data(), nullptr, vExpParamErr_cpscRatio.data());
    grExpParam_cpscRatio->SetTitle("Exponential Parameter vs H content (Capture/Scatter Ratio);H content (ppm);Exp decay parameter");
    grExpParam_cpscRatio->SetMarkerStyle(20);
    grExpParam_cpscRatio->SetMarkerColor(kBlue + 1);
    grExpParam_cpscRatio->SetLineColor(kBlue + 1);
    grExpParam_cpscRatio->Draw("APL");
    vCan.push_back(cExpParam_cpscRatio);

    TCanvas *ccpPosZ_TNcut = new TCanvas("cCpposZ_TNcut", "Capture Position Z by H ppm (TN cut)", 800, 600);
    TLegend *legcpPosZ_TNcut = new TLegend(0.55, 0.7, 0.88, 0.88);
    legcpPosZ_TNcut->SetNColumns(2);
    for (int i = 0; i < vHistcpPosZ_TNcut.size(); ++i)
    {
        vHistcpPosZ_TNcut[i]->SetLineColor(colorPalette[i]);
        vHistcpPosZ_TNcut[i]->SetLineWidth(2);
        // vHistcpPosZ_TNcut[i]->SetMaximum(40);
        vHistcpPosZ_TNcut[i]->SetMinimum(0);
        vHistcpPosZ_TNcut[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpPosZ_TNcut->AddEntry(vHistcpPosZ_TNcut[i], Form("%s", folder[i].Data()), "l");
    }
    legcpPosZ_TNcut->Draw();
    vCan.push_back(ccpPosZ_TNcut);

    // === フィットあり版 (TN cut, 別キャンバス) ===
    vector<TH1F *> vHistcpPosZ_TNcut_fit;
    vector<double> vExpParam_TNcut, vExpParamErr_TNcut, vPpm_TNcut;
    vector<TF1 *> vfExpFits_TNcut;
    TCanvas *ccpPosZ_TNcut_fit = new TCanvas("cCpposZ_TNcut_fit", "Capture Position Z by H ppm (TN cut, with fit)", 800, 600);
    TLegend *legcpPosZ_TNcut_fit = new TLegend(0.55, 0.7, 0.88, 0.88);
    legcpPosZ_TNcut_fit->SetNColumns(2);
    for (int i = 0; i < vHistcpPosZ_TNcut.size(); ++i)
    {
        TH1F *h_fit = (TH1F *)vHistcpPosZ_TNcut[i]->Clone(Form("h1_cpposZ_TNcut_fit_%d", i));
        h_fit->SetLineColor(colorPalette[i]);
        h_fit->SetLineWidth(2);
        h_fit->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpPosZ_TNcut_fit->AddEntry(h_fit, Form("%s", folder[i].Data()), "l");

        TF1 *fExp = new TF1(Form("fExp_TNcut_%d", i), "[0]*exp([1]*x)+[2]", 20, 100);
        fExp->SetParameters(h_fit->GetMaximum(), -0.05, 1.0);
        fExp->SetLineColor(kRed);
        fExp->SetLineWidth(2);
        fExp->SetNpx(1000);
        h_fit->Fit(fExp, "RQ+");
        fExp->Draw("same");

        vfExpFits_TNcut.push_back(fExp);
        double p1 = fExp->GetParameter(1);
        double p1_err = fExp->GetParError(1);
        vExpParam_TNcut.push_back(p1);
        vExpParamErr_TNcut.push_back(p1_err);

        TString ppmStr = folder[i];
        ppmStr.ReplaceAll("ppm", "");
        vPpm_TNcut.push_back(ppmStr.Atof());

        vHistcpPosZ_TNcut_fit.push_back(h_fit);
    }
    legcpPosZ_TNcut_fit->Draw();
    vCan.push_back(ccpPosZ_TNcut_fit);

    // === 指数パラメータ vs 含水率 (TN cut) ===
    TCanvas *cExpParam_TNcut = new TCanvas("cExpParam_TNcut", "Exponential Parameter vs H content (TN cut)", 800, 600);
    cExpParam_TNcut->SetLogx();
    TGraphErrors *grExpParam_TNcut = new TGraphErrors(vPpm_TNcut.size(), vPpm_TNcut.data(), vExpParam_TNcut.data(), nullptr, vExpParamErr_TNcut.data());
    grExpParam_TNcut->SetTitle("Exponential Parameter vs H content (TN cut);H content (ppm);Exp decay parameter");
    grExpParam_TNcut->SetMarkerStyle(20);
    grExpParam_TNcut->SetMarkerColor(kBlue + 1);
    grExpParam_TNcut->SetLineColor(kBlue + 1);
    grExpParam_TNcut->Draw("APL");
    vCan.push_back(cExpParam_TNcut);

    TCanvas *ccpscRatioZ_TNcut = new TCanvas("cCpScRatioZ_TNcut", "Capture/Scatter Ratio by Z (per H ppm)", 800, 600);
    TLegend *legcpscRatioZ_TNcut = new TLegend(0.55, 0.7, 0.88, 0.88);
    legcpscRatioZ_TNcut->SetNColumns(2);
    vector<TH1F *> vHistcpscRatioZ_TNcut;
    for (int i = 0; i < vHistcpPosZ_TNcut.size(); ++i)
    {
        TH1F *h1_cpScRatioZ_TNcut = (TH1F *)vHistcpPosZ_TNcut[i]->Clone(Form("h1_cpScRatioZ_%d", i));
        h1_cpScRatioZ_TNcut->Divide(vHistscPosZ[i]);
        h1_cpScRatioZ_TNcut->SetTitle("Capture/Scatter Ratio;Z (mm);Capture / Scatter");
        h1_cpScRatioZ_TNcut->SetLineColor(colorPalette[i]);
        h1_cpScRatioZ_TNcut->SetLineWidth(2);
        h1_cpScRatioZ_TNcut->SetMaximum(4);
        h1_cpScRatioZ_TNcut->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpscRatioZ_TNcut->AddEntry(h1_cpScRatioZ_TNcut, Form("%s", folder[i].Data()), "l");
        vHistcpscRatioZ_TNcut.push_back(h1_cpScRatioZ_TNcut);
    }
    legcpscRatioZ_TNcut->Draw();
    vCan.push_back(ccpscRatioZ_TNcut);

    // === フィットあり版 (Capture/Scatter Ratio TN cut, 別キャンバス) ===
    vector<TH1F *> vHistcpscRatioZ_TNcut_fit;
    vector<double> vExpParam_cpsc, vExpParamErr_cpsc, vPpm_cpsc;
    vector<TF1 *> vfExpFits_cpsc;
    TCanvas *ccpscRatioZ_TNcut_fit = new TCanvas("cCpScRatioZ_TNcut_fit", "Capture/Scatter Ratio by Z (per H ppm, TN cut, with fit)", 800, 600);
    TLegend *legcpscRatioZ_TNcut_fit = new TLegend(0.55, 0.7, 0.88, 0.88);
    legcpscRatioZ_TNcut_fit->SetNColumns(2);
    for (int i = 0; i < vHistcpscRatioZ_TNcut.size(); ++i)
    {
        TH1F *h_fit = (TH1F *)vHistcpscRatioZ_TNcut[i]->Clone(Form("h1_cpScRatioZ_TNcut_fit_%d", i));
        h_fit->SetLineColor(colorPalette[i]);
        h_fit->SetLineWidth(2);
        h_fit->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpscRatioZ_TNcut_fit->AddEntry(h_fit, Form("%s", folder[i].Data()), "l");

        TF1 *fExp = new TF1(Form("fExp_cpsc_%d", i), "[0]*exp([1]*x)+[2]", 20, 100);
        fExp->SetParameters(h_fit->GetMaximum(), -0.05, 1.0);
        fExp->SetLineColor(kRed);
        fExp->SetLineWidth(2);
        fExp->SetNpx(1000);
        h_fit->Fit(fExp, "RQ+");
        fExp->Draw("same");

        vfExpFits_cpsc.push_back(fExp);
        double p1 = fExp->GetParameter(1);
        double p1_err = fExp->GetParError(1);
        vExpParam_cpsc.push_back(p1);
        vExpParamErr_cpsc.push_back(p1_err);

        TString ppmStr = folder[i];
        ppmStr.ReplaceAll("ppm", "");
        vPpm_cpsc.push_back(ppmStr.Atof());

        vHistcpscRatioZ_TNcut_fit.push_back(h_fit);
    }
    legcpscRatioZ_TNcut_fit->Draw();
    vCan.push_back(ccpscRatioZ_TNcut_fit);

    // === 指数パラメータ vs 含水率 (Capture/Scatter Ratio TN cut) ===
    TCanvas *cExpParam_cpsc = new TCanvas("cExpParam_cpsc", "Exponential Parameter vs H content (Capture/Scatter Ratio, TN cut)", 800, 600);
    cExpParam_cpsc->SetLogx();
    // cExpParam_cpsc->SetLogy();
    TGraphErrors *grExpParam_cpsc = new TGraphErrors(vPpm_cpsc.size(), vPpm_cpsc.data(), vExpParam_cpsc.data(), nullptr, vExpParamErr_cpsc.data());
    grExpParam_cpsc->SetTitle("Exponential Parameter vs H content (Capture/Scatter Ratio, TN cut);H content (ppm);Exp decay parameter");
    grExpParam_cpsc->SetMarkerStyle(20);
    grExpParam_cpsc->SetMarkerColor(kBlue + 1);
    grExpParam_cpsc->SetLineColor(kBlue + 1);
    grExpParam_cpsc->Draw("APL");
    vCan.push_back(cExpParam_cpsc);

    TCanvas *ccpPosZ_TNcut_ratio = new TCanvas("cCpposZ_TNcut_ratio", "Capture Position Z by H ppm (Ratio)", 800, 600);
    TLegend *legcpPosZ_TNcut_ratio = new TLegend(0.55, 0.12, 0.88, 0.3);
    legcpPosZ_TNcut_ratio->SetNColumns(2);
    TH1F *h1_cpposZ_TNcut_0ppm = (TH1F *)vHistcpPosZ_TNcut[0]->Clone("h1_cpposZ_TNcut_ratio");
    vector<TH1F *> vHistcpPosZ_TNcut_ratio;
    for (int i = 0; i < vHistcpPosZ_TNcut.size(); ++i)
    {
        TH1F *h1_cpposZ_TNcut_ratio = (TH1F *)vHistcpPosZ_TNcut[i]->Clone(Form("h1_cpposZ_TNcut_ratio_%d", i));
        h1_cpposZ_TNcut_ratio->Divide(h1_cpposZ_TNcut_0ppm);
        vHistcpPosZ_TNcut_ratio.push_back(h1_cpposZ_TNcut_ratio);
        vHistcpPosZ_TNcut_ratio[i]->SetTitle("Capture Position Z Distribution Ratio (TN cut / 5e-7 eV <);Z (mm);Ratio to 0 ppm");
        vHistcpPosZ_TNcut_ratio[i]->SetLineColor(colorPalette[i]);
        vHistcpPosZ_TNcut_ratio[i]->SetLineWidth(2);
        vHistcpPosZ_TNcut_ratio[i]->SetMaximum(1.2);
        vHistcpPosZ_TNcut_ratio[i]->SetMinimum(0);
        vHistcpPosZ_TNcut_ratio[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legcpPosZ_TNcut_ratio->AddEntry(vHistcpPosZ_TNcut_ratio[i], Form("%s", folder[i].Data()), "l");
    }
    legcpPosZ_TNcut_ratio->Draw();
    vCan.push_back(ccpPosZ_TNcut_ratio);

    TCanvas *cscPosZ = new TCanvas("cScposZ", "Scatter Position Z by H ppm", 800, 600);
    TLegend *legscPosZ = new TLegend(0.55, 0.7, 0.88, 0.88);
    legscPosZ->SetNColumns(2);
    for (int i = 0; i < vHistscPosZ.size(); ++i)
    {
        vHistscPosZ[i]->SetLineColor(colorPalette[i]);
        vHistscPosZ[i]->SetLineWidth(2);
        vHistscPosZ[i]->SetMaximum(10);
        vHistscPosZ[i]->SetMinimum(2);
        vHistscPosZ[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legscPosZ->AddEntry(vHistscPosZ[i], Form("%s", folder[i].Data()), "l");
    }
    legscPosZ->Draw();
    vCan.push_back(cscPosZ);

    TCanvas *scPosZ_ratio = new TCanvas("cScposZ_ratio", "Scatter Position Z by H ppm (Ratio)", 800, 600);
    TLegend *legscPosZ_ratio = new TLegend(0.55, 0.12, 0.88, 0.3);
    legscPosZ_ratio->SetNColumns(2);
    TH1F *h1_scposZ_0ppm = (TH1F *)vHistscPosZ[0]->Clone("h1_scposZ_ratio");
    vector<TH1F *> vHistscPosZ_ratio;
    for (int i = 0; i < vHistscPosZ.size(); ++i)
    {
        TH1F *h1_scposZ_ratio = (TH1F *)vHistscPosZ[i]->Clone(Form("h1_scposZ_ratio_%d", i));
        h1_scposZ_ratio->Divide(h1_scposZ_0ppm);
        vHistscPosZ_ratio.push_back(h1_scposZ_ratio);
        vHistscPosZ_ratio[i]->SetTitle("Scatter Position Z Distribution Ratio;Z (mm);Ratio to 0 ppm");
        vHistscPosZ_ratio[i]->SetLineColor(colorPalette[i]);
        vHistscPosZ_ratio[i]->SetLineWidth(2);
        vHistscPosZ_ratio[i]->SetMaximum(1.2);
        vHistscPosZ_ratio[i]->SetMinimum(0.4);
        vHistscPosZ_ratio[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
        legscPosZ_ratio->AddEntry(vHistscPosZ_ratio[i], Form("%s", folder[i].Data()), "l");
    }
    legscPosZ_ratio->Draw();
    vCan.push_back(scPosZ_ratio);

    // save PDF
    if (true)
    {
        TString fPdfOut = "../fig/DetectPosition_Summary.pdf";
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
