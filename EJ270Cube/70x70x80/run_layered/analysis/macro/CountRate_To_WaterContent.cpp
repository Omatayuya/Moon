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
#include <regex>

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
#include <TGraphSmooth.h>

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
    bool useFidcut = true;                                // Use fiducial cut for capture count rate
    constexpr Double_t EJ270HalfWidth = 35;               // mm
    constexpr double sideCut = 5;                         // mm
    const double fidHalfWidth = EJ270HalfWidth - sideCut; // mm

    bool useCd = true;                                                                                // Use fiducial cut for capture count rate
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

    // run_layered: 一様含水率と同時にプロットする含水層の厚みと、横軸の定義
    const TString targetThickStr = "0.05";  // m, folder名 box_10x10x{thickness}m の thickness 部分と一致させる
    const double uniformRefThickness = 1.5; // m, 一様含水率リファレンス層の厚み (DetectPosition_Summary.cpp と同じ)
    bool useEquivPpm = true;                // true: 横軸を一様換算 ppm (ppm_layer x thickness / 1.5 m), false: 含水層の ppm

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

    // 一様含水率 (run/) と、厚み targetThickStr の含水層 (run_layered/) の結果ファイル一覧
    vector<TString> vPath;
    vector<TString> vFileDepthStr; // 一様含水率は "uniform"
    vector<double> vFileDepth;     // 一様含水率は -1
    vector<double> vFilePpm;       // 横軸に使う含水率 (ppm)

    vector<TString> uniformFolder = {"0ppm", "10ppm", "20ppm", "50ppm", "100ppm", "200ppm", "500ppm", "1000ppm", "2000ppm", "5000ppm", "10000ppm"};
    for (const auto &f : uniformFolder)
    {
        TString ppmStr = f(0, f.Length() - 3);
        vPath.push_back("../../../run/" + f + "/results.root");
        vFileDepthStr.push_back("uniform");
        vFileDepth.push_back(-1);
        vFilePpm.push_back(ppmStr.Atof());
    }

    {
        ifstream ifsFolder("../folders.list");
        if (!ifsFolder)
        {
            cerr << "--- ../folders.list not found. Run run_layered/setupDirs.sh first." << endl;
            return;
        }
        // folder名 "box_10x10x{thickness}m_depth_{depth}m_H_{ppm}ppm"
        std::regex folderNameRe("box_10x10x([0-9.]+)m_depth_([0-9.]+)m_H_([0-9.]+)ppm");
        string line;
        while (getline(ifsFolder, line))
        {
            std::smatch match;
            if (line.empty() || !std::regex_search(line, match, folderNameRe))
                continue;
            if (TString(match[1].str()) != targetThickStr)
                continue;

            double thickness = stod(match[1].str());
            double ppmLayer = stod(match[3].str());
            vPath.push_back("../../" + TString(line) + "/results.root");
            vFileDepthStr.push_back(match[2].str());
            vFileDepth.push_back(stod(match[2].str()));
            vFilePpm.push_back(useEquivPpm ? ppmLayer * thickness / uniformRefThickness : ppmLayer);
        }
    }

    vector<TH1F *> vHistcpPosZ;
    vector<TH1F *> vHistcpPosZ_TNcut;
    vector<TH1F *> vHistscPosZ;
    vector<double> vEqTime;

    for (int folderID = 0; folderID < vPath.size(); folderID++)
    {

        TFile *fin = TFile::Open(vPath[folderID]);
        if (!fin || fin->IsZombie())
        {
            cout << "Failed to open " << vPath[folderID] << endl;
            return;
        }
        cout << "Opened " << vPath[folderID] << endl;

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

            if (string(fpname) == "neutron" && useCd)
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
                if (useFidcut && (abs(fPrePosX) > fidHalfWidth || abs(fPrePosY) > fidHalfWidth))
                {
                    continue; // Skip this step if outside the fiducial volume
                }

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
        h1_cpposZ->SetDirectory(nullptr);
        h1_cpposZ_TNcut->SetDirectory(nullptr);
        h1_scposZ->SetDirectory(nullptr);

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

        delete DetectorchamberTree;
        fin->Close();
        delete fin;
    }

    // === 各ファイルの検出層ごとの capture レートと、全 scatter に対する比率 ===
    vector<double> vSensThick = {0, 5, 10, 20, 40}; // mm、DetectPosition.cpp と同じ流儀
    const int nSensThick = vSensThick.size() - 1;
    vector<TString> zRegionLabels;
    for (int z = 0; z < nSensThick; ++z)
        zRegionLabels.push_back(Form("%.0f #leq Z < %.0f mm", vSensThick[z], vSensThick[z + 1]));

    vector<double> vSc, vScErr;                                                 // [file] 全 scatter レート (s^-1)
    vector<vector<double>> vCpLayer(nSensThick), vCpLayerErr(nSensThick);       // [layer][file] capture レート (s^-1)
    vector<vector<double>> vRatioLayer(nSensThick), vRatioLayerErr(nSensThick); // [layer][file] capture_layer / scatter_total
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
            double ratioErr = (ratio > 0) ? ratio * sqrt(pow(capErr / capSum, 2) + pow(scatterErr / scatterSum, 2)) : 0.0;
            vRatioLayer[z].push_back(ratio);
            vRatioLayerErr[z].push_back(ratioErr);
        }
    }

    const size_t iZero = 0; // run/0ppm (一様・乾燥) を全曲線共通の 0 ppm 点とする

    // === 曲線の定義: [0] 一様含水率, [1..] 含水層の深さごと。各曲線は 含水率 -> ファイル番号 ===
    map<double, TString> depthMap; // 深さ (数値順) -> 深さ文字列
    for (size_t i = 0; i < vPath.size(); ++i)
    {
        if (vFileDepth[i] >= 0)
            depthMap[vFileDepth[i]] = vFileDepthStr[i];
    }
    const int nDepth = depthMap.size();

    vector<TString> vCurveKey = {"uniform"};
    vector<TString> vCurveLabel = {"Uniform"};
    for (const auto &depthEntry : depthMap)
    {
        vCurveKey.push_back(depthEntry.second);
        vCurveLabel.push_back("depth " + depthEntry.second + " m");
    }

    vector<map<double, size_t>> vCurvePoints(vCurveKey.size());
    for (size_t c = 0; c < vCurveKey.size(); ++c)
    {
        vCurvePoints[c][0.0] = iZero;
        for (size_t i = 0; i < vPath.size(); ++i)
        {
            if (vFileDepthStr[i] == vCurveKey[c])
                vCurvePoints[c][vFilePpm[i]] = i;
        }
    }

    // 横軸の範囲とタイトル
    double xmin = numeric_limits<double>::infinity();
    double xmax = 2e3;
    for (size_t i = 0; i < vFilePpm.size(); ++i)
    {
        if (vFilePpm[i] > 0)
            xmin = min(xmin, vFilePpm[i] / 2);
        // xmax = max(xmax, vFilePpm[i] * 2);
    }
    const TString xTitle = useEquivPpm ? Form("H content (ppm, %.1f m uniform equivalent)", uniformRefThickness) : "H content (ppm)";
    const TString layerTag = "uniform vs layer " + targetThickStr + " m";

    gStyle->SetPalette(kRainBow);
    int nColors = gStyle->GetNumberOfColors();

    vector<TCanvas *> vCan;

    // --- 全 scatter レート vs 含水率 ---
    {
        TMultiGraph *mgSc = new TMultiGraph();
        mgSc->SetTitle("Scatter count rate (total, " + layerTag + ");" + xTitle + ";Count rate (s^{-1})");
        TLegend *legSc = new TLegend(0.15, 0.15, 0.60, 0.40);
        legSc->SetNColumns(2);

        for (size_t c = 0; c < vCurveKey.size(); ++c)
        {
            vector<double> x, y, yErr;
            for (const auto &p : vCurvePoints[c])
            {
                x.push_back(p.first);
                y.push_back(vSc.at(p.second));
                yErr.push_back(vScErr.at(p.second));
            }

            bool isUniform = (c == 0);
            int d = (int)c - 1;
            Color_t col = isUniform ? kBlack : gStyle->GetColorPalette(nDepth > 1 ? static_cast<int>((0.1 + 0.8 * d / (nDepth - 1)) * (nColors - 1)) : 0);
            TGraphErrors *gr = new TGraphErrors(x.size(), x.data(), y.data(), 0, yErr.data());
            gr->SetMarkerColor(col);
            gr->SetLineColor(col);
            gr->SetLineWidth(isUniform ? 3 : 2);
            gr->SetMarkerStyle(isUniform ? 20 : 21 + (d % 9));
            mgSc->Add(gr, "PL");
            legSc->AddEntry(gr, vCurveLabel[c], "lp");
        }

        TCanvas *cSc = new TCanvas("cSc", "Scatter count rate vs H content", 800, 600);
        cSc->SetLogx();
        mgSc->Draw("A");
        mgSc->GetXaxis()->SetLimits(xmin, xmax);
        legSc->Draw();
        vCan.push_back(cSc);
    }

    // --- 検出層ごと: capture / 比率 / 0 ppm 規格化 / 分離に必要な観測時間 ---
    for (int z = 0; z < nSensThick; ++z)
    {
        TString titleTag = zRegionLabels[z];

        TMultiGraph *mgCp = new TMultiGraph();
        mgCp->SetTitle("Capture count rate (" + titleTag + ");" + xTitle + ";Count rate (s^{-1})");
        TLegend *legCp = new TLegend(0.15, 0.15, 0.60, 0.40);
        legCp->SetNColumns(2);

        TMultiGraph *mgCp_0ppm = new TMultiGraph();
        mgCp_0ppm->SetTitle("Capture count rate normalized to 0 ppm (" + titleTag + ");" + xTitle + ";Relative count rate");
        TLegend *legCp_0ppm = new TLegend(0.15, 0.15, 0.60, 0.38);
        legCp_0ppm->SetNColumns(2);

        TMultiGraph *mgRatio = new TMultiGraph();
        mgRatio->SetTitle("Count rate ratio capture / scatter total (" + titleTag + ");" + xTitle + ";Count rate ratio");
        TLegend *legRatio = new TLegend(0.15, 0.15, 0.60, 0.38);
        legRatio->SetNColumns(2);

        TMultiGraph *mgRatio_0ppm = new TMultiGraph();
        mgRatio_0ppm->SetTitle("Count rate ratio normalized to 0 ppm (" + titleTag + ");" + xTitle + ";Relative count rate ratio");
        TLegend *legRatio_0ppm = new TLegend(0.15, 0.15, 0.60, 0.38);
        legRatio_0ppm->SetNColumns(2);

        TMultiGraph *mgSigTime = new TMultiGraph();
        mgSigTime->SetTitle(Form("Observation time for %.0f#sigma separation from 0 ppm (%s);%s;Observation time (s)", nSigma, titleTag.Data(), xTitle.Data()));
        TLegend *legSigTime = new TLegend(0.15, 0.15, 0.88, 0.38);
        legSigTime->SetNColumns(2);

        // 0 ppm 点 (iZero) の値
        const double yCp0 = vCpLayer[z].at(iZero), yCp0Err = vCpLayerErr[z].at(iZero);
        const double y0 = vRatioLayer[z].at(iZero), y0Err = vRatioLayerErr[z].at(iZero);
        const double rateSc0 = vSc.at(iZero);

        for (size_t c = 0; c < vCurveKey.size(); ++c)
        {
            vector<double> x, yCp, yCpErr, yCpNorm, yCpNormErr, yRatio, yRatioErr, yNorm, yNormErr;
            vector<double> vTimePpm, vTimeReq;
            for (const auto &p : vCurvePoints[c])
            {
                const size_t j = p.second;
                x.push_back(p.first);

                yCp.push_back(vCpLayer[z].at(j));
                yCpErr.push_back(vCpLayerErr[z].at(j));
                double cpNorm = (yCp0 > 0) ? vCpLayer[z].at(j) / yCp0 : 0.0;
                yCpNorm.push_back(cpNorm);
                yCpNormErr.push_back((cpNorm > 0) ? cpNorm * sqrt(pow(vCpLayerErr[z].at(j) / vCpLayer[z].at(j), 2) + pow(yCp0Err / yCp0, 2)) : 0.0);

                yRatio.push_back(vRatioLayer[z].at(j));
                yRatioErr.push_back(vRatioLayerErr[z].at(j));
                double ratioNorm = (y0 > 0) ? vRatioLayer[z].at(j) / y0 : 0.0;
                yNorm.push_back(ratioNorm);
                yNormErr.push_back((ratioNorm > 0) ? ratioNorm * sqrt(pow(vRatioLayerErr[z].at(j) / vRatioLayer[z].at(j), 2) + pow(y0Err / y0, 2)) : 0.0);

                // ratioErr(T)^2 = ratio^2 (1/rateCap + 1/rateSc) / T   (ポアソン統計、T は観測時間 s)
                double ratioX = vRatioLayer[z].at(j);
                double rateCapX = vCpLayer[z].at(j);
                double rateScX = vSc.at(j);
                double delta = fabs(ratioX - y0);
                if (j == iZero || y0 <= 0 || yCp0 <= 0 || rateSc0 <= 0 || delta <= 0 || ratioX <= 0 || rateCapX <= 0)
                    continue;
                double C0 = pow(y0, 2) * (1.0 / rateSc0 + 1.0 / yCp0);
                double CX = pow(ratioX, 2) * (1.0 / rateScX + 1.0 / rateCapX);
                vTimePpm.push_back(p.first);
                vTimeReq.push_back(pow(nSigma, 2) * (C0 + CX) / pow(delta, 2)); // s
            }

            bool isUniform = (c == 0);
            int d = (int)c - 1;
            Color_t col = isUniform ? kGray : gStyle->GetColorPalette(nDepth > 1 ? static_cast<int>((0.1 + 0.8 * d / (nDepth - 1)) * (nColors - 1)) : 0);
            int mstyle = isUniform ? 20 : 21 + (d % 9);
            int lstyle = isUniform ? 2 : 1;
            int lwidth = 2;

            TGraphErrors *grCp = new TGraphErrors(x.size(), x.data(), yCp.data(), 0, yCpErr.data());
            grCp->SetMarkerColor(col);
            grCp->SetLineColor(col);
            grCp->SetLineWidth(lwidth);
            grCp->SetLineStyle(lstyle);
            grCp->SetMarkerStyle(mstyle);
            mgCp->Add(grCp, "PL");
            legCp->AddEntry(grCp, vCurveLabel[c], "lp");

            TGraphErrors *grCpNorm = new TGraphErrors(x.size(), x.data(), yCpNorm.data(), 0, yCpNormErr.data());
            grCpNorm->SetMarkerColor(col);
            grCpNorm->SetLineColor(col);
            grCpNorm->SetLineWidth(lwidth);
            grCpNorm->SetLineStyle(lstyle);
            grCpNorm->SetMarkerStyle(mstyle);
            mgCp_0ppm->Add(grCpNorm, "PL");
            legCp_0ppm->AddEntry(grCpNorm, vCurveLabel[c], "lp");

            TGraphErrors *grRatio = new TGraphErrors(x.size(), x.data(), yRatio.data(), 0, yRatioErr.data());
            grRatio->SetMarkerColor(col);
            grRatio->SetLineColor(col);
            grRatio->SetLineWidth(lwidth);
            grRatio->SetLineStyle(lstyle);
            grRatio->SetMarkerStyle(mstyle);
            mgRatio->Add(grRatio, "PL");
            legRatio->AddEntry(grRatio, vCurveLabel[c], "lp");

            TGraphErrors *grNorm = new TGraphErrors(x.size(), x.data(), yNorm.data(), 0, yNormErr.data());
            grNorm->SetMarkerColor(col);
            grNorm->SetLineColor(col);
            grNorm->SetLineWidth(lwidth);
            grNorm->SetLineStyle(lstyle);
            grNorm->SetMarkerStyle(mstyle);
            mgRatio_0ppm->Add(grNorm, "PL");
            legRatio_0ppm->AddEntry(grNorm, vCurveLabel[c], "lp");

            if (!vTimePpm.empty())
            {
                TGraph *grTime = new TGraph(vTimePpm.size(), vTimePpm.data(), vTimeReq.data());
                grTime->SetMarkerColor(col);
                grTime->SetLineColor(col);
                grTime->SetLineWidth(lwidth);
                grTime->SetLineStyle(lstyle);
                grTime->SetMarkerStyle(mstyle);
                mgSigTime->Add(grTime, "PL");
                legSigTime->AddEntry(grTime, vCurveLabel[c], "lp");
            }
        }

        TCanvas *cCp = new TCanvas(Form("cCp_Z%d", z), "Capture count rate vs H content", 800, 600);
        cCp->SetLogx();
        cCp->SetLogy();
        mgCp->Draw("A");
        mgCp->GetXaxis()->SetLimits(xmin, xmax);
        legCp->Draw();
        vCan.push_back(cCp);

        TCanvas *cCp_0ppm = new TCanvas(Form("cCp_0ppm_Z%d", z), "Capture count rate (normalized to 0 ppm)", 800, 600);
        cCp_0ppm->SetLogx();
        mgCp_0ppm->SetMinimum(0);
        mgCp_0ppm->Draw("A");
        mgCp_0ppm->GetXaxis()->SetLimits(xmin, xmax);
        legCp_0ppm->Draw();
        vCan.push_back(cCp_0ppm);

        TCanvas *cRatio = new TCanvas(Form("cRatio_Z%d", z), "Count rate ratio vs H content", 800, 600);
        cRatio->SetLogx();
        mgRatio->SetMinimum(0);
        mgRatio->Draw("A");
        mgRatio->GetXaxis()->SetLimits(xmin, xmax);
        legRatio->Draw();
        vCan.push_back(cRatio);

        TCanvas *cRatio_0ppm = new TCanvas(Form("cRatio_0ppm_Z%d", z), "Count rate ratio (normalized to 0 ppm)", 800, 600);
        cRatio_0ppm->SetLogx();
        mgRatio_0ppm->SetMinimum(0);
        mgRatio_0ppm->Draw("A");
        mgRatio_0ppm->GetXaxis()->SetLimits(xmin, xmax);
        legRatio_0ppm->Draw();
        vCan.push_back(cRatio_0ppm);

        if (mgSigTime->GetListOfGraphs())
        {
            TCanvas *cSigTime = new TCanvas(Form("cSigTime_Z%d", z), "Observation time for Nsigma separation vs H content", 800, 600);
            cSigTime->SetLogx();
            cSigTime->SetLogy();
            cSigTime->SetGridx(0);
            cSigTime->SetGridy(0);
            mgSigTime->Draw("A");
            mgSigTime->GetXaxis()->SetLimits(xmin, xmax);
            legSigTime->Draw();

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
            vCan.push_back(cSigTime);
        }
    }

    // save PDF
    if (true)
    {
        TString fPdfOut = "../fig/CountRate_To_WaterContent_uniform_vs_layer" + targetThickStr + "m" + (useEquivPpm ? "_equivPpm" : "_layerPpm") + ".pdf";
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
        cout << "--- " << fPdfOut << " was created." << endl;
    }

    return;
}
