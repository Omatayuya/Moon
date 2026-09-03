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

void DetectPosition()
{
    // vector<TString> folder = {"0ppm", "10ppm", "20ppm", "50ppm", "100ppm", "200ppm", "500ppm", "1000ppm", "2000ppm", "5000ppm", "10000ppm"};
    vector<TString> folder = {"0ppm"};

    for (int folderID = 0; folderID < folder.size(); folderID++)
    {

        // ==================================================================

        const string TargetParticle = "neutron";
        // const Double_t irrArea = pow(450, 2);                    // irradiation surface area (cm^2) for Proton, Helium
        const Double_t DetectorOffsetZ = 460;   // Detector offset in Z (mm)
        const Double_t irrArea = 600 * 600;     // irradiation surface area (cm^2) for Proton, Helium
        constexpr Double_t EJ270HalfWidth = 35; // EJ270 width (mm)

        constexpr double sideCut = 5; // mm, 20mm以内の範囲でのみカウントする

        // energy window
        constexpr double scatterEdepLow = 1.0; // MeV
        // constexpr double scatterEdepHigh = 3.0; // MeV
        constexpr double captureEdepLow = 4.5;  // MeV
        constexpr double captureEdepHigh = 5.0; // MeV

        // Thermal neutron cut (109Cd)
        constexpr double TNEnergyCut = 5e-7; // MeV

        constexpr double fidHalfWidth = EJ270HalfWidth - sideCut; // mm, 50mmの検出器のうち、20mm以内の範囲を除いた30mmの範囲でカウントする

        vector<double> cdPlanePosXY = {
            10 + DetectorOffsetZ,
            70 + DetectorOffsetZ,
        }; // XY平面 (法線: Z軸) の位置 [mm]
        vector<double> cdPlanePosZX = {-EJ270HalfWidth + sideCut, EJ270HalfWidth - sideCut}; // ZX平面 (法線: Y軸) の位置 [mm]
        vector<double> cdPlanePosZY = {-EJ270HalfWidth + sideCut, EJ270HalfWidth - sideCut}; // ZY平面 (法線: X軸) の位置 [mm]

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

        set<int> killedEvents;

        int nHits = HitTree->GetEntries();
        for (int i = 0; i < nHits; ++i)
        {
            HitTree->GetEntry(i);
            if (killedEvents.count(eventID) > 0)
            {
                continue;
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

            if (string(fpname) == "neutron")
            {
                bool cdAbsorbed = false;
                for (double zPlane : cdPlanePosXY) // XY平面 (法線: Z軸)
                {
                    if ((fPrePosZ - zPlane) * (fPostPosZ - zPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                    {
                        cdAbsorbed = true;
                    }
                }
                for (double yPlane : cdPlanePosZX) // ZX平面 (法線: Y軸)
                {
                    if ((fPrePosY - yPlane) * (fPostPosY - yPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                    {
                        cdAbsorbed = true;
                    }
                }
                for (double xPlane : cdPlanePosZY) // ZY平面 (法線: X軸)
                {
                    if ((fPrePosX - xPlane) * (fPostPosX - xPlane) < 0.0 && fPreKinE <= TNEnergyCut)
                    {
                        cdAbsorbed = true;
                    }
                }
                if (cdAbsorbed)
                {
                    killedEvents.insert(eventID);
                    continue;
                }
            }

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

        cout << "Virtual Cd layer absorbed: " << killedEvents.size() << " / " << nEvents << " events" << endl;

        for (const auto &entry : DetectorchamberMap)
        {
            eventChamberID = entry.first;
            chamberEventData = entry.second;
            DetectorchamberTree->Fill();
        }
        DetectorchamberMap.clear();
        const int entries = DetectorchamberTree->GetEntries();

        // ==================================================================
        int BinWidthZ = 5; // mm
        int minZ = 0;      // mm
        int maxZ = 80;     // mm
        int nBinsZ = (maxZ - minZ) / BinWidthZ;
        int BinWidthXY = 5; // mm
        int minXY = -35;    // mm
        int maxXY = 35;     // mm
        int nBinsXY = (maxXY - minXY) / BinWidthXY;

        // Match vH_ip_theta's log-uniform energy binning (see B2RunAction.cc)
        constexpr int nEnergyBins = 200;
        constexpr double energyMin = 1e-10; // MeV
        constexpr double energyMax = 1e4;   // MeV
        double energyBins[nEnergyBins + 1];
        for (int i = 0; i <= nEnergyBins; ++i)
            energyBins[i] = energyMin * pow(energyMax / energyMin, static_cast<double>(i) / nEnergyBins);

        // Thermal, epithermala and fast neutron bins
        constexpr int nEBins = 4;
        const double primEnergyEdges[nEBins + 1] = {0.0, 5e-7, 1e-3, 1.0, numeric_limits<double>::infinity()};
        const TString primEnergyLabels[nEBins] = {"E < 0.5 eV", "0.5 eV #leq E < 1 keV", "1 keV #leq E < 1 MeV", "E #geq 1 MeV"};
        const int primEnergyColors[nEBins] = {kOrange + 8, kGreen - 7, kGreen + 2, kBlue};

        // Sensor thickness bins
        vector<double> vSensThick = {0, 5, 10, 20, 40, 60, 70, 75, 80};
        const int nSensThick = vSensThick.size() - 1;
        vector<TString> zRegionLabels;
        for (int z = 0; z < nSensThick; ++z)
        {
            zRegionLabels.push_back(Form("%.0f <= Z < %.0f ", vSensThick[z], vSensThick[z + 1]));
        }
        constexpr double xFidCutLow = -(EJ270HalfWidth - sideCut); // mm, フィデューシャルカット下限
        constexpr double xFidCutHigh = EJ270HalfWidth - sideCut;   // mm, フィデューシャルカット上限
        // ==================================================================
        TH2D *h2_cpposZ = new TH2D("h2_cpposZ", Form("Capture Position Z; Energy (MeV); Z (mm); Count rate"), nEnergyBins, energyBins, nBinsZ, minZ, maxZ);
        TH1F *h1_cpposZ = new TH1F("h1_cpposZ", Form("Capture Position Z ;Z (mm);Count rate (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        vector<TH1F *> vh1_cpposZ_byE(nEBins);
        for (int e = 0; e < nEBins; ++e)
        {
            TString hname = Form("h1_cpposZ_E%d", e);
            vh1_cpposZ_byE[e] = new TH1F(hname,
                                         Form("Capture Position Z (%s);Z (mm);Count rate (s^{-1} %d mm^{-1})", primEnergyLabels[e].Data(), BinWidthZ),
                                         nBinsZ, minZ, maxZ);
            vh1_cpposZ_byE[e]->SetLineColor(primEnergyColors[e]);
        }
        TH2D *h2_cpposXY = new TH2D("h2_cpposXY", Form("Capture Position XY Distribution;X (mm);Y (mm);Count rate (s^{-1} %d #times %d mm^{-2})", BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        vector<TH2D *> vh2_cpposY(nSensThick);
        vector<TH1D *> vh1_cpposY(nSensThick);
        vector<vector<TH1D *>> vvh1_cpposY_byE(nSensThick, vector<TH1D *>(nEBins));
        vector<TH1D *> vh1_cpSpectrum(nSensThick);
        vector<TH1D *> vh1_cpSpectrum_ind(nSensThick);
        TArrayI savedPalette = TColor::GetPalette();
        gStyle->SetPalette(kRainBow);
        for (int z = 0; z < nSensThick; ++z)
        {
            TString hname_cpposEY = Form("h2_cpposEY_Z%d", z);
            vh2_cpposY[z] = new TH2D(hname_cpposEY,
                                     Form("Capture Position Y (%s, %.0f < X < %.0f );Energy (MeV);Y (mm);Count rate",
                                          zRegionLabels[z].Data(), xFidCutLow, xFidCutHigh),
                                     nEnergyBins, energyBins, nBinsXY, minXY, maxXY);

            TString hname_cpposY = Form("h1_cpposY_Z%d", z);
            vh1_cpposY[z] = new TH1D(hname_cpposY,
                                     Form("Capture Position Y (%s, %.0f < X < %.0f );Y (mm);Count rate (s^{-1} %d mm^{-1})",
                                          zRegionLabels[z].Data(), xFidCutLow, xFidCutHigh, BinWidthXY),
                                     nBinsXY, minXY, maxXY);

            TString hname_cpSpectrum = Form("h1_cpposEfid_Z%d", z);
            vh1_cpSpectrum[z] = new TH1D(hname_cpSpectrum,
                                         Form("Capture Energy (%s, %.0f < X < %.0f );Energy (MeV);Count rate (s^{-1})",
                                              zRegionLabels[z].Data(), xFidCutLow, xFidCutHigh),
                                         nEnergyBins, energyBins);
            vh1_cpSpectrum[z]->SetLineColor(TColor::GetColorPalette(
                static_cast<int>(1.0 * z * (TColor::GetNumberOfColors() - 1) / (nSensThick - 1))));

            TString hname_cpSpectrum_ind = Form("h1_cpposEfid_Z%d_ind", z);
            vh1_cpSpectrum_ind[z] = new TH1D(hname_cpSpectrum_ind,
                                             Form("Capture Energy (%s, %.0f < X < %.0f );Energy (MeV);Count rate (s^{-1})",
                                                  zRegionLabels[z].Data(), xFidCutLow, xFidCutHigh),
                                             nEnergyBins, energyBins);

            for (int e = 0; e < nEBins; ++e)
            {
                TString hname_cpposY_byE = Form("h1_cpposY_Z%d_E%d", z, e);
                vvh1_cpposY_byE[z][e] = new TH1D(hname_cpposY_byE,
                                                Form("Capture Position Y (%s, %s);Y (mm);Count rate (s^{-1} %d mm^{-1})",
                                                     zRegionLabels[z].Data(), primEnergyLabels[e].Data(), BinWidthXY),
                                                nBinsXY, minXY, maxXY);
                vvh1_cpposY_byE[z][e]->SetLineColor(primEnergyColors[e]);
            }
        }
        gStyle->SetPalette(savedPalette.GetSize(), savedPalette.GetArray());
        TH1F *h1_scposZ = new TH1F("h1_scposZ", Form("Scatter Position Z Distribution;Z (mm);Count rate (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        TH2D *h2_scposXY = new TH2D("h2_scposXY", Form("Scatter Position XY Distribution;X (mm);Y (mm);Count rate (s^{-1} %d#times%d mm^{-2})", BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);

        // vvHist: vector of vector of pairs of histogram and legend label
        using HistLegPair = pair<TH1 *, TString>;
        vector<vector<HistLegPair>> vvHist;
        vvHist.push_back({{h2_cpposZ, "Capture Position Z"}});
        {
            h1_cpposZ->SetLineColor(kBlack);
            vector<HistLegPair> group;
            group.push_back({h1_cpposZ, "Total"});
            for (int e = 0; e < nEBins; ++e)
                group.push_back({vh1_cpposZ_byE[e], primEnergyLabels[e]});
            vvHist.push_back(group);
        }
        vvHist.push_back({{h2_cpposXY, "Capture Position XY"}});
        for (int z = 0; z < nSensThick; ++z)
        {
            vvHist.push_back({{vh2_cpposY[z], Form("Capture Position EY (%s)", zRegionLabels[z].Data())}});

            vh1_cpposY[z]->SetLineColor(kBlack);
            vector<HistLegPair> group;
            group.push_back({vh1_cpposY[z], "Total"});
            for (int e = 0; e < nEBins; ++e)
            {
                group.push_back({vvh1_cpposY_byE[z][e], primEnergyLabels[e]});
            }
            vvHist.push_back({{vh1_cpSpectrum_ind[z], Form("Capture Energy (%s, %.0f < X < %.0f)", zRegionLabels[z].Data(), xFidCutLow, xFidCutHigh)}});
            vvHist.push_back(group);
        }
        {
            vh1_cpSpectrum[0]->SetTitle("Capture Energy by Z Region ;Energy (MeV);Count rate (s^{-1} bin^{-1})");
            vector<HistLegPair> group;
            for (int z = 0; z < nSensThick; ++z)
            {
                group.push_back({vh1_cpSpectrum[z], zRegionLabels[z]});
            }
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
                h2_cpposZ->Fill(eventChamberID.primEnergy, capturePosZ);

                // Z(position-axis) projection
                h1_cpposZ->Fill(capturePosZ);
                for (int e = 0; e < nEBins; ++e)
                {
                    if (eventChamberID.primEnergy >= primEnergyEdges[e] && eventChamberID.primEnergy < primEnergyEdges[e + 1])
                    {
                        vh1_cpposZ_byE[e]->Fill(capturePosZ);
                        break;
                    }
                }

                h2_cpposXY->Fill(capturePosX, capturePosY);

                // separate by thickness of sensitive layer
                for (int z = 0; z < nSensThick; ++z)
                {
                    if (capturePosZ >= vSensThick[z] && capturePosZ < vSensThick[z + 1])
                    {
                        if (capturePosX > xFidCutLow && capturePosX < xFidCutHigh)
                        {
                            vh2_cpposY[z]->Fill(eventChamberID.primEnergy, capturePosY);

                            // Y(position-axis) projection
                            vh1_cpposY[z]->Fill(capturePosY);
                            for (int e = 0; e < nEBins; ++e)
                            {
                                if (eventChamberID.primEnergy >= primEnergyEdges[e] && eventChamberID.primEnergy < primEnergyEdges[e + 1])
                                {
                                    vvh1_cpposY_byE[z][e]->Fill(capturePosY);
                                }
                            }
                            // Energy projection
                            vh1_cpSpectrum[z]->Fill(eventChamberID.primEnergy);
                            vh1_cpSpectrum_ind[z]->Fill(eventChamberID.primEnergy);
                        }
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

        for (int i = 0; i < vvHist.size(); ++i)
        {
            auto &vHist = vvHist[i];
            int nvHist = vHist.size();
            vCan.push_back(new TCanvas(Form("c%d", i + 2), vHist[0].first->GetTitle(), 800, 600));
            TLegend *legend = (nvHist > 1) ? new TLegend(0.55, 0.77, 0.88, 0.93) : nullptr;
            if (nvHist > 6)
                legend->SetNColumns(2);

            double yMax = 0;
            for (auto &entry : vHist)
            {
                auto *h = entry.first;
                h->Scale(1.0 / eqTime); // Convert to cps
                yMax = max(yMax, entry.first->GetMaximum());
            }

            bool hasTotal = (nvHist > 1 && vHist[0].second == "Total");
            double groupTotalInt = (nvHist > 1)
                                       ? vHist[0].first->Integral(0, vHist[0].first->GetNbinsX() + 1)
                                       : 0.0;

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
                    else if (std::find(vh2_cpposY.begin(), vh2_cpposY.end(), h2) != vh2_cpposY.end())
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
                    bool isEnergyProj = TString(h1->GetName()).BeginsWith("h1_cpposEfid_Z");
                    if (isEnergyProj)
                    {
                        gPad->SetLogx();
                        gPad->SetLogy();
                        h1->SetMaximum(yMax * 2);
                        h1->SetMinimum(5e-3);
                    }
                    else
                    {
                        h1->SetMaximum(yMax * 1.3);
                        h1->SetMinimum(0);
                    }
                    h1->Draw(j == 0 ? "HIST E" : "HIST E SAME");
                }

                if (legend)
                {
                    bool showFrac = (hasTotal && j > 0 && groupTotalInt > 0);
                    double frac = (j > 0 && groupTotalInt > 0)
                                      ? h->Integral(0, h->GetNbinsX() + 1) / groupTotalInt * 100.0
                                      : 0.0;
                    legend->AddEntry(h, showFrac ? Form("%s (%.1f%%)", title.Data(), frac) : title, "l");
                }
            }
            if (legend)
                legend->Draw();
        }

        // save PDF
        if (true)
        {
            TString fPdfOut = "../fig/" + folder[folderID] + "_DetectPosition_sideCut" + Form("%.0f", sideCut) + "mmt.pdf";
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
