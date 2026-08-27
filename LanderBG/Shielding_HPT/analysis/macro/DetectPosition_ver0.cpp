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
    TString strPhi = {"W_0"};
    // vector<TString> vStrPrimary = {"01_H", "02_He"};
    vector<TString> vStrPrimary = {"01-H"};
    const size_t nGCR = vStrPrimary.size(); // number of GCR nuclide

    // GCR flux
    vector<Double_t> vGcrFlux; // GCR 2pi flux (cm^-2 s^-1) @ 10MeV/n-100GeV/n
    constexpr Double_t gcrFactor = 4.;
    if (strPhi == "W_0")
    {
        vGcrFlux.emplace_back(4.93906 / gcrFactor); // 01-H
        vGcrFlux.emplace_back(0.45538 / gcrFactor); // 02-He
    }
    else if (strPhi == "W_118.5")
    {
        vGcrFlux.emplace_back(1.44867 / gcrFactor); // 01-H
        vGcrFlux.emplace_back(0.18490 / gcrFactor); // 02-He
    }
    else
    {
        cerr << "--- Phi error" << endl;
        return;
    }
    // if (vGcrFlux.size() != nGCR)
    // {
    //     cerr << "--- GCR nuclide number error" << endl;
    //     return;
    // }
    const Double_t irrArea = 5000 * 5000;

    for (int folderID = 0; folderID < vStrPrimary.size(); folderID++)
    {

        // ==================================================================

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

        for (auto axis : {"X", "Y", "Z"})
        {
            gStyle->SetLabelFont(62, axis);
            gStyle->SetTitleFont(62, axis);
            gStyle->SetTitleOffset(1.1, axis); // 軸タイトルのオフセット
            gStyle->SetLabelSize(0.04, axis);  // 目盛り数字のサイズ
            gStyle->SetTitleSize(0.04, axis);  // 軸タイトルのサイズ
        }
        gStyle->SetTextFont(62);
        gStyle->SetTitleFont(62, "");

        gStyle->SetPadGridX(true);
        gStyle->SetPadGridY(true);
        // gStyle->SetPalette(kRainBow);
        gStyle->SetOptStat(0);

        // ==================================================================

        TFile *fin = TFile::Open("../../Sim01/" + vStrPrimary[folderID] + "/" + strPhi + "/results.root");
        // TFile *fin = TFile::Open("../../Sim01/02-He/" + vStrPrimary[folderID] + "/results.root");
        if (!fin || fin->IsZombie())
        {
            cout << "Failed to open results.root" << endl;
            return;
        }

        TTree *HitTree = (TTree *)fin->Get("scorer");

        /*Input ParticleData*/
        TH1F *hPrimEnergy = (TH1F *)fin->Get("PrimaryEnergy");
        int nEvents = hPrimEnergy->GetEntries();

        const Double_t eqTime = static_cast<Double_t>(nEvents) / (vGcrFlux.at(folderID) * irrArea);

        // const Double_t eqTime = nEvents / (irrArea * GCRFlux);

        cout << "Total number of events in InputDataTree: " << nEvents << endl;
        cout << "Equivalent time: " << eqTime << " s" << endl;

        /* Output Data*/
        struct ChamberEventData
        {
            double neutronEntranceEnergy = 0.0;
            double edepSum = 0.0;
            double cpEdepSum = 0.0;
            double scEdepSum = 0.0;
            double cpPosX = 0.0, cpPosY = 0.0, cpPosZ = 0.0;
            double scPosX = 0.0, scPosY = 0.0, scPosZ = 0.0;
            double cpTriggerTime = numeric_limits<double>::infinity();
            double scTriggerTime = numeric_limits<double>::infinity();
            bool captureflag = false;
            bool scatterflag = false;
        };
        map<int, ChamberEventData> eventMap; // key: eventID

        // Must match EventAction::ParticleID() / ProcessID() in source/src/EventAction.cc
        constexpr int ID_proton = 0, ID_neutron = 1, ID_triton = 4, ID_alpha = 6;
        constexpr int ID_neutronInelastic = 16, ID_hadElastic = 31;

        int eventID = 0;
        vector<int> *parentID = nullptr, *parentParticleID = nullptr, *particleID = nullptr;
        vector<int> *stepProcessID = nullptr, *creatorProcessID = nullptr;
        vector<double> *parentEntranceEnergy = nullptr, *preEnergy = nullptr, *postEnergy = nullptr, *edep = nullptr, *gTime = nullptr;
        vector<double> *prePosX = nullptr, *prePosY = nullptr, *prePosZ = nullptr;
        vector<double> *postPosX = nullptr, *postPosY = nullptr, *postPosZ = nullptr;

        HitTree->SetBranchAddress("eventID", &eventID);
        HitTree->SetBranchAddress("parentEntranceEnergy", &parentEntranceEnergy);
        HitTree->SetBranchAddress("parentID", &parentID);
        HitTree->SetBranchAddress("parentParticleID", &parentParticleID);
        HitTree->SetBranchAddress("particleID", &particleID);
        HitTree->SetBranchAddress("stepProcessID", &stepProcessID);
        HitTree->SetBranchAddress("creatorProcessID", &creatorProcessID);
        HitTree->SetBranchAddress("preEnergy", &preEnergy);
        HitTree->SetBranchAddress("postEnergy", &postEnergy);
        HitTree->SetBranchAddress("edep", &edep);
        HitTree->SetBranchAddress("gTime", &gTime);
        HitTree->SetBranchAddress("prePosX", &prePosX);
        HitTree->SetBranchAddress("prePosY", &prePosY);
        HitTree->SetBranchAddress("prePosZ", &prePosZ);
        HitTree->SetBranchAddress("postPosX", &postPosX);
        HitTree->SetBranchAddress("postPosY", &postPosY);
        HitTree->SetBranchAddress("postPosZ", &postPosZ);

        // Get neutron entrance energy from the "neutronEntrance" tree
        TTree *EntranceTree = (TTree *)fin->Get("neutronEntrance");
        vector<double> *entranceEnergy = nullptr; // MeV, EJ270に入射した中性子全数(反応の有無によらず)
        EntranceTree->SetBranchAddress("energy", &entranceEnergy);

        int nEntries = HitTree->GetEntries();
        for (int i = 0; i < nEntries; ++i)
        {
            HitTree->GetEntry(i);

            auto &acc = eventMap[eventID];
            size_t nHits = parentID->size();

            for (size_t j = 0; j < nHits; ++j)
            {
                acc.neutronEntranceEnergy = parentEntranceEnergy->at(j);
                acc.edepSum += edep->at(j);

                /* screening conditions for Capture & Scatter event*/
                const bool isCaptureParticle =
                    (particleID->at(j) == ID_alpha || particleID->at(j) == ID_triton);
                const bool isScatterParticle = (particleID->at(j) == ID_proton);

                // 親トラックがEJ-270に入射した時点でneutronだったか(proton-proton散乱を除外するため)
                const bool parentIsNeutron = (parentParticleID->at(j) == ID_neutron);

                const bool isCaptureDepositStep =
                    (edep->at(j) > 0.0) &&
                    isCaptureParticle &&
                    (creatorProcessID->at(j) == ID_neutronInelastic); // 6Li(n,a)t は neutronInelastic 側の反応モデル
                const bool isScatterDepositStep =
                    (edep->at(j) > 0.0) &&
                    isScatterParticle &&
                    parentIsNeutron &&
                    (creatorProcessID->at(j) == ID_hadElastic);

                if (isCaptureDepositStep)
                {
                    acc.cpEdepSum += edep->at(j);

                    if (acc.captureflag == false && acc.cpEdepSum > captureEdepLow)
                    {
                        acc.captureflag = true;
                        acc.cpTriggerTime = min(acc.cpTriggerTime, gTime->at(j));
                        acc.cpPosX = prePosX->at(j); // already in mm
                        acc.cpPosY = prePosY->at(j); // already in mm
                        acc.cpPosZ = prePosZ->at(j); // already in mm
                    }
                }

                if (isScatterDepositStep)
                {
                    acc.scEdepSum += edep->at(j);
                    if (acc.scatterflag == false && acc.scEdepSum > scatterEdepLow)
                    {
                        acc.scatterflag = true;
                        acc.scTriggerTime = min(acc.scTriggerTime, gTime->at(j));
                        acc.scPosX = prePosX->at(j); // already in mm
                        acc.scPosY = prePosY->at(j); // already in mm
                        acc.scPosZ = prePosZ->at(j); // already in mm
                    }
                }
            }
        }

        int BinWidthZ = 5; // mm
        int minZ = 0;      // mm
        int maxZ = 200;    // mm
        int nBinsZ = (maxZ - minZ) / BinWidthZ;
        int BinWidthXY = 5; // mm
        int minXY = -75;    // mm
        int maxXY = 75;     // mm
        int nBinsXY = (maxXY - minXY) / BinWidthXY;

        // Match vH_ip_theta's log-uniform energy binning (see B2RunAction.cc)
        constexpr int nEnergyBins = 500;
        constexpr double energyMin = 1e-10; // MeV
        constexpr double energyMax = 1e4;   // MeV
        double energyBins[nEnergyBins + 1];
        for (int i = 0; i <= nEnergyBins; ++i)
            energyBins[i] = energyMin * pow(energyMax / energyMin, static_cast<double>(i) / nEnergyBins);

        TH1F *h1_neutronEntranceEnergy = new TH1F("h1_neutronEntranceEnergy", "Neutron Entrance Energy;Energy (MeV);Counts (s^{-1})", nEnergyBins, energyBins);
        TH1F *h1_cpposZ = new TH1F("h1_cpposZ", Form("Capture Position Z ;Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        TH1F *h1_cpposZ_TNcut = new TH1F("h1_cpposZ_TNcut", Form("Capture Position Z  (TN cut / %.1f eV <);Z (mm);Counts (s^{-1} %d mm^{-1})", TNEnergyCut * 1e6, BinWidthZ), nBinsZ, minZ, maxZ);
        TH1F *h1_cpposZ_sidecut = new TH1F("h1_cpposZ_sidecut", Form("Capture Position Z  (Side 10 mmt cut);Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        TH2D *h2_cpposXY = new TH2D("h2_cpposXY", Form("Capture Position XY Distribution;X (mm);Y (mm);Counts (s^{-1} %d #times %d mm^{-2})", BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        TH2D *h2_cpposXY_TNcut = new TH2D("h2_cpposXY_TNcut", Form("Capture Position XY Distribution (TN cut / %.1f eV <);X (mm);Y (mm);Counts (s^{-1} %d #times %d mm^{-2})", TNEnergyCut * 1e6, BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        TH1F *h1_scposZ = new TH1F("h1_scposZ", Form("Scatter Position Z Distribution;Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
        TH2D *h2_scposXY = new TH2D("h2_scposXY", Form("Scatter Position XY Distribution;X (mm);Y (mm);Counts (s^{-1} %d #times %d mm^{-2})", BinWidthXY, BinWidthXY), nBinsXY, minXY, maxXY, nBinsXY, minXY, maxXY);
        vector<TH1 *> vHist;
        // vHist.push_back(h1_neutronEntranceEnergy);
        // vHist.push_back(h1_cpposZ);
        vHist.push_back(h1_cpposZ_TNcut);
        // vHist.push_back(h1_cpposZ_sidecut);
        vHist.push_back(h2_cpposXY);
        vHist.push_back(h2_cpposXY_TNcut);
        vHist.push_back(h1_scposZ);
        vHist.push_back(h2_scposXY);

        vector<TH1F *> vH_cpposZ_byE(nEBins);
        for (int e = 0; e < nEBins; ++e)
        {
            TString hname = Form("h1_cpposZ_E%d", e);
            vH_cpposZ_byE[e] = new TH1F(hname,
                                        Form("Capture Position Z (%s);Z (mm);Counts (s^{-1} %d mm^{-1})", primEnergyLabels[e].Data(), BinWidthZ),
                                        nBinsZ, minZ, maxZ);
        }

        for (const auto &kv : eventMap)
        {
            const auto &chamberEventData = kv.second;

            double capturePosZ = chamberEventData.cpPosZ + 100; // Convert m to mm and subtract detector offset
            double capturePosX = chamberEventData.cpPosX;       // Convert m to mm
            double capturePosY = chamberEventData.cpPosY;       // Convert m to mm
            if (chamberEventData.captureflag)
            {
                h1_cpposZ->Fill(capturePosZ);
                h2_cpposXY->Fill(capturePosX, capturePosY);
                if (chamberEventData.neutronEntranceEnergy > TNEnergyCut) // TN cut
                {
                    h1_cpposZ_TNcut->Fill(capturePosZ);
                    h2_cpposXY_TNcut->Fill(capturePosX, capturePosY);
                }

                if (capturePosX < -(50 - 10) || capturePosX > (50 - 10) || capturePosY < -(50 - 10) || capturePosY > (50 - 10)) // Side cut
                {
                    h1_cpposZ_sidecut->Fill(capturePosZ);
                }

                // neutronEntranceEnergyで分岐してFill
                for (int e = 0; e < nEBins; ++e)
                {
                    if (chamberEventData.neutronEntranceEnergy >= primEnergyEdges[e] && chamberEventData.neutronEntranceEnergy < primEnergyEdges[e + 1])
                    {
                        vH_cpposZ_byE[e]->Fill(capturePosZ);
                        break;
                    }
                }
            }

            double scatterPosZ = chamberEventData.scPosZ + 100; // Convert m to mm and subtract detector offset
            double scatterPosX = chamberEventData.scPosX;       // Convert m to mm
            double scatterPosY = chamberEventData.scPosY;       // Convert
            if (chamberEventData.scatterflag)
            {
                h1_scposZ->Fill(scatterPosZ);
                h2_scposXY->Fill(scatterPosX, scatterPosY);
            }
        }

        // neutronEntranceEnergy
        int nEntranceEntries = EntranceTree->GetEntries();
        for (int i = 0; i < nEntranceEntries; ++i)
        {
            EntranceTree->GetEntry(i);
            for (double e : *entranceEnergy)
                h1_neutronEntranceEnergy->Fill(e);
        }

        // Draw
        vector<TCanvas *> vCan;
        TCanvas *cEntranceEnergy = new TCanvas("cEntranceEnergy", "Neutron Entrance Energy", 800, 600);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        h1_neutronEntranceEnergy->Scale(1.0 / eqTime); // 他のヒストグラムと同様にcps換算
        h1_neutronEntranceEnergy->Draw("HIST E");
        vCan.push_back(cEntranceEnergy);

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
        h1_cpposZ->SetMinimum(1e-5);
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

        for (int i = 0; i < vHist.size(); ++i)
        {
            vCan.push_back(new TCanvas(Form("c%d", i + 2), vHist[i]->GetTitle(), 800, 600));
            if (auto h2 = dynamic_cast<TH2 *>(vHist[i]))
            {
                h2->Scale(1.0 / eqTime); // Convert to cps
                h2->Draw("COLZ");
            }
            else if (auto h1 = dynamic_cast<TH1 *>(vHist[i]))
            {
                gPad->SetGridx();
                gPad->SetGridy();
                h1->Scale(1.0 / eqTime); // Convert to cps
                h1->Draw("HIST E");
            }
        }

        // save PDF
        if (true)
        {
            TString fPdfOut = "../fig/" + vStrPrimary[folderID] + "_" + strPhi + "_DetectPosition.pdf";
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
