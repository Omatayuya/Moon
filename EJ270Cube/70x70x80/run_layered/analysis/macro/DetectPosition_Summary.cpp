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
#include <regex>
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
#include <TSystem.h>

void DetectPosition_Summary()
{
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

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
    gStyle->SetPadLeftMargin(0.15);

    vector<TString> folder;
    {
        ifstream ifsFolder("../folders.list");
        if (!ifsFolder)
        {
            cerr << "--- ../folders.list not found. Run run_layered/setupDirs.sh first." << endl;
            return;
        }
        string line;
        while (getline(ifsFolder, line))
        {
            if (!line.empty())
                folder.push_back(TString(line));
        }
    }

    // folder名 "box_10x10x{thickness}m_depth_{depth}m_H_{ppm}ppm" を解析する正規表現
    std::regex folderNameRe("box_10x10x([0-9.]+)m_depth_([0-9.]+)m_H_([0-9.]+)ppm");

    struct FolderInfo
    {
        TString name;
        TString ppmStr;
        TString thicknessStr;
        TString depthStr;
        double ppm = 0;
        double thickness = 0;
        double depth = 0;
    };

    // (ppmStr, thicknessStr) ごとに depth 違いのフォルダをまとめる = by_H_total の各グループに対応
    map<pair<string, string>, vector<FolderInfo>> groupMap;
    for (const auto &f : folder)
    {
        std::smatch match;
        std::string fs = f.Data();
        if (!std::regex_search(fs, match, folderNameRe))
        {
            cerr << "Warning: failed to parse folder name: " << fs << endl;
            continue;
        }

        FolderInfo fi;
        fi.name = f;
        fi.thicknessStr = match[1].str();
        fi.depthStr = match[2].str();
        fi.ppmStr = match[3].str();
        fi.thickness = std::stod(fi.thicknessStr.Data());
        fi.depth = std::stod(fi.depthStr.Data());
        fi.ppm = std::stod(fi.ppmStr.Data());

        groupMap[{fi.ppmStr.Data(), fi.thicknessStr.Data()}].push_back(fi);
    }

    const double uniformRefThickness = 1.5; // m, 一様含水率リファレンス層の厚み

    for (auto &group : groupMap)
    {
        vector<FolderInfo> &depthList = group.second;
        sort(depthList.begin(), depthList.end(), [](const FolderInfo &a, const FolderInfo &b)
             { return a.depth < b.depth; });

        const FolderInfo &sample = depthList.front();
        TString groupLabel = "H_" + sample.ppmStr + "ppm_thickness_" + sample.thicknessStr + "m";

        // 一様含水率の参照ppmを計算: ppm_layer x thickness_layer = ppm_uniform x uniformRefThickness
        double ppmUniform = sample.ppm * sample.thickness / uniformRefThickness;
        int ppmUniformInt = static_cast<int>(std::lround(ppmUniform));
        TString uniformFolder = Form("%dppm", ppmUniformInt);

        cout << "=== " << groupLabel << " (uniform ref: " << uniformFolder << ") ===" << endl;

        // このグループで読み込む結果ファイル一覧 (depth違い + 一様含水率リファレンスを末尾に追加)
        vector<TString> resultsPathList;
        vector<TString> labelList;
        for (size_t i = 0; i < depthList.size(); ++i)
        {
            resultsPathList.push_back("../../" + depthList[i].name + "/results.root");
            labelList.push_back("depth=" + depthList[i].depthStr + "m");
        }
        const size_t refFileIdx = resultsPathList.size();
        resultsPathList.push_back("../../../run/" + uniformFolder + "/results.root");
        labelList.push_back(Form("uniform %s (ref)", uniformFolder.Data()));

        vector<TH1F *> vHistcpPosZ, vHistcpPosZ_TNcut, vHistscPosZ;
        vector<TString> legendLabels;
        int refIndex = -1;

        for (size_t fileIdx = 0; fileIdx < resultsPathList.size(); ++fileIdx)
        {
            const Double_t DetectorOffsetZ = 460;   // Detector offset in Z (mm)
            const Double_t irrArea = 600 * 600;     // irradiation surface area (cm^2)
            constexpr Double_t EJ270HalfWidth = 35; // EJ270 width (mm)

            constexpr double scatterEdepLow = 1.0;  // MeV
            constexpr double scatterEdepHigh = 3.0; // MeV
            constexpr double captureEdepLow = 4.5;  // MeV
            constexpr double captureEdepHigh = 5.0; // MeV

            constexpr double TNEnergyCut = 5e-7; // MeV (Thermal neutron cut, 109Cd)

            bool useFidcut = true;
            constexpr double sideCut = 5;                             // mm
            constexpr double fidHalfWidth = EJ270HalfWidth - sideCut; // mm

            bool useCd = true;
            vector<double> cdPlanePosXY = {
                10 + DetectorOffsetZ,
                20 + DetectorOffsetZ,
                40 + DetectorOffsetZ,
                60 + DetectorOffsetZ,
                70 + DetectorOffsetZ,
            }; // XY平面 (法線: Z軸) の位置 [mm]
            vector<double> cdPlanePosZX = {-EJ270HalfWidth + sideCut, EJ270HalfWidth - sideCut}; // ZX平面 (法線: Y軸)
            vector<double> cdPlanePosZY = {-EJ270HalfWidth + sideCut, EJ270HalfWidth - sideCut}; // ZY平面 (法線: X軸)

            TFile *fin = TFile::Open(resultsPathList[fileIdx]);
            if (!fin || fin->IsZombie())
            {
                cerr << "Failed to open " << resultsPathList[fileIdx] << endl;
                if (fin)
                {
                    fin->Close();
                    delete fin;
                }
                continue;
            }

            TTree *HitTree = (TTree *)fin->Get("Hit");
            TTree *RunInfoTree = (TTree *)fin->Get("RunInfo");
            double moonNeutronFlux = 0.0;
            RunInfoTree->SetBranchAddress("TotalFlux", &moonNeutronFlux);
            RunInfoTree->GetEntry(0);
            const Double_t MoonNeutronFlux = moonNeutronFlux;

            TH1F *hPrimEnergy = (TH1F *)fin->Get("hPrimEnergy");
            int nEvents = hPrimEnergy->GetEntries();
            const Double_t eqTime = nEvents / (irrArea * MoonNeutronFlux);

            struct ChamberEventData
            {
                double edepSum = 0.0;
                double cpEdepSum = 0.0;
                double scEdepSum = 0.0;
                double cpPosZ = 0.0;
                double scPosZ = 0.0;
                double cpTriggerTime = numeric_limits<double>::infinity();
                double scTriggerTime = numeric_limits<double>::infinity();
                bool captureflag = false;
                bool scatterflag = false;
            };

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
            };

            map<EventChamberID, ChamberEventData> DetectorchamberMap;

            double primEnergy;
            int eventID;
            char collection[256], fpname[256], fCProc[256];
            double fPreKinE, fEdep, fGTime;
            double fPrePosX, fPrePosY, fPrePosZ, fPostPosX, fPostPosY, fPostPosZ;

            HitTree->SetBranchAddress("eventID", &eventID);
            HitTree->SetBranchAddress("primEnergy", &primEnergy);
            HitTree->SetBranchAddress("collection", collection);
            HitTree->SetBranchAddress("Pname", fpname);
            HitTree->SetBranchAddress("CProc", fCProc);
            HitTree->SetBranchAddress("PreKinE", &fPreKinE);
            HitTree->SetBranchAddress("Edep", &fEdep);
            HitTree->SetBranchAddress("GTime", &fGTime);
            HitTree->SetBranchAddress("PrePosX", &fPrePosX);
            HitTree->SetBranchAddress("PrePosY", &fPrePosY);
            HitTree->SetBranchAddress("PrePosZ", &fPrePosZ);
            HitTree->SetBranchAddress("PostPosX", &fPostPosX);
            HitTree->SetBranchAddress("PostPosY", &fPostPosY);
            HitTree->SetBranchAddress("PostPosZ", &fPostPosZ);

            const vector<string> captureStepParticles = {"alpha", "triton"};
            const vector<string> scatterStepParticles = {"proton"};

            set<int> killedEvents;

            int nHits = HitTree->GetEntries();
            for (int i = 0; i < nHits; ++i)
            {
                HitTree->GetEntry(i);
                if (killedEvents.count(eventID) > 0)
                    continue;

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
                        killedEvents.insert(eventID);
                        continue;
                    }
                }

                if (string(collection) == "TrackerHitsCollection")
                {
                    if (useFidcut && (fabs(fPrePosX) > fidHalfWidth || fabs(fPrePosY) > fidHalfWidth))
                        continue;

                    EventChamberID key{eventID, primEnergy, string(collection)};
                    auto &acc = DetectorchamberMap[key];
                    acc.edepSum += fEdep;

                    const bool isCaptureParticle =
                        find(captureStepParticles.begin(), captureStepParticles.end(), fpname) != captureStepParticles.end();
                    const bool isScatterParticle =
                        find(scatterStepParticles.begin(), scatterStepParticles.end(), fpname) != scatterStepParticles.end();

                    const bool isCaptureDepositStep =
                        (fEdep > 0.0) && isCaptureParticle && (string(fCProc) == "neutronInelastic");
                    const bool isScatterDepositStep =
                        (fEdep > 0.0) && isScatterParticle && (string(fCProc) == "hadElastic");

                    if (isCaptureDepositStep)
                    {
                        acc.cpEdepSum += fEdep;
                        if (!acc.captureflag && acc.cpEdepSum > captureEdepLow && acc.cpEdepSum < captureEdepHigh)
                        {
                            acc.captureflag = true;
                            acc.cpTriggerTime = min(acc.cpTriggerTime, fGTime);
                            acc.cpPosZ = fPrePosZ;
                        }
                    }

                    if (isScatterDepositStep)
                    {
                        acc.scEdepSum += fEdep;
                        if (!acc.scatterflag && acc.scEdepSum > scatterEdepLow && acc.scEdepSum < scatterEdepHigh)
                        {
                            acc.scatterflag = true;
                            acc.scTriggerTime = min(acc.scTriggerTime, fGTime);
                            acc.scPosZ = fPrePosZ;
                        }
                    }
                }
            }

            cout << resultsPathList[fileIdx] << ": Virtual Cd layer absorbed " << killedEvents.size() << " / " << nEvents << " events" << endl;

            int BinWidthZ = 5; // mm
            int minZ = 0;      // mm
            int maxZ = 80;     // mm
            int nBinsZ = (maxZ - minZ) / BinWidthZ;

            TH1F *h1_cpposZ = new TH1F(Form("h1_cpposZ_%s_%zu", groupLabel.Data(), fileIdx), Form("Capture Position Z Distribution;Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
            TH1F *h1_cpposZ_TNcut = new TH1F(Form("h1_cpposZ_TNcut_%s_%zu", groupLabel.Data(), fileIdx), Form("Capture Position Z Distribution (TN cut / %.1f eV <);Z (mm);Counts (s^{-1} %d mm^{-1})", TNEnergyCut * 1e6, BinWidthZ), nBinsZ, minZ, maxZ);
            TH1F *h1_scposZ = new TH1F(Form("h1_scposZ_%s_%zu", groupLabel.Data(), fileIdx), Form("Scatter Position Z Distribution;Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
            h1_cpposZ->SetDirectory(nullptr);
            h1_cpposZ_TNcut->SetDirectory(nullptr);
            h1_scposZ->SetDirectory(nullptr);

            for (const auto &entry : DetectorchamberMap)
            {
                const auto &id = entry.first;
                const auto &data = entry.second;

                if (data.captureflag)
                {
                    double capturePosZ = data.cpPosZ - DetectorOffsetZ;
                    h1_cpposZ->Fill(capturePosZ);
                    if (id.primEnergy > TNEnergyCut)
                        h1_cpposZ_TNcut->Fill(capturePosZ);
                }
                if (data.scatterflag)
                {
                    double scatterPosZ = data.scPosZ - DetectorOffsetZ;
                    h1_scposZ->Fill(scatterPosZ);
                }
            }

            h1_cpposZ->Scale(1.0 / eqTime);
            h1_cpposZ_TNcut->Scale(1.0 / eqTime);
            h1_scposZ->Scale(1.0 / eqTime);

            vHistcpPosZ.push_back(h1_cpposZ);
            vHistcpPosZ_TNcut.push_back(h1_cpposZ_TNcut);
            vHistscPosZ.push_back(h1_scposZ);
            legendLabels.push_back(labelList[fileIdx]);
            if (fileIdx == refFileIdx)
                refIndex = (int)vHistcpPosZ.size() - 1;

            fin->Close();
            delete fin;
        }

        const bool hasRef = (refIndex >= 0);
        if (!hasRef)
            cerr << "Warning: uniform reference folder not found: ../../../run/" << uniformFolder << "/results.root" << endl;

        const int nCurves = (int)vHistcpPosZ.size();
        if (nCurves == 0)
            continue;
        int nColors = gStyle->GetNumberOfColors();
        vector<Color_t> colorPalette;
        for (int i = 0; i < nCurves; ++i)
            colorPalette.push_back(gStyle->GetColorPalette(i * nColors / nCurves));

        vector<TCanvas *> vCan;

        // Capture Position Z
        {
            TCanvas *c = new TCanvas("cCpposZ_" + groupLabel, "Capture Position Z (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.6, 0.88, 0.88);
            leg->SetNColumns(2);
            for (int i = 0; i < nCurves; ++i)
            {
                bool isRef = (hasRef && i == refIndex);
                vHistcpPosZ[i]->SetLineColor(isRef ? kBlack : colorPalette[i]);
                vHistcpPosZ[i]->SetLineStyle(isRef ? 2 : 1);
                vHistcpPosZ[i]->SetLineWidth(2);
                vHistcpPosZ[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
                leg->AddEntry(vHistcpPosZ[i], legendLabels[i], "l");
            }
            leg->Draw();
            vCan.push_back(c);
        }

        // Capture Position Z (Ratio to uniform reference)
        if (hasRef)
        {
            TCanvas *c = new TCanvas("cCpposZ_ratio_" + groupLabel, "Capture Position Z Ratio (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.12, 0.88, 0.35);
            leg->SetNColumns(2);
            for (int i = 0; i < nCurves; ++i)
            {
                TH1F *hRatio = (TH1F *)vHistcpPosZ[i]->Clone(Form("h1_cpposZ_ratio_%s_%d", groupLabel.Data(), i));
                hRatio->Divide(vHistcpPosZ[refIndex]);
                hRatio->SetTitle(";Z (mm);Ratio to uniform ref");
                hRatio->SetMinimum(0);
                hRatio->SetMaximum(2.0);
                hRatio->Draw(i == 0 ? "HIST E" : "HIST E SAME");
                leg->AddEntry(hRatio, legendLabels[i], "l");
            }
            leg->Draw();
            vCan.push_back(c);
        }

        // Capture Position Z (TN cut)
        {
            TCanvas *c = new TCanvas("cCpposZ_TNcut_" + groupLabel, "Capture Position Z TN cut (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.6, 0.88, 0.88);
            leg->SetNColumns(2);
            for (int i = 0; i < nCurves; ++i)
            {
                bool isRef = (hasRef && i == refIndex);
                vHistcpPosZ_TNcut[i]->SetLineColor(isRef ? kBlack : colorPalette[i]);
                vHistcpPosZ_TNcut[i]->SetLineStyle(isRef ? 2 : 1);
                vHistcpPosZ_TNcut[i]->SetLineWidth(2);
                vHistcpPosZ_TNcut[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
                leg->AddEntry(vHistcpPosZ_TNcut[i], legendLabels[i], "l");
            }
            leg->Draw();
            vCan.push_back(c);
        }

        // Capture Position Z TN cut (Ratio to uniform reference)
        if (hasRef)
        {
            TCanvas *c = new TCanvas("cCpposZ_TNcut_ratio_" + groupLabel, "Capture Position Z TN cut Ratio (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.12, 0.88, 0.35);
            leg->SetNColumns(2);
            for (int i = 0; i < nCurves; ++i)
            {
                TH1F *hRatio = (TH1F *)vHistcpPosZ_TNcut[i]->Clone(Form("h1_cpposZ_TNcut_ratio_%s_%d", groupLabel.Data(), i));
                hRatio->Divide(vHistcpPosZ_TNcut[refIndex]);
                hRatio->SetTitle(";Z (mm);Ratio to uniform ref");
                hRatio->SetMinimum(0);
                hRatio->SetMaximum(2.0);
                hRatio->Draw(i == 0 ? "HIST E" : "HIST E SAME");
                leg->AddEntry(hRatio, legendLabels[i], "l");
            }
            leg->Draw();
            vCan.push_back(c);
        }

        // Scatter Position Z
        {
            TCanvas *c = new TCanvas("cScposZ_" + groupLabel, "Scatter Position Z (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.6, 0.88, 0.88);
            leg->SetNColumns(2);
            for (int i = 0; i < nCurves; ++i)
            {
                bool isRef = (hasRef && i == refIndex);
                vHistscPosZ[i]->SetLineColor(isRef ? kBlack : colorPalette[i]);
                vHistscPosZ[i]->SetLineStyle(isRef ? 2 : 1);
                vHistscPosZ[i]->SetLineWidth(2);
                vHistscPosZ[i]->Draw(i == 0 ? "HIST E" : "HIST E SAME");
                leg->AddEntry(vHistscPosZ[i], legendLabels[i], "l");
            }
            leg->Draw();
            vCan.push_back(c);
        }

        // Scatter Position Z (Ratio to uniform reference)
        if (hasRef)
        {
            TCanvas *c = new TCanvas("cScposZ_ratio_" + groupLabel, "Scatter Position Z Ratio (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.12, 0.88, 0.35);
            leg->SetNColumns(2);
            for (int i = 0; i < nCurves; ++i)
            {
                TH1F *hRatio = (TH1F *)vHistscPosZ[i]->Clone(Form("h1_scposZ_ratio_%s_%d", groupLabel.Data(), i));
                hRatio->Divide(vHistscPosZ[refIndex]);
                hRatio->SetTitle(";Z (mm);Ratio to uniform ref");
                hRatio->SetMinimum(0);
                hRatio->SetMaximum(2.0);
                hRatio->Draw(i == 0 ? "HIST E" : "HIST E SAME");
                leg->AddEntry(hRatio, legendLabels[i], "l");
            }
            leg->Draw();
            vCan.push_back(c);
        }

        // save PDF (by_H_total の各グループフォルダの中に保存)
        TString outDir = "../fig/by_H_total/" + groupLabel;
        gSystem->mkdir(outDir, kTRUE);
        TString fPdfOut = outDir + "/" + groupLabel + "_Summary.pdf";
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

        for (auto *c : vCan)
            delete c;
        vCan.clear();
    }

    return;
}
