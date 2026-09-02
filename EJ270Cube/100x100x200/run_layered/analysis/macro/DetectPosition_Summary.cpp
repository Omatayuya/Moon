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
#include <TSystem.h>

struct PosHistSet
{
    TH1F *cpPosZ = nullptr;
    TH1F *cpPosZ_TNcut = nullptr;
    TH1F *scPosZ = nullptr;
};

// results.root 1つを読み込み、Capture/Scatter Position Z のヒストグラムを作る
PosHistSet BuildPosHist(const TString &resultsRootPath, const TString &label)
{
    PosHistSet hs;

    const Double_t DetectorOffsetZ = 400;  // Detector offset in Z (mm)
    const Double_t irrArea = 600 * 600;    // irradiation surface area (cm^2)
    constexpr double scatterEdepLow = 1.0; // MeV
    constexpr double captureEdepLow = 4.5; // MeV
    constexpr double TNEnergyCut = 5e-7;   // MeV (Thermal neutron cut, 109Cd)

    TFile *fin = TFile::Open(resultsRootPath);
    if (!fin || fin->IsZombie())
    {
        cerr << "Failed to open " << resultsRootPath << endl;
        return hs;
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
    char chamberNb[256], fpname[256], fCProc[256];
    double fEdep, fGTime;
    double fPrePosX, fPrePosY, fPrePosZ;

    HitTree->SetBranchAddress("eventID", &eventID);
    HitTree->SetBranchAddress("primEnergy", &primEnergy);
    HitTree->SetBranchAddress("chamberNb", chamberNb);
    HitTree->SetBranchAddress("Pname", fpname);
    HitTree->SetBranchAddress("CProc", fCProc);
    HitTree->SetBranchAddress("Edep", &fEdep);
    HitTree->SetBranchAddress("GTime", &fGTime);
    HitTree->SetBranchAddress("PrePosX", &fPrePosX);
    HitTree->SetBranchAddress("PrePosY", &fPrePosY);
    HitTree->SetBranchAddress("PrePosZ", &fPrePosZ);

    const vector<string> captureStepParticles = {"alpha", "triton"};
    const vector<string> scatterStepParticles = {"proton"};

    int nHits = HitTree->GetEntries();
    for (int i = 0; i < nHits; ++i)
    {
        HitTree->GetEntry(i);
        if (string(chamberNb) != "target1")
            continue;

        EventChamberID key{eventID, primEnergy, string(chamberNb)};
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
            if (!acc.captureflag && acc.cpEdepSum > captureEdepLow)
            {
                acc.captureflag = true;
                acc.cpTriggerTime = min(acc.cpTriggerTime, fGTime);
                acc.cpPosZ = fPrePosZ;
            }
        }

        if (isScatterDepositStep)
        {
            acc.scEdepSum += fEdep;
            if (!acc.scatterflag && acc.scEdepSum > scatterEdepLow)
            {
                acc.scatterflag = true;
                acc.scTriggerTime = min(acc.scTriggerTime, fGTime);
                acc.scPosZ = fPrePosZ;
            }
        }
    }

    int BinWidthZ = 5; // mm
    int minZ = 0;      // mm
    int maxZ = 200;    // mm
    int nBinsZ = (maxZ - minZ) / BinWidthZ;

    hs.cpPosZ = new TH1F("h1_cpposZ_" + label, Form("Capture Position Z Distribution;Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
    hs.cpPosZ_TNcut = new TH1F("h1_cpposZ_TNcut_" + label, Form("Capture Position Z Distribution (TN cut / %.1f eV <);Z (mm);Counts (s^{-1} %d mm^{-1})", TNEnergyCut * 1e6, BinWidthZ), nBinsZ, minZ, maxZ);
    hs.scPosZ = new TH1F("h1_scposZ_" + label, Form("Scatter Position Z Distribution;Z (mm);Counts (s^{-1} %d mm^{-1})", BinWidthZ), nBinsZ, minZ, maxZ);
    hs.cpPosZ->SetDirectory(nullptr);
    hs.cpPosZ_TNcut->SetDirectory(nullptr);
    hs.scPosZ->SetDirectory(nullptr);

    for (const auto &entry : DetectorchamberMap)
    {
        const auto &id = entry.first;
        const auto &data = entry.second;

        if (data.captureflag)
        {
            double capturePosZ = data.cpPosZ - DetectorOffsetZ;
            hs.cpPosZ->Fill(capturePosZ);
            if (id.primEnergy > TNEnergyCut)
                hs.cpPosZ_TNcut->Fill(capturePosZ);
        }
        if (data.scatterflag)
        {
            double scatterPosZ = data.scPosZ - DetectorOffsetZ;
            hs.scPosZ->Fill(scatterPosZ);
        }
    }

    hs.cpPosZ->Scale(1.0 / eqTime);
    hs.cpPosZ_TNcut->Scale(1.0 / eqTime);
    hs.scPosZ->Scale(1.0 / eqTime);

    fin->Close();
    delete fin;

    return hs;
}

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

void DetectPosition_Summary()
{
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

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

    // folders.list を読み込み、box_10x10x{thickness}m_depth_{depth}m_H_{ppm}ppm を解析
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

    std::regex folderNameRe("box_10x10x([0-9.]+)m_depth_([0-9.]+)m_H_([0-9.]+)ppm");

    // (ppmStr, thicknessStr) ごとに depth 違いのフォルダをまとめる = by_H_total の各フォルダに対応
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

        // 一様含水率の参照ppmを計算: ppm_layer x thickness_layer = ppm_uniform x 1.5
        double ppmUniform = sample.ppm * sample.thickness / uniformRefThickness;
        int ppmUniformInt = static_cast<int>(std::lround(ppmUniform));
        TString uniformFolder = Form("%dppm", ppmUniformInt);

        cout << "=== " << groupLabel << " (uniform ref: " << uniformFolder << ") ===" << endl;

        vector<TH1F *> vHistcpPosZ, vHistcpPosZ_TNcut, vHistscPosZ;
        vector<TString> legendLabels;

        for (size_t i = 0; i < depthList.size(); ++i)
        {
            const FolderInfo &fi = depthList[i];
            TString resultsPath = "../../" + fi.name + "/results.root";
            TString label = Form("%s_d%zu", groupLabel.Data(), i);
            PosHistSet hs = BuildPosHist(resultsPath, label);
            if (!hs.cpPosZ)
                continue;
            vHistcpPosZ.push_back(hs.cpPosZ);
            vHistcpPosZ_TNcut.push_back(hs.cpPosZ_TNcut);
            vHistscPosZ.push_back(hs.scPosZ);
            legendLabels.push_back("depth=" + fi.depthStr + "m");
        }

        // 一様含水率リファレンス
        TString uniformResultsPath = "../../../run/" + uniformFolder + "/results.root";
        PosHistSet hsRef = BuildPosHist(uniformResultsPath, groupLabel + "_ref");
        bool hasRef = (hsRef.cpPosZ != nullptr);
        int refIndex = -1;
        if (hasRef)
        {
            refIndex = (int)vHistcpPosZ.size();
            vHistcpPosZ.push_back(hsRef.cpPosZ);
            vHistcpPosZ_TNcut.push_back(hsRef.cpPosZ_TNcut);
            vHistscPosZ.push_back(hsRef.scPosZ);
            legendLabels.push_back(Form("uniform %s (ref)", uniformFolder.Data()));
        }
        else
        {
            cerr << "Warning: uniform reference folder not found: " << uniformResultsPath << endl;
        }

        const int nCurves = (int)vHistcpPosZ.size();
        int nColors = gStyle->GetNumberOfColors();
        vector<Color_t> colorPalette;
        for (int i = 0; i < nCurves; ++i)
            colorPalette.push_back(gStyle->GetColorPalette(i * nColors / nCurves));

        auto styleAndDraw = [&](vector<TH1F *> &hists, const char *drawOpt)
        {
            for (int i = 0; i < nCurves; ++i)
            {
                bool isRef = (i == refIndex);
                hists[i]->SetLineColor(isRef ? kBlack : colorPalette[i]);
                hists[i]->SetLineStyle(isRef ? 2 : 1);
                hists[i]->SetLineWidth(2);
                hists[i]->Draw(i == 0 ? drawOpt : (TString(drawOpt) + " SAME").Data());
            }
        };

        auto drawRatioCanvas = [&](vector<TH1F *> &hists, const TString &canvasName, const TString &canvasTitle, const TString &yAxisTitle) -> TCanvas *
        {
            if (!hasRef)
                return nullptr;
            TCanvas *c = new TCanvas(canvasName, canvasTitle, 800, 600);
            TLegend *leg = new TLegend(0.55, 0.12, 0.88, 0.35);
            leg->SetNColumns(2);
            for (int i = 0; i < nCurves; ++i)
            {
                TH1F *hRatio = (TH1F *)hists[i]->Clone(canvasName + Form("_ratio_%d", i));
                hRatio->Divide(hists[refIndex]);
                hRatio->SetTitle(";Z (mm);" + yAxisTitle);
                hRatio->SetMinimum(0);
                hRatio->SetMaximum(2.0);
                hRatio->Draw(i == 0 ? "HIST E" : "HIST E SAME");
                leg->AddEntry(hRatio, legendLabels[i], "l");
            }
            leg->Draw();
            return c;
        };

        vector<TCanvas *> vCan;

        // Capture Position Z
        {
            TCanvas *c = new TCanvas("cCpposZ_" + groupLabel, "Capture Position Z (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.6, 0.88, 0.88);
            leg->SetNColumns(2);
            styleAndDraw(vHistcpPosZ, "HIST E");
            for (int i = 0; i < nCurves; ++i)
                leg->AddEntry(vHistcpPosZ[i], legendLabels[i], "l");
            leg->Draw();
            vCan.push_back(c);
        }
        // Capture Position Z (Ratio to uniform reference)
        if (TCanvas *c = drawRatioCanvas(vHistcpPosZ, "cCpposZ_ratio_" + groupLabel, "Capture Position Z Ratio (" + groupLabel + ")", "Ratio to uniform ref"))
            vCan.push_back(c);

        // Capture Position Z (TN cut)
        {
            TCanvas *c = new TCanvas("cCpposZ_TNcut_" + groupLabel, "Capture Position Z TN cut (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.6, 0.88, 0.88);
            leg->SetNColumns(2);
            styleAndDraw(vHistcpPosZ_TNcut, "HIST E");
            for (int i = 0; i < nCurves; ++i)
                leg->AddEntry(vHistcpPosZ_TNcut[i], legendLabels[i], "l");
            leg->Draw();
            vCan.push_back(c);
        }
        // Capture Position Z TN cut (Ratio to uniform reference)
        if (TCanvas *c = drawRatioCanvas(vHistcpPosZ_TNcut, "cCpposZ_TNcut_ratio_" + groupLabel, "Capture Position Z TN cut Ratio (" + groupLabel + ")", "Ratio to uniform ref"))
            vCan.push_back(c);

        // Scatter Position Z
        {
            TCanvas *c = new TCanvas("cScposZ_" + groupLabel, "Scatter Position Z (" + groupLabel + ")", 800, 600);
            TLegend *leg = new TLegend(0.55, 0.6, 0.88, 0.88);
            leg->SetNColumns(2);
            styleAndDraw(vHistscPosZ, "HIST E");
            for (int i = 0; i < nCurves; ++i)
                leg->AddEntry(vHistscPosZ[i], legendLabels[i], "l");
            leg->Draw();
            vCan.push_back(c);
        }
        // Scatter Position Z (Ratio to uniform reference)
        if (TCanvas *c = drawRatioCanvas(vHistscPosZ, "cScposZ_ratio_" + groupLabel, "Scatter Position Z Ratio (" + groupLabel + ")", "Ratio to uniform ref"))
            vCan.push_back(c);

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
