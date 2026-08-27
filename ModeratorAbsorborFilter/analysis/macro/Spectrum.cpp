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

namespace
{
    constexpr Double_t neutronYield = 0.116; // Cf-252 neutron yield (s^-1 Bq^-1)

    const int nBins = 500;
    const double eMin = 1e-10;
    const double eMax = 1e2;
    double log_width = (log(eMax) - log(eMin)) / nBins;

    double gNormTh = 1.0;
    double gNormEp = 1.0;
    double gNormFa = 1.0;

    constexpr double kThermalUpper = 1e-7;   // MeV
    constexpr double kEpithermalUpper = 0.1; // MeV
    constexpr double kTwoOverSqrtPi = 1.1283791670955126;

    enum class FastModel
    {
        Fission,
        Evaporation,
        Gaussian
    };

    struct SpectrumShape
    {
        double T0{2.53e-8}; // Thermal parameter (MeV)

        double Ed{7.07e-8}; // Epithermal cutoff (MeV)
        double beta{1};     // Epithermal cutoff (MeV)
        double b{0.03};     // Epithermal slope deviation from 1/E

        double alpha{0.5}; // Fission fast exponent
        double beta2{1};   // Fission fast temperature (MeV)

        double Te{1.0}; // Evaporation temperature (MeV)

        double Em{2.0};    // Gaussian mean (MeV)
        double sigma{0.5}; // Gaussian width (MeV)

        FastModel fastModel{FastModel::Fission};
    };

    SpectrumShape gShape{};

    double PhiThermal(double E, double T0)
    {
        return E / pow(T0, 2) * exp(-E / T0);
    }

    double PhiEpithermal(double E, double Ed, double beta, double b)
    {
        const double lowCut = 1.0 - exp(-pow(E / Ed, 2));
        const double shape = pow(E, b - 1);
        const double highCut = exp(-E / beta);
        return lowCut * shape * highCut;
    }

    double PhiFast(double E, const SpectrumShape &shape)
    {
        switch (shape.fastModel)
        {
        case FastModel::Fission:
        {
            const double pref = std::pow(E, shape.alpha);
            return pref * std::exp(-E / shape.beta2);
        }
        case FastModel::Evaporation:
        {
            const double x = E / shape.Te;
            return x * x * std::exp(-x);
        }
        case FastModel::Gaussian:
        default:
        {
            const double invSigma = 1.0 / (shape.sigma * std::sqrt(2.0 * M_PI));
            const double d = E - shape.Em;
            return invSigma * std::exp(-0.5 * (d * d) / (shape.sigma * shape.sigma));
        }
        }
    }

    void UpdateNormalizations()
    {
        TF1 fThLog("fThLog", [](double *x, double *)
                   {
        double u = x[0];
        double E = exp(u);
        return PhiThermal(E, gShape.T0) * E; }, log(eMin), log(eMax), 0);

        gNormTh = fThLog.Integral(log(eMin), log(eMax));

        TF1 fEpLog("fEpLog", [](double *x, double *)
                   {
        double u = x[0];
        double E = exp(u);
        return PhiEpithermal(E, gShape.Ed, gShape.beta, gShape.b) * E; }, log(eMin), log(eMax), 0);

        gNormEp = fEpLog.Integral(log(eMin), log(eMax));

        TF1 fFa("fFa", [](double *x, double *)
                { return PhiFast(x[0], gShape); }, eMin, eMax, 0);
        gNormFa = fFa.Integral(eMin, eMax);
    }

    double PhiThermal_normalized(double E)
    {
        return PhiThermal(E, gShape.T0) / gNormTh;
    }
    double PhiEpithermal_normalized(double E)
    {
        return PhiEpithermal(E, gShape.Ed, gShape.beta, gShape.b) / gNormEp;
    }
    double PhiFast_normalized(double E)
    {
        return PhiFast(E, gShape) / gNormFa;
    }

    // 真のスペクトルモデル (図の式: phi = Pth*phi_th + Pe*phi_e + Pf*phi_f)
    // params[0]: P_th, params[1]: P_e, params[2]: P_f
    double ModelFunction(double *x, const double *params)
    {
        // UpdateNormalizations();
        const double energy = x[0];
        const double norm = params[0];
        const double Pth = params[1];
        const double Pe = params[2];
        // const double Pf = params[3];
        const double Pf = 1.0 - Pth - Pe;

        const double th = (gNormTh > 0.0) ? PhiThermal_normalized(energy) : 0.0;
        const double ep = (gNormEp > 0.0) ? PhiEpithermal_normalized(energy) : 0.0;
        const double fa = (gNormFa > 0.0) ? PhiFast_normalized(energy) : 0.0;

        return norm * (Pth * th + Pe * ep + Pf * fa);
    }

}

void Spectrum()
{
    vector<string> folder = {"build"};

    for (int folderID = 0; folderID < folder.size(); folderID++)
    {

        // ==================================================================

        const string TargetParticle = "neutron";
        constexpr Double_t neutronYield = 5.72086; // Moon neutron yield (cm^-2 s^-1)
        constexpr double S = 2500; // Input area (cm^2)

        double totalFlux = neutronYield * S; // Total neutron flux (s^-1)
        double Todsdt = 1/totalFlux; 

        double xmin = 5e-10;
        double xmax = 1e4;
        const int xbins = 500;

        double xmin_na = 1e5;
        double xmax_na = 1e7;
        const int xbins_na = 100;

        bool isCaptureEdepWindow = true;
        double na_capture_lower = 4.78; // MeV
        double na_capture_upper = 10;   // MeV

        double Threshold = 0.001; // MeV

        auto AcceptDirection = [](const string &chamber, double dx, double dy, double dz)
        {
            if (chamber == "target1001")
                return dx > 0.0;
            if (chamber == "target1002")
                return dz < 0.0;
            if (chamber == "target1003")
                return dz > 0.0;
            if (chamber == "target1004")
                return dy > 0.0;
            if (chamber == "target1005")
                return dy < 0.0;
            if (chamber == "target1006")
                return dx < 0.0;
            return false;
        };

        for (auto axis : {"X", "Y", "Z"})
        {
            gStyle->SetLabelFont(62, axis);
            gStyle->SetTitleFont(62, axis);
            gStyle->SetTitleOffset(1.2, axis);  // 軸タイトルのオフセット
            gStyle->SetLabelSize(0.04, axis); // 目盛り数字のサイズ
            gStyle->SetTitleSize(0.04, axis); // 軸タイトルのサイズ
        }
        gStyle->SetTextFont(62);
        gStyle->SetTitleFont(62, "");

        gStyle->SetOptStat(0);
        // ==================================================================

        string filepath_ini = "../../" + folder[folderID] + "/initial_particle.csv";
        string filepath_op = "../../" + folder[folderID] + "/output.csv";
        ifstream ifs_op(filepath_op.c_str());
        ifstream ifs_ini(filepath_ini.c_str());
        string str;

        if (ifs_op.fail() || ifs_ini.fail())
        {
            cout << "Failed to open file." << endl;
            if (ifs_op.fail())
            {
                cout << "Output file: " << filepath_op << " not found." << endl;
            }
            if (ifs_ini.fail())
            {
                cout << "Initial particle file: " << filepath_ini << " not found." << endl;
            }
            return;
        }

        double *bins = new double[xbins + 1];
        for (int i = 0; i <= xbins; ++i)
        {
            bins[i] = xmin * pow(xmax / xmin, static_cast<double>(i) / xbins);
        }

        /*Input ParticleData*/
        struct InputEventData
        {
            int eventID = 0;
            double initialEnergy = 0.0;
            double PosX = 0.0;
            double PosY = 0.0;
            double PosZ = 0.0;
        };

        InputEventData InputData;
        TTree *InputDataTree = new TTree("InputDataTree", "event-level tree from csvfile");
        InputDataTree->Branch("eventID", &InputData.eventID, "eventID/I");
        InputDataTree->Branch("initialEnergy", &InputData.initialEnergy, "initialEnergy/D");
        InputDataTree->Branch("PosX", &InputData.PosX, "PosX/D");
        InputDataTree->Branch("PosY", &InputData.PosY, "PosY/D");
        InputDataTree->Branch("PosZ", &InputData.PosZ, "PosZ/D");

        TH1F *h1_ip = new TH1F("h1_ip", "Moon Neutron Spectrum; Energy (MeV); Energy #times Intensity", xbins, bins);

        // 全入射粒子の球面位置分布（分母）: phi vs cos(theta) 等面積射影
        int phi_bins = 10;       // 5-degree bins
        int cos_theta_bins = 18; // 0.05 bins in cos(theta)
        TH2F *h2_origin_all = new TH2F("h2_origin_all",
                                       "Input Particle Distribution;#phi (rad);cos(#theta)",
                                       phi_bins, -M_PI, M_PI,
                                       cos_theta_bins, -1.0, 1.0);

        while (getline(ifs_ini, str))
        {
            stringstream ss(str);
            string item;

            getline(ss, item, ',');
            InputData.eventID = stoi(item);
            getline(ss, item, ',');
            InputData.initialEnergy = stod(item);
            getline(ss, item, ',');
            InputData.PosX = stod(item);
            getline(ss, item, ',');
            InputData.PosY = stod(item);
            getline(ss, item, ',');
            InputData.PosZ = stod(item);
            InputDataTree->Fill();

            h1_ip->Fill(InputData.initialEnergy);

            double r = sqrt(InputData.PosX * InputData.PosX +
                            InputData.PosY * InputData.PosY +
                            InputData.PosZ * InputData.PosZ);
            if (r > 0)
            {
                double phi = atan2(InputData.PosY, InputData.PosX);
                double cosTheta = InputData.PosZ / r;
                h2_origin_all->Fill(phi, cosTheta);
            }
        }

        InputDataTree->BuildIndex("eventID");
        int nEvents = InputDataTree->GetEntries();

        /* Output Data*/

        // struct ChamberEventData
        // {
        //     int eventID = 0;
        //     string chamberNb;
        //     double edepSum = 0.0;
        //     double captureEdepSum = 0.0;
        //     double scatterEdepSum = 0.0;
        //     double TriggerTime = numeric_limits<double>::infinity();
        // };

        // ChamberEventData branchData;
        // TTree *DetectorchamberTree = new TTree("DetectorchamberTree", "event-chamber-level tree from csvfile");
        // DetectorchamberTree->Branch("eventID", &branchData.eventID, "eventID/I");
        // DetectorchamberTree->Branch("chamberNb", &branchData.chamberNb);
        // DetectorchamberTree->Branch("edepSum", &branchData.edepSum, "edepSum/D");
        // DetectorchamberTree->Branch("captureEdepSum", &branchData.captureEdepSum, "captureEdepSum/D");
        // DetectorchamberTree->Branch("scatterEdepSum", &branchData.scatterEdepSum, "scatterEdepSum/D");
        // DetectorchamberTree->Branch("TriggerTime", &branchData.TriggerTime, "TriggerTime/D");

        // // (eventID, chamberNb) をキーに集計するmap
        // map<pair<int, string>, ChamberEventData> chamberAccumulator;

        // const vector<string> captureStepParticles = {"alpha", "triton"};
        // // const vector<string> scatterStepParticles = {"proton", "C12", "O16", "N14"};
        // const vector<string> scatterStepParticles = {"proton"};
        // constexpr int primaryTrackID = 1;

        // while (getline(ifs_op, str))
        // {
        //     stringstream ss(str);
        //     string item;
        //     getline(ss, item, ',');
        //     int eventID = stoi(item);
        //     getline(ss, item, ',');
        //     string chamberNb = item;
        //     getline(ss, item, ',');
        //     int fPID = stoi(item);
        //     getline(ss, item, ',');
        //     int fPPID = stoi(item);
        //     getline(ss, item, ',');
        //     string fpname = item;
        //     getline(ss, item, ',');
        //     string fPreProc = item;
        //     getline(ss, item, ',');
        //     string fPostProc = item;
        //     getline(ss, item, ',');
        //     string fCProc = item;
        //     getline(ss, item, ',');
        //     string fCModel = item;
        //     getline(ss, item, ',');
        //     double fPreKinE = stod(item);
        //     getline(ss, item, ',');
        //     double fPostKinE = stod(item);
        //     getline(ss, item, ',');
        //     double fEdep = stod(item);
        //     getline(ss, item, ',');
        //     double fGTime = stod(item);
        //     getline(ss, item, ',');
        //     double fLTime = stod(item);
        //     getline(ss, item, ',');
        //     double fPTime = stod(item);
        //     getline(ss, item, ',');
        //     double fDTime = stod(item);
        //     getline(ss, item, ',');
        //     double fPrePosx = stod(item);
        //     getline(ss, item, ',');
        //     double fPrePosy = stod(item);
        //     getline(ss, item, ',');
        //     double fPrePosz = stod(item);
        //     getline(ss, item, ',');
        //     double fPosPosx = stod(item);
        //     getline(ss, item, ',');
        //     double fPosPosy = stod(item);
        //     getline(ss, item, ',');
        //     double fPosPosz = stod(item);
        //     getline(ss, item, ',');
        //     double fStepLength = stod(item);

        //     double fEdepQ = 0.0;     // Initialize quenched energy deposit
        //     const double kB = 0.012; // Birks' constant (mm/MeV)
        //     const double S = 1.0;    // Scintillation efficiency

        //     if (fStepLength > 0. && fEdep > 0.)
        //     {
        //         double dedx = fEdep / fStepLength; // [MeV/mm]
        //         double quenchingFactor = 1. / (1. + kB * dedx);
        //         fEdepQ = S * fEdep * quenchingFactor;
        //     }

        //     const double preKinEMeV = fPreKinE / 1e6; // Convert eV to MeV
        //     const double edepMeV = fEdep / 1e6;      // Convert eV to MeV

        //     if (chamberNb.rfind("target2", 0) == 0)
        //     {
        //         auto &acc = chamberAccumulator[{eventID, chamberNb}];
        //         acc.eventID = eventID;
        //         acc.chamberNb = chamberNb;
        //         acc.edepSum += edepMeV;

        //         /* screening conditions for Capture & Scatter event*/
        //         const bool isCaptureParticle =
        //             find(captureStepParticles.begin(), captureStepParticles.end(), fpname) != captureStepParticles.end();
        //         const bool isScatterParticle =
        //             find(scatterStepParticles.begin(), scatterStepParticles.end(), fpname) != scatterStepParticles.end();

        //         const bool isPrimaryNeutronSecondary = (fPPID == primaryTrackID);

        //         const bool isCaptureDepositStep =
        //             (edepMeV > 0.0) &&
        //             isCaptureParticle &&
        //             (fCProc == "neutronInelastic");
        //         const bool isScatterDepositStep =
        //             (edepMeV > 0.0) &&
        //             isScatterParticle &&
        //             isPrimaryNeutronSecondary &&
        //             (fCProc == "hadElastic");

        //         if (isCaptureDepositStep)
        //         {
        //             acc.captureEdepSum += edepMeV;
        //             acc.TriggerTime = min(acc.TriggerTime, fGTime);
        //         }

        //         if (isScatterDepositStep)
        //         {
        //             acc.scatterEdepSum += edepMeV;
        //             acc.TriggerTime = min(acc.TriggerTime, fGTime);
        //         }
        //     }
        // }
        // ifs_op.close();

        // // (eventID, chamberNb) ごとに1エントリ Fill
        // for (const auto &[key, data] : chamberAccumulator)
        // {
        //     branchData = data;
        //     DetectorchamberTree->Fill();
        //     // cout << "  chamberNb=" << branchData.chamberNb
        //     //      << "  scatterEdep=" << branchData.scatterEdepSum
        //     //      << "  captureEdep=" << branchData.captureEdepSum << endl;
        // }
        // const Long64_t DetectorchamberEntries = DetectorchamberTree->GetEntries();

        // /* Analysis */
        // // Clear accumulators for potential reuse
        // chamberAccumulator.clear();

        // const string referenceCaptureChamber = "target2001";

        // struct Eventinfo
        // {
        //     double initialEnergy;
        //     double TriggerTime;
        // };
        // map<int, Eventinfo> hasScatterEventmap;
        // map<int, Eventinfo> hasEdepEventmap;

        // for (Long64_t entry = 0; entry < DetectorchamberEntries; ++entry)
        // {
        //     DetectorchamberTree->GetEntry(entry);
        //     const int eventID = branchData.eventID;

        //     InputDataTree->GetEntryWithIndex(eventID);
        //     double initialEnergy = InputData.initialEnergy;

        //     const bool hasScatter =
        //         branchData.chamberNb == referenceCaptureChamber &&
        //         branchData.scatterEdepSum > Threshold;

        //     const bool hasEdep =
        //         branchData.chamberNb != referenceCaptureChamber &&
        //         (branchData.scatterEdepSum > Threshold || branchData.captureEdepSum > 0);

        //     if (hasScatter)
        //     {
        //         hasScatterEventmap.emplace(eventID, Eventinfo{initialEnergy, branchData.TriggerTime});
        //     }
        //     if (hasEdep)
        //     {
        //         hasEdepEventmap.emplace(eventID, Eventinfo{initialEnergy, branchData.TriggerTime});
        //     }
        // }

        // cout << "[Debug] hasScatterEventmap size: " << hasScatterEventmap.size() << endl;
        // cout << "[Debug] hasEdepEventmap size:   " << hasEdepEventmap.size() << endl;

        // set<int> ScatterAndEdepEventSet;
        // for (const auto &scatterEvent : hasScatterEventmap)
        // {
        //     const int &eventID = scatterEvent.first;
        //     const double &scatterTriggerTime = scatterEvent.second.TriggerTime;

        //     auto edepIt = hasEdepEventmap.find(eventID);
        //     if (edepIt == hasEdepEventmap.end())
        //         continue; // scatter のみでedepなし

        //     const double &edepTriggerTime = edepIt->second.TriggerTime;

        //     if (scatterTriggerTime < edepTriggerTime)
        //     {
        //         ScatterAndEdepEventSet.insert(eventID);
        //     }
        // }
        // cout << "Number of events with both scatter and edep: " << ScatterAndEdepEventSet.size() << endl;

        // // 選別イベントの球面位置分布
        // TH2F *h2_origin_selected = new TH2F("h2_origin_selected",
        //                                     "Detected Events Distribution;#phi (rad);cos(#theta)",
        //                                     phi_bins, -M_PI, M_PI,
        //                                     cos_theta_bins, -1.0, 1.0);

        // for (const int evID : ScatterAndEdepEventSet)
        // {
        //     InputDataTree->GetEntryWithIndex(evID);
        //     double r = sqrt(InputData.PosX * InputData.PosX +
        //                     InputData.PosY * InputData.PosY +
        //                     InputData.PosZ * InputData.PosZ);
        //     if (r > 0)
        //     {
        //         double phi = atan2(InputData.PosY, InputData.PosX);
        //         double cosTheta = InputData.PosZ / r;
        //         h2_origin_selected->Fill(phi, cosTheta);
        //     }
        // }

        TCanvas *c1 = new TCanvas("c1", "Moon Neutron Spectrum; Energy (MeV); Energy #times Intensity", 800, 600);
        gPad->SetLogx();
        gPad->SetLogy();
        h1_ip->Draw("HIST");
        h1_ip->Scale(Todsdt);
        c1->SaveAs("../image/InputParticleEnergy.pdf");

        // // 検出効率の位置依存性マップ（選別 / 全体）
        // TH2F *h2_origin_eff = (TH2F *)h2_origin_selected->Clone("h2_origin_eff");
        // h2_origin_eff->SetTitle("Detection Efficiency on Sphere;#phi (rad);cos(#theta)");
        // h2_origin_eff->Divide(h2_origin_all);

        // TH1D *pj = h2_origin_selected->ProjectionY("pj", 1, phi_bins, "e");

        // gStyle->SetOptStat(0);

        // TCanvas *c2 = new TCanvas("c2", "Initial Position Heatmap", 1800, 500);
        // c2->Divide(3, 1);

        // c2->cd(1);
        // h2_origin_all->Draw("COLZ");

        // c2->cd(2);
        // h2_origin_selected->Draw("COLZ");

        // c2->cd(3);
        // pj->Draw("HIST");

        // c2->SaveAs("../image/DelayedCoincidence_PositionHeatmap.pdf");
    }

    return;
}
