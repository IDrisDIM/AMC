// ============================================================
// compare_dcm_ld368.C
// Compare DCM, AMC mode1 (d=2.5, V0=0.3) with ld=3 (ALADIN),
// ld=6 (Corrected ALADIN), ld=8 (hard GS+ALADIN).
//
// Usage:
//   root -l -q 'compare_dcm_ld368.C()'
//
// Produces 4 comparison PNGs (legend only on the first pad):
//   - Z3_kinematics_DCM_ld368.png
//   - Z2_kinematics_DCM_ld368.png
//   - fragmentation_moments_DCM_ld368.png
//   - fragmentation_DCM_ld368.png
// ============================================================

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <vector>
#include <iostream>
#include <algorithm>

// ======================= Стиль легенды =======================
static void StyleLeg(TLegend *leg)
{
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite, 0.8);
    leg->SetTextSize(0.032);
    leg->SetTextFont(42);
}

// ======================= Общий Y-диапазон для TH1 =======================
static void GetCommonYRange(TH1 *h1, TH1 *h2, TH1 *h3, TH1 *h4,
                            double &ymin, double &ymax,
                            bool logScale = false)
{
    auto minPos = [](TH1 *h) -> double {
        double m = 1e30;
        if (!h) return m;
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double c = h->GetBinContent(i);
            if (c > 0 && c < m) m = c;
        }
        return m;
    };
    double minPosAll = 1e30;
    double maxAll = 0;
    for (TH1 *h : {h1, h2, h3, h4}) {
        if (!h) continue;
        minPosAll = std::min(minPosAll, minPos(h));
        maxAll = std::max(maxAll, h->GetMaximum());
    }
    if (minPosAll > 1e29) minPosAll = 1e-6;
    if (maxAll <= 0) maxAll = 1.0;

    if (logScale) {
        ymin = minPosAll * 0.1;
        ymax = maxAll * 2.0;
    } else {
        ymin = 0.0;
        ymax = maxAll * 1.15;
    }
    if (ymin >= ymax) ymin = ymax * 1e-3;
}

// ======================= Общий Y-диапазон для TProfile =======================
static void GetCommonProfileYRange(TProfile *p1, TProfile *p2, TProfile *p3, TProfile *p4,
                                   double &ymin, double &ymax,
                                   bool crossZero = true)
{
    double vmin = 1e30, vmax = -1e30;
    auto scan = [&](TProfile *p) {
        if (!p) return;
        for (int i = 1; i <= p->GetNbinsX(); ++i) {
            if (p->GetBinEntries(i) < 1) continue;
            double v = p->GetBinContent(i);
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
        }
    };
    scan(p1); scan(p2); scan(p3); scan(p4);
    if (vmin > vmax) { vmin = -1.0; vmax = 1.0; }
    if (crossZero) { vmin = std::min(vmin, 0.0); vmax = std::max(vmax, 0.0); }
    double range = vmax - vmin;
    double pad = 0.10 * range;
    if (pad < 1e-6) pad = 0.05;
    ymin = vmin - pad;
    ymax = vmax + pad;
}

// ======================= Рисование 1D с общим диапазоном =======================
static void Compare1D(std::vector<TH1*> hists,
                      const std::vector<TString>& names,
                      const char* title, bool logY, bool drawLegend)
{
    if (hists.empty()) return;

    Color_t colors[10]  = {kRed+2, kBlue+2, kGreen+2, kMagenta+2, kOrange+2,
                           kCyan+2, kViolet+2, kPink+2, kTeal+2, kSpring+2};
    Style_t markers[10] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

    for (size_t i = 0; i < hists.size(); ++i) {
        TH1* h = hists[i];
        h->SetLineColor(colors[i % 10]);
        h->SetMarkerColor(colors[i % 10]);
        h->SetMarkerStyle(markers[i % 10]);
        h->SetMarkerSize(1.);
        h->SetLineWidth(1);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
        h->SetTitle(title);
    }

    // Общий Y-диапазон
    double ymin, ymax;
    TH1 *a = hists.size() > 0 ? hists[0] : nullptr;
    TH1 *b = hists.size() > 1 ? hists[1] : nullptr;
    TH1 *c = hists.size() > 2 ? hists[2] : nullptr;
    TH1 *d = hists.size() > 3 ? hists[3] : nullptr;
    GetCommonYRange(a, b, c, d, ymin, ymax, logY);

    if (logY) gPad->SetLogy();
    hists[0]->GetYaxis()->SetRangeUser(ymin, ymax);

    for (size_t i = 0; i < hists.size(); ++i) {
        if (i == 0) hists[i]->Draw("HIST E1");
        else        hists[i]->Draw("HIST E1 SAME");
    }

    if (drawLegend) {
        TLegend* leg = new TLegend(0.36, 0.55, 0.68, 0.88);
        for (size_t i = 0; i < hists.size(); ++i)
            leg->AddEntry(hists[i], names[i], "lp");
        StyleLeg(leg);
        leg->Draw();
    }
}

// ======================= Рисование TProfile с общим диапазоном =======================
static void CompareProfile(std::vector<TProfile*> profs,
                           const std::vector<TString>& names,
                           const char* title, bool drawLegend, bool crossZero)
{
    if (profs.empty()) return;

    Color_t colors[10]  = {kRed+2, kBlue+2, kGreen+2, kMagenta+2, kOrange+2,
                           kCyan+2, kViolet+2, kPink+2, kTeal+2, kSpring+2};
    Style_t markers[10] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

    for (size_t i = 0; i < profs.size(); ++i) {
        TProfile* p = profs[i];
        p->SetMarkerStyle(markers[i % 10]);
        p->SetMarkerColor(colors[i % 10]);
        p->SetLineColor(colors[i % 10]);
        p->SetMarkerSize(1.2);
        p->SetLineWidth(2);
        p->SetTitle(title);
    }

    double ymin, ymax;
    TProfile *a = profs.size() > 0 ? profs[0] : nullptr;
    TProfile *b = profs.size() > 1 ? profs[1] : nullptr;
    TProfile *c = profs.size() > 2 ? profs[2] : nullptr;
    TProfile *d = profs.size() > 3 ? profs[3] : nullptr;
    GetCommonProfileYRange(a, b, c, d, ymin, ymax, crossZero);
    profs[0]->GetYaxis()->SetRangeUser(ymin, ymax);

    for (size_t i = 0; i < profs.size(); ++i) {
        if (i == 0) profs[i]->Draw("PE1");
        else        profs[i]->Draw("PE1 SAME");
    }

    if (crossZero) {
        double xmin = profs[0]->GetXaxis()->GetXmin();
        double xmax = profs[0]->GetXaxis()->GetXmax();
        TLine* l = new TLine(xmin, 0, xmax, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
    }

    if (drawLegend) {
        TLegend* leg = new TLegend(0.36, 0.55, 0.68, 0.88);
        for (size_t i = 0; i < profs.size(); ++i)
            leg->AddEntry(profs[i], names[i], "lp");
        StyleLeg(leg);
        leg->Draw();
    }
}

// ======================= Основная функция =======================
void compare_dcm_ld368()
{
    // ==================== USER CONFIGURATION ====================
    std::vector<TString> filePaths = {
        "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch/DCM_mch.root",
        "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch/MODE1LD3_mch.root",
        "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch/MODE1LD6_mch.root",
        "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch/MODE1HYB23_mch.root",
    };

    std::vector<TString> modelNames = {
        "DCM Xe124+Xe130 (Mch>0)",
        "AMC ld=3 (ALADIN) (Mch>0)",
        "AMC ld=6 (Corrected ALADIN) (Mch>0)",
        "AMC ld=8 (GS+ALADIN) (Mch>0)",
    };

    TString outDir = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch";
    // ============================================================

    if (filePaths.size() != modelNames.size()) {
        std::cerr << "ERROR: filePaths and modelNames must have the same size!" << std::endl;
        return;
    }
    const int N = (int)filePaths.size();
    if (N == 0) { std::cerr << "ERROR: empty file list!" << std::endl; return; }

    gSystem->mkdir(outDir, kTRUE);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);

    // ==================== CANVAS 1: Z3 kinematics ====================
    {
        TCanvas* c = new TCanvas("cZ3hs", "Z=3 kinematics comparison", 1600, 1200);
        c->Divide(2, 2);

        const char* hnames[4]  = {"hY_Z3", "hEta_Z3", "hP_Z3", "hPt_Z3"};
        const char* titles[4]  = {
            "Rapidity (Z=3); y; Normalized",
            "Pseudorapidity (Z=3); #eta; Normalized",
            "|p| (Z=3); |p| (GeV/c); Normalized",
            "p_{T} (Z=3); p_{T} (GeV/c); Normalized"
        };

        for (int pad = 0; pad < 4; ++pad) {
            c->cd(pad + 1);
            std::vector<TH1*> hists;
            std::vector<TString> namesFound;
            for (int i = 0; i < N; ++i) {
                TFile* f = TFile::Open(filePaths[i], "READ");
                if (!f || f->IsZombie()) {
                    std::cerr << "Cannot open " << filePaths[i] << std::endl;
                    continue;
                }
                TH1* h = (TH1*)f->Get(hnames[pad]);
                if (h) { h->SetDirectory(0); hists.push_back(h); namesFound.push_back(modelNames[i]); }
                else std::cerr << "No " << hnames[pad] << " in " << filePaths[i] << std::endl;
                f->Close();
            }
            Compare1D(hists, namesFound, titles[pad], false, pad == 0);
        }
        c->SaveAs(outDir + "/Z3_kinematics_DCM_ld368.png");
        
    }

    // ==================== CANVAS 2: Z2 kinematics ====================
    {
        TCanvas* c = new TCanvas("cZ2hs", "Z=2 kinematics comparison", 1600, 1200);
        c->Divide(2, 2);

        const char* hnames[4]  = {"hY_Z2", "hEta_Z2", "hP_Z2", "hPt_Z2"};
        const char* titles[4]  = {
            "Rapidity (Z=2); y; Normalized",
            "Pseudorapidity (Z=2); #eta; Normalized",
            "|p| (Z=2); |p| (GeV/c); Normalized",
            "p_{T} (Z=2); p_{T} (GeV/c); Normalized"
        };

        for (int pad = 0; pad < 4; ++pad) {
            c->cd(pad + 1);
            std::vector<TH1*> hists;
            std::vector<TString> namesFound;
            for (int i = 0; i < N; ++i) {
                TFile* f = TFile::Open(filePaths[i], "READ");
                if (!f || f->IsZombie()) continue;
                TH1* h = (TH1*)f->Get(hnames[pad]);
                if (h) { h->SetDirectory(0); hists.push_back(h); namesFound.push_back(modelNames[i]); }
                f->Close();
            }
            Compare1D(hists, namesFound, titles[pad], false, pad == 0);
        }
        c->SaveAs(outDir + "/Z2_kinematics_DCM_ld368.png");
        
    }

    // ==================== CANVAS 3: Fragmentation moments ====================
    {
        TCanvas* c = new TCanvas("cMomhs", "Frag moments DCM vs ld3/6/8", 1800, 1200);
        c->Divide(3, 2);

        const char* pnames[6] = {
            "pN3_vs_Zb2", "pN3_vs_Zb3", "pNfrag_vs_Zb2",
            "pZ2_vs_Zb2", "pA_vs_Zb2", "pZ_vs_Zb2"
        };
        const char* titles[6] = {
            "<N(Z=3)> vs Zb2; Zb2; <N(Z=3)>",
            "<N(Z=3)> vs Zb3; Zb3; <N(Z=3)>",
            "<N_{frag}> vs Zb2; Zb2; <N_{frag}> (Z#geq2, |y|>0.5)",
            "<Z^{2}> vs Zb2; Zb2; <Z^{2}>(Z#geq2, |y|>0.5)",
            "<A> vs Zb2; Zb2; <A> (Z#geq2, |y|>0.5)",
            "<Z> vs Zb2; Zb2; <Z> (Z#geq2, |y|>0.5)"
        };

        for (int pad = 0; pad < 6; ++pad) {
            c->cd(pad + 1);
            std::vector<TProfile*> profs;
            std::vector<TString> namesFound;
            for (int i = 0; i < N; ++i) {
                TFile* f = TFile::Open(filePaths[i], "READ");
                if (!f || f->IsZombie()) continue;
                TProfile* p = (TProfile*)f->Get(pnames[pad]);
                if (p) { p->SetDirectory(0); profs.push_back(p); namesFound.push_back(modelNames[i]); }
                else std::cerr << "No " << pnames[pad] << " in " << filePaths[i] << std::endl;
                f->Close();
            }
            // crossZero=false для <N>, <A>, <Z> — они всегда положительны
            bool crossZero = (pad == 3 || pad == 5);  // только для <Z^2> и <Z>
            CompareProfile(profs, namesFound, titles[pad], pad == 0, crossZero);
        }
        c->SaveAs(outDir + "/fragmentation_moments_DCM_ld368.png");
        
    }

    // ==================== CANVAS 4: Fragmentation ====================
    {
        TCanvas* c = new TCanvas("cFraths", "Fragmentation DCM vs ld3/6/8", 1800, 1200);
        c->Divide(3, 2);

        // Пады 1-2: 1D (лог. шкала)
        {
            const char* hnames[2] = {"hA_frag", "hZ_frag"};
            const char* titles[2] = {
                "Mass number (Z#geq2, |y|>0.5); A; Normalized",
                "Charge (Z#geq2, |y|>0.5); Z; Normalized"
            };
            for (int pad = 0; pad < 2; ++pad) {
                c->cd(pad + 1);
                std::vector<TH1*> hists;
                std::vector<TString> namesFound;
                for (int i = 0; i < N; ++i) {
                    TFile* f = TFile::Open(filePaths[i], "READ");
                    if (!f || f->IsZombie()) continue;
                    TH1* h = (TH1*)f->Get(hnames[pad]);
                    if (h) { h->SetDirectory(0); hists.push_back(h); namesFound.push_back(modelNames[i]); }
                    f->Close();
                }
                Compare1D(hists, namesFound, titles[pad], true, pad == 0);
            }
        }

        // Пады 3-6: TProfile
        {
            const char* pnames[4] = {
                "pZmax_vs_Zb2", "pZmax_vs_Zb3",
                "pN2_vs_Zb2", "pN2_vs_Zb3"
            };
            const char* titles[4] = {
                "Z_{max} vs Zb2; Zb2; <Z_{max}> (|y|>0.5)",
                "Z_{max} vs Zb3; Zb3; <Z_{max}> (|y|>0.5)",
                "<N(Z=2)> vs Zb2; Zb2; <N(Z=2)> (|y|>0.5)",
                "<N(Z=2)> vs Zb3; Zb3; <N(Z=2)>"
            };
            for (int pad = 0; pad < 4; ++pad) {
                c->cd(pad + 3);
                std::vector<TProfile*> profs;
                std::vector<TString> namesFound;
                for (int i = 0; i < N; ++i) {
                    TFile* f = TFile::Open(filePaths[i], "READ");
                    if (!f || f->IsZombie()) continue;
                    TProfile* p = (TProfile*)f->Get(pnames[pad]);
                    if (p) { p->SetDirectory(0); profs.push_back(p); namesFound.push_back(modelNames[i]); }
                    f->Close();
                }
                CompareProfile(profs, namesFound, titles[pad], false, false);
            }
        }
        c->SaveAs(outDir + "/fragmentation_DCM_ld368.png");
       
    }

    std::cout << "\nComparison plots saved in: " << outDir << std::endl;
    std::cout << "  - Z3_kinematics_DCM_ld368.png"         << std::endl;
    std::cout << "  - Z2_kinematics_DCM_ld368.png"         << std::endl;
    std::cout << "  - fragmentation_moments_DCM_ld368.png" << std::endl;
    std::cout << "  - fragmentation_DCM_ld368.png"         << std::endl;
}


#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TString.h>
#include <cmath>
#include <iostream>

void mf_transition_compare(double lowBound = 3.0, double upBound = 5.0)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // ---- Параметры ----
    const int N = 400;
    const double E_min = 0.0;
    const double E_max = 15.0;

    // ---- График 1: Hyperbolic tangent (как в оригинале AAMCC) ----
    // aE = 1/(2*(upBound - lowBound))
    // E0 = (upBound + lowBound)/2
    // transF = 0.5*tanh((E - E0)/aE) + 0.5
    double aE_tanh = 1.0 / (2.0 * (upBound - lowBound));
    double E0_tanh = (upBound + lowBound) / 2.0;

    TGraph* g_tanh = new TGraph(N + 1);
    g_tanh->SetName("g_tanh");
    g_tanh->SetTitle(Form("MF transition probability, window [%.1f, %.1f] MeV/nucleon",
                          lowBound, upBound));

    for (int i = 0; i <= N; ++i) {
        double E = E_min + (E_max - E_min) * i / N;
        double transF = 0.5 * std::tanh((E - E0_tanh) / aE_tanh) + 0.5;
        g_tanh->SetPoint(i, E, transF);
    }
    g_tanh->SetLineColor(kRed + 2);
    //g_tanh->SetLineWidth(2.5);
    g_tanh->SetLineStyle(1);

    // ---- График 2: Linear interpolation ----
    // transF = 0 при E < lowBound
    // transF = 1 при E > upBound
    // transF линейно между
    TGraph* g_lin = new TGraph(N + 1);
    g_lin->SetName("g_lin");

    for (int i = 0; i <= N; ++i) {
        double E = E_min + (E_max - E_min) * i / N;
        double transF;
        if (E <= lowBound)      transF = 0.0;
        else if (E >= upBound)  transF = 1.0;
        else                    transF = (E - lowBound) / (upBound - lowBound);
        g_lin->SetPoint(i, E, transF);
    }
    g_lin->SetLineColor(kBlue + 2);
    //g_lin->SetLineWidth(2.5);
    g_lin->SetLineStyle(2);

    // ---- График 3: Резкий порог (для справки) ----
    // transF = 0 при E < E0
    // transF = 1 при E > E0
    TGraph* g_step = new TGraph(N + 1);
    g_step->SetName("g_step");
    for (int i = 0; i <= N; ++i) {
        double E = E_min + (E_max - E_min) * i / N;
        double transF = (E < E0_tanh) ? 0.0 : 1.0;
        g_step->SetPoint(i, E, transF);
    }
    g_step->SetLineColor(kGray + 2);
    //g_step->SetLineWidth(1.5);
    g_step->SetLineStyle(3);

    // ---- Канвас ----
    TCanvas* c = new TCanvas("c_mf_transition",
                             Form("MF transition: [%.1f, %.1f] MeV/nucleon", lowBound, upBound),
                             1200, 800);

    g_tanh->GetXaxis()->SetTitle("E*/A (MeV/nucleon)");
    g_tanh->GetYaxis()->SetTitle("P(MF)");
    g_tanh->GetXaxis()->SetTitleSize(0.045);
    g_tanh->GetYaxis()->SetTitleSize(0.045);
    g_tanh->GetXaxis()->SetTitleOffset(1.1);
    g_tanh->GetYaxis()->SetTitleOffset(1.2);
    g_tanh->GetYaxis()->SetRangeUser(-0.05, 1.05);
    g_tanh->GetXaxis()->SetRangeUser(E_min, E_max);

    g_tanh->Draw("AL");
    g_lin->Draw("L SAME");
    g_step->Draw("L SAME");

    // ---- Вертикальные линии на границах окна ----
    TLine* l_low = new TLine(lowBound, -0.05, lowBound, 1.05);
    l_low->SetLineStyle(3);
    l_low->SetLineColor(kGray + 1);
    l_low->Draw("same");

    TLine* l_up = new TLine(upBound, -0.05, upBound, 1.05);
    l_up->SetLineStyle(3);
    l_up->SetLineColor(kGray + 1);
    l_up->Draw("same");

    // ---- Легенда ----
    TLegend* leg = new TLegend(0.55, 0.15, 0.88, 0.35);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite, 0.85);
    leg->SetTextSize(0.032);
    leg->AddEntry(g_tanh, Form("tanh (aE=%.3f, E0=%.2f)", aE_tanh, E0_tanh), "l");
    leg->AddEntry(g_lin,  "linear interpolation", "l");
    leg->AddEntry(g_step, "sharp step", "l");
    leg->Draw();

    // ---- Дополнительная информация ----
    // Ширина переходной зоны
    double deltaE_tanh_90 = 2.0 * aE_tanh * std::atanh(0.9); // ±90% диапазон
    double deltaE_lin_100 = upBound - lowBound;

    // Смещение центра tanh
    double center_tanh = E0_tanh;
    double center_lin  = 0.5 * (lowBound + upBound);

    std::cout << "---- MF transition parameters ----" << std::endl;
    std::cout << "Window: lowBound = " << lowBound << " MeV/n, upBound = " << upBound << " MeV/n" << std::endl;
    std::cout << "tanh: aE = " << aE_tanh << ", E0 = " << E0_tanh << std::endl;
    std::cout << "tanh: 10-90% width = " << 2.0 * aE_tanh * std::atanh(0.8)
              << " MeV/n (approximately)" << std::endl;
    std::cout << "linear: full transition width = " << deltaE_lin_100 << " MeV/n" << std::endl;
    std::cout << "linear: center = " << center_lin << " MeV/n" << std::endl;

    c->SaveAs(Form("mf_transition_low%.0f_up%.0f.png", lowBound, upBound));
    std::cout << "Saved: mf_transition_low"
              << lowBound << "_up" << upBound << ".png" << std::endl;
}
