// ============================================================
// compare_models.C
// Compare histograms from different AAMCC output files.
//
// Usage:
//   root -l -q 'compare_models.C()'
//   or from ROOT prompt: .x compare_models.C
//
// Produces 4 comparison PNGs (legend only on the first pad):
//   - Z3_kinematics_comparison.png
//   - Z2_kinematics_comparison.png
//   - fragmentation_moments_comparison.png
//   - fragmentation_comparison.png
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
static void GetCommonYRange(TH1 *h1, TH1 *h2, TH1 *h3,
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
    for (TH1 *h : {h1, h2, h3}) {
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
static void GetCommonProfileYRange(TProfile *p1, TProfile *p2, TProfile *p3,
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
    scan(p1); scan(p2); scan(p3);
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
        h->SetMarkerSize(1.1);
        h->SetLineWidth(2);
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
        h->SetTitle(title);
    }

    // Общий Y-диапазон (учитываем до 3 кривых)
    double ymin, ymax;
    TH1 *a = hists.size() > 0 ? hists[0] : nullptr;
    TH1 *b = hists.size() > 1 ? hists[1] : nullptr;
    TH1 *c = hists.size() > 2 ? hists[2] : nullptr;
    GetCommonYRange(a, b, c, ymin, ymax, logY);

    if (logY) gPad->SetLogy();
    hists[0]->GetYaxis()->SetRangeUser(ymin, ymax);

    for (size_t i = 0; i < hists.size(); ++i) {
        if (i == 0) hists[i]->Draw("HIST E1");
        else        hists[i]->Draw("HIST E1 SAME");
    }

    if (drawLegend) {
        TLegend* leg = new TLegend(0.58, 0.62, 0.92, 0.88);
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
    GetCommonProfileYRange(a, b, c, ymin, ymax, crossZero);
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
        TLegend* leg = new TLegend(0.58, 0.62, 0.92, 0.88);
        for (size_t i = 0; i < profs.size(); ++i)
            leg->AddEntry(profs[i], names[i], "lp");
        StyleLeg(leg);
        leg->Draw();
    }
}

// ======================= Основная функция =======================
void compare_models()
{
    // ==================== USER CONFIGURATION ====================
    std::vector<TString> filePaths = {
        //"/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/XeXe_densEn1/Out_XeXe_EnergyDens1_mcini_.root",
        "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/XeXe_densEn7/Out_XeXe_EnergyDens7_mcini_.root",
        "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/XeXe_densEn4dmin21/OutUnigen_XeXe_dmin2p1_mcini_.root",
        "/home/oem/WORKLAST/DCMSMMmodel/output_hist.root",
    };

    std::vector<TString> modelNames = {
        "EnergyDens7",
        "EnergyDens4 2.1",
        "DCM",
    };

    TString outDir = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/ComparePlots";
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
        TCanvas* c = new TCanvas("cZ3cmp", "Z=3 kinematics comparison", 1600, 1200);
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
            for (int i = 0; i < N; ++i) {
                TFile* f = TFile::Open(filePaths[i], "READ");
                if (!f || f->IsZombie()) {
                    std::cerr << "Cannot open " << filePaths[i] << std::endl;
                    continue;
                }
                TH1* h = (TH1*)f->Get(hnames[pad]);
                if (h) { h->SetDirectory(0); hists.push_back(h); }
                else std::cerr << "No " << hnames[pad] << " in " << filePaths[i] << std::endl;
                f->Close();
            }
            Compare1D(hists, modelNames, titles[pad], false, pad == 0);
        }
        c->SaveAs(outDir + "/Z3_kinematics_comparison.png");
        
    }

    // ==================== CANVAS 2: Z2 kinematics ====================
    {
        TCanvas* c = new TCanvas("cZ2cmp", "Z=2 kinematics comparison", 1600, 1200);
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
            for (int i = 0; i < N; ++i) {
                TFile* f = TFile::Open(filePaths[i], "READ");
                if (!f || f->IsZombie()) continue;
                TH1* h = (TH1*)f->Get(hnames[pad]);
                if (h) { h->SetDirectory(0); hists.push_back(h); }
                f->Close();
            }
            Compare1D(hists, modelNames, titles[pad], false, pad == 0);
        }
        c->SaveAs(outDir + "/Z2_kinematics_comparison.png");
        
    }

    // ==================== CANVAS 3: Fragmentation moments ====================
    {
        TCanvas* c = new TCanvas("cFragMomcmp", "Fragmentation moments comparison", 1800, 1200);
        c->Divide(3, 2);

        const char* pnames[6] = {
            "pN3_vs_Zb2", "pN3_vs_Zb3", "pNfrag_vs_Zb2",
            "pZ2_vs_Zb2", "pA_vs_Zb2", "pZ_vs_Zb2"
        };
        const char* titles[6] = {
            "<N(Z=3)> vs Zb2; Zb2; <N(Z=3)>",
            "<N(Z=3)> vs Zb3; Zb3; <N(Z=3)>",
            "<N_{frag}> vs Zb2; Zb2; <N_{frag}> (Z#geq2)",
            "<Z^{2}> vs Zb2; Zb2; <Z^{2}>(Z#geq2)",
            "<A> vs Zb2; Zb2; <A> (Z#geq2)",
            "<Z> vs Zb2; Zb2; <Z> (Z#geq2)"
        };

        for (int pad = 0; pad < 6; ++pad) {
            c->cd(pad + 1);
            std::vector<TProfile*> profs;
            for (int i = 0; i < N; ++i) {
                TFile* f = TFile::Open(filePaths[i], "READ");
                if (!f || f->IsZombie()) continue;
                TProfile* p = (TProfile*)f->Get(pnames[pad]);
                if (p) { p->SetDirectory(0); profs.push_back(p); }
                else std::cerr << "No " << pnames[pad] << " in " << filePaths[i] << std::endl;
                f->Close();
            }
            // crossZero=false для <N>, <A>, <Z> — они всегда положительны
            bool crossZero = (pad == 3 || pad == 5);  // только для <Z^2> и <Z>
            CompareProfile(profs, modelNames, titles[pad], pad == 0, crossZero);
        }
        c->SaveAs(outDir + "/fragmentation_moments_comparison.png");
        
    }

    // ==================== CANVAS 4: Fragmentation ====================
    {
        TCanvas* c = new TCanvas("cFragcmp", "Fragmentation comparison", 1800, 1200);
        c->Divide(3, 2);

        // Пады 1-2: 1D (лог. шкала)
        {
            const char* hnames[2] = {"hA_frag", "hZ_frag"};
            const char* titles[2] = {
                "Mass number (Z#geq2); A; Normalized",
                "Charge (Z#geq2); Z; Normalized"
            };
            for (int pad = 0; pad < 2; ++pad) {
                c->cd(pad + 1);
                std::vector<TH1*> hists;
                for (int i = 0; i < N; ++i) {
                    TFile* f = TFile::Open(filePaths[i], "READ");
                    if (!f || f->IsZombie()) continue;
                    TH1* h = (TH1*)f->Get(hnames[pad]);
                    if (h) { h->SetDirectory(0); hists.push_back(h); }
                    f->Close();
                }
                Compare1D(hists, modelNames, titles[pad], true, pad == 0);
            }
        }

        // Пады 3-6: TProfile
        {
            const char* pnames[4] = {
                "pZmax_vs_Zb2", "pZmax_vs_Zb3",
                "pN2_vs_Zb2", "pN2_vs_Zb3"
            };
            const char* titles[4] = {
                "Z_{max} vs Zb2; Zb2; <Z_{max}>",
                "Z_{max} vs Zb3; Zb3; <Z_{max}>",
                "<N(Z=2)> vs Zb2; Zb2; <N(Z=2)>",
                "<N(Z=2)> vs Zb3; Zb3; <N(Z=2)>"
            };
            for (int pad = 0; pad < 4; ++pad) {
                c->cd(pad + 3);
                std::vector<TProfile*> profs;
                for (int i = 0; i < N; ++i) {
                    TFile* f = TFile::Open(filePaths[i], "READ");
                    if (!f || f->IsZombie()) continue;
                    TProfile* p = (TProfile*)f->Get(pnames[pad]);
                    if (p) { p->SetDirectory(0); profs.push_back(p); }
                    f->Close();
                }
                CompareProfile(profs, modelNames, titles[pad], false, false);
            }
        }
        c->SaveAs(outDir + "/fragmentation_comparison.png");
       
    }

    std::cout << "\nComparison plots saved in: " << outDir << std::endl;
    std::cout << "  - Z3_kinematics_comparison.png"         << std::endl;
    std::cout << "  - Z2_kinematics_comparison.png"         << std::endl;
    std::cout << "  - fragmentation_moments_comparison.png" << std::endl;
    std::cout << "  - fragmentation_comparison.png"         << std::endl;
}