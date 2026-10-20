// ============================================================
// compare_dcm_sm88.C
// Compare DCM and AMC mode 1 (d=2.5, V0=0.3, ld=88 Smooth Hybrid GS+ALADIN).
//
// Usage:
//   root -l -q 'compare_dcm_sm88.C()'
//
// Produces 4 comparison PNGs (legend only on the first pad):
//   - Z3_kinematics_DCM_vs_sm88.png
//   - Z2_kinematics_DCM_vs_sm88.png
//   - fragmentation_moments_DCM_vs_sm88.png
//   - fragmentation_DCM_vs_sm88.png
// ============================================================

#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <vector>
#include <iostream>
#include <algorithm>

static void StyleLeg(TLegend *leg)
{
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite, 0.8);
    leg->SetTextSize(0.032);
    leg->SetTextFont(42);
}

static void GetCommonYRange(TH1 *h1, TH1 *h2,
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
    for (TH1 *h : {h1, h2}) {
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

static void GetCommonProfileYRange(TProfile *p1, TProfile *p2,
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
    scan(p1); scan(p2);
    if (vmin > vmax) { vmin = -1.0; vmax = 1.0; }
    if (crossZero) { vmin = std::min(vmin, 0.0); vmax = std::max(vmax, 0.0); }
    double range = vmax - vmin;
    double pad = 0.10 * range;
    if (pad < 1e-6) pad = 0.05;
    ymin = vmin - pad;
    ymax = vmax + pad;
}

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

    double ymin, ymax;
    TH1 *a = hists.size() > 0 ? hists[0] : nullptr;
    TH1 *b = hists.size() > 1 ? hists[1] : nullptr;
    GetCommonYRange(a, b, ymin, ymax, logY);

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
    GetCommonProfileYRange(a, b, ymin, ymax, crossZero);
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

void compare_dcm_sm88()
{
    // ==================== USER CONFIGURATION ====================
    std::vector<TString> filePaths = {
        "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch/DCM_mch.root",
        "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch/MODE1SM88_mch.root",
    };

    std::vector<TString> modelNames = {
        "DCM Xe124+Xe130 (Mch>0)",
        "AMC mode1 d=2.5 V0=0.3 ld=88 (Smooth GS+ALADIN) (Mch>0)",
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
        TCanvas* c = new TCanvas("cZ3hyb", "Z=3 kinematics DCM vs sm88", 1600, 1200);
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
                if (!f || f->IsZombie()) continue;
                TH1* h = (TH1*)f->Get(hnames[pad]);
                if (h) { h->SetDirectory(0); hists.push_back(h); namesFound.push_back(modelNames[i]); }
                f->Close();
            }
            Compare1D(hists, namesFound, titles[pad], false, pad == 0);
        }
        c->SaveAs(outDir + "/Z3_kinematics_DCM_vs_sm88.png");
    }

    // ==================== CANVAS 2: Z2 kinematics ====================
    {
        TCanvas* c = new TCanvas("cZ2hyb", "Z=2 kinematics DCM vs sm88", 1600, 1200);
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
        c->SaveAs(outDir + "/Z2_kinematics_DCM_vs_sm88.png");
    }

    // ==================== CANVAS 3: Fragmentation moments ====================
    {
        TCanvas* c = new TCanvas("cMomhyb", "Fragmentation moments DCM vs sm88", 1800, 1200);
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
            bool crossZero = (pad == 3 || pad == 5);
            CompareProfile(profs, namesFound, titles[pad], pad == 0, crossZero);
        }
        c->SaveAs(outDir + "/fragmentation_moments_DCM_vs_sm88.png");
    }

    // ==================== CANVAS 4: Fragmentation ====================
    {
        TCanvas* c = new TCanvas("cFrathyb", "Fragmentation DCM vs sm88", 1800, 1200);
        c->Divide(3, 2);

        // pads 1-2: 1D (log scale)
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

        // pads 3-6: TProfile
        {
            const char* pnames[4] = {
                "pZmax_vs_Zb2", "pZmax_vs_Zb3",
                "pN2_vs_Zb2", "pN2_vs_Zb3"
            };
            const char* titles[4] = {
                "Z_{max} vs Zb2; Zb2; <Z_{max}> (|y|>0.5)",
                "Z_{max} vs Zb3; Zb3; <Z_{max}> (|y|>0.5)",
                "<N(Z=2)> vs Zb2; Zb2; <N(Z=2)> (|y|>0.5)",
                "<N(Z=2)> vs Zb3; Zb3; <N(Z=2)> (|y|>0.5)"
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
        c->SaveAs(outDir + "/fragmentation_DCM_vs_sm88.png");
    }

    std::cout << "\nWrote: Z3_kinematics_DCM_vs_sm88.png, Z2_kinematics_DCM_vs_sm88.png, "
              << "fragmentation_moments_DCM_vs_sm88.png, fragmentation_DCM_vs_sm88.png" << std::endl;
}
