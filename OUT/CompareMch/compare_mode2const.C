// Compare mode 2 (2 passes, constant V0, growing dmin) at d=2.5 and d=3.0
// against DCM Xe124+Xe130 and the mode-1 reference, using Mch>0 selection.
//
// Runs: V0=0.3, ld=2 (Gaimard-Schmidt), G4 de-excitation.
//
// Outputs:
//   OUT/CompareMch/mode2const_table.txt
//   OUT/CompareMch/mode2const_Z.png
//   OUT/CompareMch/mode2const_ratio.png

R__LOAD_LIBRARY(/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/build/libMcIniData.so)

#include "/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/include/UEvent.h"
#include "/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/include/UParticle.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TStyle.h"
#include "TMath.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static const char *kDCMFile = "/home/oem/WORKLAST/DCMSMMmodel/scr/XeXe_3p8/XeXe_event100.out";

struct Stats {
    long nsel = 0;
    double nfrag = 0;
    double sumZ = 0;
    double z35 = 0;
    double z45 = 0;
    double z50 = 0;
};

struct Variant {
    std::string label;
    std::string file;
    Color_t color;
    Style_t marker;
    Style_t line;
};

static void GetAZ_pdg(int pdg, double &A, double &Z)
{
    int ap = std::abs(pdg);
    if (ap > 1000000000) {
        A = (ap / 10) % 1000;
        Z = (ap / 10000) % 1000;
        if (pdg < 0) Z = -Z;
        return;
    }
    TParticlePDG *pp = TDatabasePDG::Instance()->GetParticle(pdg);
    double q = pp ? pp->Charge() / 3.0 : 0.0;
    Z = (std::fabs(q) < 0.5) ? 0.0 : (q > 0 ? 1.0 : -1.0);
    A = (ap == 2212 || ap == 2112) ? 1 : 0;
}

static bool IsHadron_pdg(int pdg)
{
    int ap = std::abs(pdg);
    if (ap > 1000000000) return true;
    return ap >= 100;
}

static TH1F *BookZ(const char *name)
{
    TH1F *h = new TH1F(name, "Fragments (Z#geq2); Z; dN/dZ per selected event", 100, 0.5, 100.5);
    h->SetDirectory(0);
    return h;
}

static TH1F *AnalyzeAMC(const char *fn, const char *name, Stats &st)
{
    TFile fin(fn);
    if (!fin.IsOpen()) { std::cerr << "Cannot open " << fn << std::endl; return nullptr; }
    TTree *t = (TTree *)fin.Get("events");
    if (!t) { std::cerr << "No events tree in " << fn << std::endl; return nullptr; }

    UEvent *ev = nullptr;
    t->SetBranchAddress("event", &ev);

    TH1F *hZ = BookZ(name);
    Long64_t n = t->GetEntries();

    for (Long64_t j = 0; j < n; ++j) {
        t->GetEntry(j);
        TClonesArray *arr = ev->GetParticleList();
        int np = arr->GetEntries();

        bool mch = false;
        for (int i = 0; i < np; ++i) {
            UParticle *p = (UParticle *)arr->At(i);
            int pdg = p->GetPdg();
            if (!IsHadron_pdg(pdg)) continue;
            double A, Z;
            GetAZ_pdg(pdg, A, Z);
            if (Z == 0) continue;
            double pz = p->Pz(), E = p->E();
            if (E > std::fabs(pz) + 1e-9) {
                double y = 0.5 * std::log((E + pz) / (E - pz));
                if (y > -0.5 && y < 0.5) { mch = true; break; }
            }
        }
        if (!mch) continue;
        ++st.nsel;

        for (int i = 0; i < np; ++i) {
            UParticle *p = (UParticle *)arr->At(i);
            double A, Z;
            GetAZ_pdg(p->GetPdg(), A, Z);
            int Zi = (int)Z;
            if (Zi < 2) continue;
            hZ->Fill(Zi);
            st.nfrag += 1;
            st.sumZ += Zi;
            if (Zi >= 35) st.z35 += 1;
            if (Zi >= 45) st.z45 += 1;
            if (Zi >= 50) st.z50 += 1;
        }
    }
    return hZ;
}

struct DcmParticle { double A, Z, S, px, py, pz, mass, extra; };
struct DcmEvent { int num; double b, bx, by; std::vector<DcmParticle> particles; };

static TH1F *AnalyzeDCM(const char *fn, const char *name, Stats &st)
{
    std::ifstream infile(fn);
    if (!infile.is_open()) { std::cerr << "Cannot open " << fn << std::endl; return nullptr; }

    TH1F *hZ = BookZ(name);
    std::string line;
    int eventNum; double b, bx, by;
    long nAll = 0;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (!(iss >> eventNum >> b >> bx >> by)) continue;
        DcmEvent ev; ev.num = eventNum; ev.b = b; ev.bx = bx; ev.by = by;
        if (!std::getline(infile, line)) break;
        int npart;
        { std::istringstream iss2(line); if (!(iss2 >> npart)) continue; }
        for (int i = 0; i < npart; ++i) {
            if (!std::getline(infile, line)) break;
            std::istringstream iss3(line);
            DcmParticle p;
            if (!(iss3 >> p.A >> p.Z >> p.S >> p.px >> p.py >> p.pz >> p.mass >> p.extra)) continue;
            ev.particles.push_back(p);
        }
        ++nAll;

        bool mch = false;
        for (auto &p : ev.particles) {
            if (p.Z == 0) continue;
            double p2 = p.px * p.px + p.py * p.py + p.pz * p.pz;
            double E = std::sqrt(p2 + p.mass * p.mass);
            if (E > std::fabs(p.pz) + 1e-9) {
                double y = 0.5 * std::log((E + p.pz) / (E - p.pz));
                if (y > -0.5 && y < 0.5) { mch = true; break; }
            }
        }
        if (!mch) continue;
        ++st.nsel;

        for (auto &p : ev.particles) {
            int Zi = (int)p.Z;
            if (Zi < 2) continue;
            hZ->Fill(Zi);
            st.nfrag += 1;
            st.sumZ += Zi;
            if (Zi >= 35) st.z35 += 1;
            if (Zi >= 45) st.z45 += 1;
            if (Zi >= 50) st.z50 += 1;
        }
    }
    std::cout << "DCM: total=" << nAll << " selected(Mch>0)=" << st.nsel << std::endl;
    return hZ;
}

// Variable Z bins used for chi2 and ratio plots.
static const int kEdges[] = {2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 101};
static const int kNEdges = sizeof(kEdges) / sizeof(kEdges[0]);

static double BinRate(const TH1 *h, long nsel, int lo, int hi)
{
    if (nsel <= 0) return 0;
    return h->Integral(lo, hi - 1) / double(nsel);
}

static long BinCount(const TH1 *h, int lo, int hi)
{
    return (long)h->Integral(lo, hi - 1);
}

static double Chi2NDF(const TH1 *hA, long nA, const TH1 *hB, long nB,
                      const char *labelA, const char *labelB, int &ndf)
{
    double chi2 = 0;
    ndf = 0;
    for (int i = 0; i < kNEdges - 1; ++i) {
        int lo = kEdges[i], hi = kEdges[i + 1];
        double NA = BinCount(hA, lo, hi);
        double NB = BinCount(hB, lo, hi);
        if (NA + NB <= 0) continue;
        double rA = NA / double(nA);
        double rB = NB / double(nB);
        double varA = NA / (double(nA) * double(nA));
        double varB = NB / (double(nB) * double(nB));
        double var = varA + varB;
        if (var <= 0) continue;
        double d = rA - rB;
        chi2 += d * d / var;
        ++ndf;
    }
    if (ndf == 0) return 0;
    std::cout << "  chi2(" << labelA << " vs " << labelB << ") = "
              << chi2 << " / ndf " << ndf << " = " << chi2 / ndf << std::endl;
    return chi2 / ndf;
}

static void PrintStats(FILE *f, const char *label, const Stats &st, double chi2ndf)
{
    fprintf(f, "%-14s %6ld  %7.3f  %6.2f  %7.3f  %7.3f  %7.3f",
            label, st.nsel, st.nfrag / st.nsel, st.sumZ / st.nfrag,
            st.z35 / st.nsel, st.z45 / st.nsel, st.z50 / st.nsel);
    if (chi2ndf >= 0) fprintf(f, "  %8.2f", chi2ndf);
    fprintf(f, "\n");
    printf("%-14s %6ld  %7.3f  %6.2f  %7.3f  %7.3f  %7.3f",
           label, st.nsel, st.nfrag / st.nsel, st.sumZ / st.nfrag,
           st.z35 / st.nsel, st.z45 / st.nsel, st.z50 / st.nsel);
    if (chi2ndf >= 0) printf("  chi2/ndf=%6.2f", chi2ndf);
    printf("\n");
}

void compare_mode2const()
{
    gStyle->SetOptStat(0);

    const std::string outdir = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch";
    const std::string base = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT";

    std::vector<Variant> vars = {
        {"mode2 d2.5", base + "/scan_mode2_v0const/d25_v030/xexe130_mcini_.root",     kRed + 1,    20, 1},
        {"mode2 d3.0", base + "/scan_mode2_v0const/d30_v030/xexe130_mcini_.root",     kBlue + 1,   21, 1},
        {"mode1 d3.0 ref", base + "/scan_ld2_dmin/d30_v030/xexe130_mcini_.root",      kGray + 2,   25, 7}
    };

    Stats stDCM;
    TH1F *hDCM = AnalyzeDCM(kDCMFile, "hZ_DCM", stDCM);
    if (!hDCM || stDCM.nsel == 0) { std::cerr << "DCM analysis failed" << std::endl; return; }

    std::vector<TH1F *> hists;
    std::vector<Stats> stats;
    hists.push_back(hDCM);
    stats.push_back(stDCM);
    for (auto &v : vars) {
        Stats st;
        TH1F *h = AnalyzeAMC(v.file.c_str(), ("hZ_" + v.label).c_str(), st);
        if (!h || st.nsel == 0) {
            std::cerr << "Failed: " << v.label << " (" << v.file << ")" << std::endl;
            hists.push_back(nullptr);
            stats.push_back(st);
            continue;
        }
        hists.push_back(h);
        stats.push_back(st);
        std::cout << v.label << ": selected=" << st.nsel << std::endl;
    }

    // chi2 vs DCM
    std::vector<double> chi2(vars.size(), -1);
    std::vector<int> ndf(vars.size(), 0);
    std::cout << "\nChi2 vs DCM (per selected event, variable Z bins):\n";
    for (size_t i = 0; i < vars.size(); ++i) {
        if (!hists[i + 1]) continue;
        chi2[i] = Chi2NDF(hists[i + 1], stats[i + 1].nsel, hDCM, stDCM.nsel,
                          vars[i].label.c_str(), "DCM", ndf[i]);
    }

    // Pairwise: mode 2 d=2.5 vs d=3.0
    int pairNdf = 0;
    double pairChi2 = -1;
    if (hists[1] && hists[2]) {
        std::cout << "\nPairwise mode2 d2.5 vs d3.0:\n";
        pairChi2 = Chi2NDF(hists[1], stats[1].nsel,
                           hists[2], stats[2].nsel,
                           "mode2 d2.5", "mode2 d3.0", pairNdf);
    }

    // ---- table ----
    FILE *ftxt = fopen((outdir + "/mode2const_table.txt").c_str(), "w");
    if (ftxt) {
        fprintf(ftxt, "Event selection: Mch>0 (charged hadron in -0.5<y<0.5)\n");
        fprintf(ftxt, "DCM: Xe-124 + Xe-130, 3.8 AGeV, 5000 events\n");
        fprintf(ftxt, "AAMCC: UrQMD input xexe130, 7895 events, mode 2 (2 passes, constant V0, growing dmin)\n\n");
        fprintf(ftxt, "mode2 d2.5: d=2.5 fm, V0=0.3, ld=2, G4\n");
        fprintf(ftxt, "mode2 d3.0: d=3.0 fm, V0=0.3, ld=2, G4\n");
        fprintf(ftxt, "Reference: mode 1, d=3.0 fm, V0=0.3, ld=2, G4\n\n");
        fprintf(ftxt, "%-14s %6s  %7s  %6s  %7s  %7s  %7s  %8s\n",
                "variant", "Nsel", "n/ev", "<Z>", "Z>=35", "Z>=45", "Z>=50", "chi2/ndf");
        PrintStats(ftxt, "DCM", stDCM, -1);
        for (size_t i = 0; i < vars.size(); ++i)
            if (hists[i + 1]) PrintStats(ftxt, vars[i].label.c_str(), stats[i + 1], chi2[i]);
        fprintf(ftxt, "\nPairwise mode2 d2.5 vs d3.0: chi2/ndf = %.2f (ndf=%d)\n",
                pairChi2, pairNdf);

        // per-bin rates and ratios
        fprintf(ftxt, "\nPer selected event rates dN/dZ (variable Z bins):\n");
        fprintf(ftxt, "%-8s %10s", "Z", "DCM");
        for (auto &v : vars) fprintf(ftxt, " %10s", v.label.c_str());
        fprintf(ftxt, "\n");
        for (int i = 0; i < kNEdges - 1; ++i) {
            int lo = kEdges[i], hi = kEdges[i + 1] - 1;
            fprintf(ftxt, "%-8s %10.5f", Form("%d-%d", lo, hi), BinRate(hDCM, stDCM.nsel, lo, hi + 1));
            for (size_t k = 0; k < vars.size(); ++k) {
                if (!hists[k + 1]) { fprintf(ftxt, " %10s", "-"); continue; }
                fprintf(ftxt, " %10.5f", BinRate(hists[k + 1], stats[k + 1].nsel, lo, hi + 1));
            }
            fprintf(ftxt, "\n");
        }

        fprintf(ftxt, "\nRatios AMC/DCM (per selected event):\n");
        fprintf(ftxt, "%-8s", "Z");
        for (auto &v : vars) fprintf(ftxt, " %10s", v.label.c_str());
        fprintf(ftxt, "\n");
        for (int i = 0; i < kNEdges - 1; ++i) {
            int lo = kEdges[i], hi = kEdges[i + 1] - 1;
            double rD = BinRate(hDCM, stDCM.nsel, lo, hi + 1);
            fprintf(ftxt, "%-8s", Form("%d-%d", lo, hi));
            for (size_t k = 0; k < vars.size(); ++k) {
                if (!hists[k + 1]) { fprintf(ftxt, " %10s", "-"); continue; }
                double rA = BinRate(hists[k + 1], stats[k + 1].nsel, lo, hi + 1);
                if (rD > 0 && rA > 0) fprintf(ftxt, " %10.3f", rA / rD);
                else fprintf(ftxt, " %10s", "n/a");
            }
            fprintf(ftxt, "\n");
        }
        fclose(ftxt);
    }

    // ---- plot 1: per-event Z distributions ----
    TCanvas *c1 = new TCanvas("c1", "Z distributions", 1400, 1000);
    c1->SetLogy();
    c1->SetLeftMargin(0.11);
    c1->SetRightMargin(0.05);

    TH1F *hD = (TH1F *)hDCM->Clone("hD_norm");
    hD->Scale(1.0 / stDCM.nsel);
    hD->SetLineColor(kBlack);
    hD->SetLineWidth(3);
    hD->GetXaxis()->SetRangeUser(0.5, 60.5);
    hD->SetMinimum(1e-7);
    hD->SetMaximum(10);
    hD->SetTitle("Fragments Z#geq2, Mch>0; Z; dN/dZ per selected event");
    hD->Draw("HIST");

    TLegend *leg = new TLegend(0.55, 0.55, 0.93, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->AddEntry(hD, "DCM Xe124+Xe130", "l");

    for (size_t i = 0; i < vars.size(); ++i) {
        if (!hists[i + 1]) continue;
        TH1F *h = (TH1F *)hists[i + 1]->Clone(("h_" + vars[i].label).c_str());
        h->Scale(1.0 / stats[i + 1].nsel);
        h->SetLineColor(vars[i].color);
        h->SetLineWidth(2);
        h->SetLineStyle(vars[i].line);
        h->Draw("HIST SAME");
        leg->AddEntry(h, vars[i].label.c_str(), "l");
    }
    leg->Draw();
    c1->SaveAs((outdir + "/mode2const_Z.png").c_str());

    // ---- plot 2: ratios to DCM ----
    TCanvas *c2 = new TCanvas("c2", "Ratio to DCM", 1400, 1000);
    c2->SetLeftMargin(0.11);
    c2->SetRightMargin(0.05);
    c2->SetBottomMargin(0.13);

    TH1D *frame = new TH1D("frame", "; Z; AMC / DCM", 14, 0, 14);
    frame->SetMinimum(0.0);
    frame->SetMaximum(2.2);
    for (int i = 1; i <= 14; ++i) {
        int idx = i - 1;
        int lo = kEdges[idx], hi = kEdges[idx + 1] - 1;
        frame->GetXaxis()->SetBinLabel(i, Form("%d-%d", lo, hi));
    }
    frame->GetYaxis()->SetTitle("AMC / DCM (per selected event)");
    frame->Draw("HIST");
    TLine *one = new TLine(0, 1, 14, 1);
    one->SetLineStyle(2);
    one->SetLineColor(kBlack);
    one->Draw();

    TLegend *leg2 = new TLegend(0.52, 0.62, 0.93, 0.92);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.028);

    for (size_t k = 0; k < vars.size(); ++k) {
        if (!hists[k + 1]) continue;
        TGraphErrors *g = new TGraphErrors();
        int ip = 0;
        for (int i = 0; i < kNEdges - 1; ++i) {
            int lo = kEdges[i], hi = kEdges[i + 1] - 1;
            double ND = BinCount(hDCM, lo, hi + 1);
            double NA = BinCount(hists[k + 1], lo, hi + 1);
            if (ND <= 0 || NA <= 0) continue;
            double rD = ND / double(stDCM.nsel);
            double rA = NA / double(stats[k + 1].nsel);
            double ratio = rA / rD;
            double er = ratio * std::sqrt(1.0 / NA + 1.0 / ND);
            g->SetPoint(ip, i + 0.5, ratio);
            g->SetPointError(ip, 0, er);
            ++ip;
        }
        g->SetMarkerStyle(vars[k].marker);
        g->SetMarkerColor(vars[k].color);
        g->SetLineColor(vars[k].color);
        g->SetLineStyle(vars[k].line);
        g->SetLineWidth(2);
        g->Draw("PL SAME");
        leg2->AddEntry(g, vars[k].label.c_str(), "pl");
    }
    leg2->Draw();
    c2->SaveAs((outdir + "/mode2const_ratio.png").c_str());

    std::cout << "\nWrote:\n  " << outdir << "/mode2const_table.txt\n  "
              << outdir << "/mode2const_Z.png\n  "
              << outdir << "/mode2const_ratio.png" << std::endl;
}
