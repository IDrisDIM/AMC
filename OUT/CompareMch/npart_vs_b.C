// ============================================================
// npart_vs_b.C
//
// Comparison AAMCC vs DCM of:
//   G1: 1D NpartA (participants from the projectile, side A)
//   G2: 2D b vs NpartA (one panel per model)
//   G3: <Sum Z_spectators> vs b (side A)
//
// Event selection: NONE (no multiplicity cut).
// A_init(projectile) = 124 in both codes (AAMCC Conditions.Mass_on_A=124,
// DCM ANUCL1=124).
//
// AAMCC (single mode Xe124+Xe130, hybrid GS+ALADIN):
//   NpartA        = A_init - sum(Aa_cl)   (spectator clusters from MST-Clusters)
//   Sum Z_spec    = sum(Za_cl)            (spectator clusters, side A)
//   b             = Glauber.impact_parameter
//   NB: Glauber.Npart/NpartA/NpartB branches are present but written as 0 in
//       this build, so participants are derived from the spectator clusters.
//
// DCM (from .r12):
//   NpartA        = 124 - A_res (A_res = projectile prefragment from .r12)
//   Sum Z_spec    = Z_res (charge of projectile prefragment)
//
// Plots:
//   G1  npartA_1d_348.png           1D NpartA (counts)
//   G2  npartA_vs_b_348.png         2D b vs NpartA (AAMCC | DCM)
//   G3  zspec_vs_b_348.png          <Sum Z_spec> vs b
//   G4  npartA_vs_b_profile_348.png <NpartA> vs b (left) + Cb(b) (right)
//   G5  npartA_vs_cb_348.png        <NpartA> vs Cb
//
// Run:
//   root -l -q 'OUT/CompareMch/npart_vs_b.C'
// ============================================================
#include "/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/include/UEvent.h"
#include "/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/include/UParticle.h"
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdio>

R__LOAD_LIBRARY(/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/build/libMcIniData.so)

static const int Ainit = 124;    // projectile (side A) in both codes
static const double yCut = 0.5;  // midrapidity cut (CM)

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

void npart_vs_b()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadRightMargin(0.13);
    gStyle->SetPadLeftMargin(0.09);

    struct Mode {
        const char* raw;    // AAMCC raw root (clusters/glauber)
        const char* mcini;  // AAMCC mcini root (particles)
        const char* label;
        Color_t color;
    };
    const int NM = 1;
    Mode modes[NM] = {
        {"OUT/scan_mode1_hyb23/d25_v030/xexe130.root",
         "OUT/scan_mode1_hyb23/d25_v030/xexe130_mcini_.root",
         "AAMCC Xe124+Xe130", kBlue+2}
    };
    const char* dcmLabel = "DCM Xe124+Xe130";
    Color_t dcmColor = kRed+2;

    const int NTOT = NM + 1;  // + DCM

    TH1F* hNpart[NTOT] = {nullptr};
    TH2F* h2b[NTOT] = {nullptr};
    TProfile* pZs[NTOT] = {nullptr};
    TProfile* pNpartB[NTOT] = {nullptr};
    TH1F* hB[NTOT] = {nullptr};   // P(b) from all events
    TProfile* pNpartCb[NTOT] = {nullptr};  // <NpartA> vs Cb

    const int nbB = 60;  // 0..15 fm, 0.25 fm
    for (int m = 0; m < NTOT; ++m) {
        hNpart[m] = new TH1F(Form("hNpart_%d", m), "", Ainit + 1, -0.5, Ainit + 0.5);
        h2b[m]    = new TH2F(Form("h2b_%d", m), "",
                             nbB, 0.0, 15.0, Ainit + 1, -0.5, Ainit + 0.5);
        pZs[m]    = new TProfile(Form("pZs_%d", m), "", nbB, 0.0, 15.0);
        pNpartB[m] = new TProfile(Form("pNpartB_%d", m), "", nbB, 0.0, 15.0);
        hB[m]     = new TH1F(Form("hB_%d", m), "", nbB, 0.0, 15.0);
        pNpartCb[m] = new TProfile(Form("pNpartCb_%d", m), "", 100, 0.0, 1.0);
    }

    // per-event storage for the <NpartA> vs Cb profile (filled after Cb is built)
    std::vector<double> evB[NTOT];
    std::vector<double> evNpart[NTOT];

    // ---- AAMCC modes ----
    for (int m = 0; m < NM; ++m) {
        TFile* fr = TFile::Open(modes[m].raw);
        TFile* fm = TFile::Open(modes[m].mcini);
        if (!fr || fr->IsZombie() || !fm || fm->IsZombie()) {
            std::cout << "Cannot open files for " << modes[m].label << std::endl;
            continue;
        }
        TTree* tg   = (TTree*)fr->Get("Glauber");
        TTree* tc   = (TTree*)fr->Get("MST-Clusters");
        TTree* tev  = (TTree*)fm->Get("events");
        if (!tg || !tc || !tev) {
            std::cout << "Trees missing for " << modes[m].label << std::endl;
            fr->Close(); fm->Close();
            continue;
        }

        float b = 0.f;
        std::vector<int>* Aa = nullptr;
        std::vector<int>* Za = nullptr;
        tg->SetBranchAddress("impact_parameter", &b);
        tc->SetBranchAddress("Aa_cl", &Aa);
        tc->SetBranchAddress("Za_cl", &Za);

        UEvent* ev = nullptr;
        tev->SetBranchAddress("event", &ev);

        Long64_t n = tg->GetEntries();
        if (tev->GetEntries() != n)
            std::cout << "WARNING: entry mismatch raw=" << n
                      << " mcini=" << tev->GetEntries() << std::endl;

        long nSel = 0;
        for (Long64_t i = 0; i < n; ++i) {
            tg->GetEntry(i);
            tc->GetEntry(i);
            tev->GetEntry(i);

            ++nSel;

            int sumAa = 0;
            double sumZ = 0;
            if (Aa) for (size_t k = 0; k < Aa->size(); ++k) sumAa += (*Aa)[k];
            if (Za) for (size_t k = 0; k < Za->size(); ++k) sumZ += (*Za)[k];

            int npartA = Ainit - sumAa;
            if (npartA < 0) npartA = 0;

            hNpart[m]->Fill(npartA);
            h2b[m]->Fill(b, npartA);
            pZs[m]->Fill(b, sumZ);
            pNpartB[m]->Fill(b, npartA);
            hB[m]->Fill(b);
            evB[m].push_back(b);
            evNpart[m].push_back(npartA);
        }
        fr->Close();
        fm->Close();
        std::cout << modes[m].label << ": " << nSel << "/" << n
                  << " events (no cut)" << std::endl;
    }

    // ---- DCM (from .r12) ----
    {
        const int m = NM;  // DCM index
        const char* dcmFile =
            "/home/oem/WORKLAST/DCMSMMmodel/scr/XeXe_3p8/XeXe_event100.r12";
        std::ifstream infile(dcmFile);
        if (!infile.is_open()) {
            std::cerr << "Cannot open " << dcmFile << std::endl;
        } else {
            std::string line;
            for (int i = 0; i < 8; ++i) std::getline(infile, line);  // text header

            long nAll = 0, nSel = 0;
            while (std::getline(infile, line)) {
                int nc12, npart, mv;
                double b, bx, by;
                {
                    std::istringstream iss(line);
                    if (!(iss >> nc12 >> npart >> mv >> b >> bx >> by)) break;
                }
                ++nAll;

                // prefragment lines: 1st = projectile (A1, Xe124)
                double pAres = 0, pZres = 0, pEexc = 0;
                bool haveProj = false;
                for (int pf = 0; pf < 2; ++pf) {
                    if (!std::getline(infile, line)) break;
                    if (pf == 0) {
                        double S, px, py, pz;
                        std::istringstream iss(line);
                        if (iss >> pAres >> pZres >> S >> pEexc
                                >> px >> py >> pz)
                            haveProj = true;
                    }
                }

                // (Mch cut removed) still consume the particle lines
                for (int i = 0; i < mv; ++i) {
                    if (!std::getline(infile, line)) break;
                }
                if (!haveProj || pAres < 0) continue;
                ++nSel;

                int npartA_DCM = Ainit - int(pAres + 0.5);
                if (npartA_DCM < 0) npartA_DCM = 0;
                hNpart[m]->Fill(npartA_DCM);
                h2b[m]->Fill(b, npartA_DCM);
                pZs[m]->Fill(b, pZres);
                pNpartB[m]->Fill(b, npartA_DCM);
                hB[m]->Fill(b);
                evB[m].push_back(b);
                evNpart[m].push_back(npartA_DCM);
            }
            infile.close();
            std::cout << dcmLabel << ": " << nSel << "/" << nAll
                      << " events (no cut)" << std::endl;
        }
    }

    for (int m = 0; m < NTOT; ++m) {
        hNpart[m]->SetLineColor(m < NM ? modes[m].color : dcmColor);
        hNpart[m]->SetLineWidth(2);
        pZs[m]->SetLineColor(m < NM ? modes[m].color : dcmColor);
        pZs[m]->SetMarkerColor(m < NM ? modes[m].color : dcmColor);
        pZs[m]->SetMarkerStyle(20 + m);
        pZs[m]->SetMarkerSize(0.8);
        pZs[m]->SetLineWidth(2);
        pNpartB[m]->SetLineColor(m < NM ? modes[m].color : dcmColor);
        pNpartB[m]->SetMarkerColor(m < NM ? modes[m].color : dcmColor);
        pNpartB[m]->SetMarkerStyle(20 + m);
        pNpartB[m]->SetMarkerSize(0.8);
        pNpartB[m]->SetLineWidth(2);
    }
    std::vector<TString> labels;
    for (int m = 0; m < NM; ++m) labels.push_back(modes[m].label);
    labels.push_back(dcmLabel);

    // Build normalized P(b) and cumulative Cb(b) = integral_0^b P(b') db'
    // Manual accumulation over visible bins (b>=0); value stored at each bin
    // is the running sum of P up to and including that bin (right-edge value).
    // GetCumulative() was wrong here (it accumulates from the underflow bin in
    // an unnormalised way, giving Cb(b~0) ~ 1 instead of ~ 0).
    TH1F* hCb[NTOT] = {nullptr};
    const int nbx = hB[0]->GetNbinsX();
    for (int m = 0; m < NTOT; ++m) {
        double N = hB[m]->Integral(1, nbx);         // visible bins only
        if (N > 0) hB[m]->Scale(1.0 / N);           // P(b) normalized to unit area
        hCb[m] = new TH1F(Form("hCb_%d", m), "",
                          nbx, hB[m]->GetXaxis()->GetXmin(),
                          hB[m]->GetXaxis()->GetXmax());
        double sum = 0.0;
        for (int i = 1; i <= nbx; ++i) {
            sum += hB[m]->GetBinContent(i);
            hCb[m]->SetBinContent(i, sum);          // Cb(b_high) = integral_0^b_high
        }
        hCb[m]->SetLineColor(m < NM ? modes[m].color : dcmColor);
        hCb[m]->SetLineWidth(2);
        pNpartCb[m]->SetLineColor(m < NM ? modes[m].color : dcmColor);
        pNpartCb[m]->SetMarkerColor(m < NM ? modes[m].color : dcmColor);
        pNpartCb[m]->SetMarkerStyle(20 + m);
        pNpartCb[m]->SetMarkerSize(0.7);
        pNpartCb[m]->SetLineWidth(2);
    }
    // <NpartA> vs Cb : map each event's b to its cumulative value.
    // Use the previous bin's cumulative (i.e. integral_0^{b_low}) so that the
    // very first bin maps to Cb ~ 0 rather than already ~ its own weight.
    for (int m = 0; m < NTOT; ++m) {
        for (size_t k = 0; k < evB[m].size(); ++k) {
            int bin = hB[m]->FindBin(evB[m][k]);
            if (bin < 1) bin = 1;
            if (bin > nbx) bin = nbx;
            double cb = (bin > 1) ? hCb[m]->GetBinContent(bin - 1) : 0.0;
            pNpartCb[m]->Fill(cb, evNpart[m][k]);
        }
    }
    std::cout << "P(b) events (no cut):";
    for (int m = 0; m < NTOT; ++m)
        std::cout << " " << labels[m] << "=" << evB[m].size();
    std::cout << std::endl;
    for (int m = 0; m < NTOT; ++m) {
        double c0 = hCb[m]->GetBinContent(1);
        double cMid = hCb[m]->GetBinContent(nbx / 2);
        double cEnd = hCb[m]->GetBinContent(nbx);
        std::cout << "Cb sanity " << labels[m]
                  << ": first-bin=" << c0
                  << " mid=" << cMid
                  << " last=" << cEnd << std::endl;
    }

    const char* outDir = "OUT/CompareMch";
    gSystem->mkdir(outDir, kTRUE);

    // ---------- G1: 1D NpartA ----------
    {
        TCanvas* cv = new TCanvas("cNpart", "cNpart", 1000, 800);
        hNpart[0]->SetTitle("NpartA (participants from projectile);NpartA;Counts");
        hNpart[0]->GetXaxis()->SetTitleSize(0.04);
        hNpart[0]->GetYaxis()->SetTitleSize(0.04);
        hNpart[0]->GetYaxis()->SetTitleOffset(1.1);
        hNpart[0]->GetXaxis()->SetLabelSize(0.035);
        hNpart[0]->GetYaxis()->SetLabelSize(0.035);

        double ymax = 0;
        for (int m = 0; m < NTOT; ++m) ymax = std::max(ymax, hNpart[m]->GetMaximum());
        hNpart[0]->SetMaximum(ymax * 1.15);

        for (int m = 0; m < NTOT; ++m)
            hNpart[m]->Draw(m == 0 ? "HIST" : "HIST SAME");

        TLegend* leg = new TLegend(0.60, 0.65, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.85);
        leg->SetTextSize(0.03);
        for (int m = 0; m < NTOT; ++m) leg->AddEntry(hNpart[m], labels[m], "l");
        leg->Draw();

        cv->SaveAs(TString(outDir) + "/npartA_1d_348.png");
        std::cout << "Saved " << outDir << "/npartA_1d_348.png" << std::endl;
    }

    // ---------- G2: 2D b vs NpartA (two panels: AAMCC ld=8 | DCM) ----------
    {
        gStyle->SetOptStat(0);
        int mArr[2] = {0, NM};  // AAMCC and DCM
        TCanvas* cv = new TCanvas("cb", "cb", 1600, 800);
        cv->Divide(2, 1);

        double zmax = 0;
        for (int j = 0; j < 2; ++j) zmax = std::max(zmax, h2b[mArr[j]]->GetMaximum());

        for (int j = 0; j < 2; ++j) {
            int m = mArr[j];
            cv->cd(j + 1);
            gPad->SetRightMargin(0.15);
            gPad->SetLeftMargin(0.10);
            h2b[m]->SetTitle(Form("b vs NpartA - %s;b (fm);NpartA", labels[m].Data()));
            h2b[m]->GetXaxis()->SetTitleSize(0.045);
            h2b[m]->GetYaxis()->SetTitleSize(0.045);
            h2b[m]->GetXaxis()->SetLabelSize(0.038);
            h2b[m]->GetYaxis()->SetLabelSize(0.038);
            h2b[m]->GetZaxis()->SetTitle("Counts");
            h2b[m]->GetZaxis()->SetTitleSize(0.04);
            if (zmax > 0) h2b[m]->SetMaximum(zmax);
            h2b[m]->Draw("COLZ");
        }
        cv->SaveAs(TString(outDir) + "/npartA_vs_b_348.png");
        std::cout << "Saved " << outDir << "/npartA_vs_b_348.png" << std::endl;
    }

    // ---------- G3: <Sum Z_spec> vs b ----------
    {
        TCanvas* cv = new TCanvas("cZs", "cZs", 1000, 800);
        pZs[0]->SetTitle("#LT #Sigma Z_{spectators} #GT vs b;b (fm);#LT #Sigma Z_{spectators} #GT");
        pZs[0]->GetXaxis()->SetTitleSize(0.04);
        pZs[0]->GetYaxis()->SetTitleSize(0.04);
        pZs[0]->GetYaxis()->SetTitleOffset(1.2);
        pZs[0]->GetXaxis()->SetLabelSize(0.035);
        pZs[0]->GetYaxis()->SetLabelSize(0.035);
        pZs[0]->GetYaxis()->SetRangeUser(0, Ainit + 2);

        for (int m = 0; m < NTOT; ++m)
            pZs[m]->Draw(m == 0 ? "E1" : "E1 SAME");

        TLegend* leg = new TLegend(0.60, 0.65, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.85);
        leg->SetTextSize(0.03);
        for (int m = 0; m < NTOT; ++m) leg->AddEntry(pZs[m], labels[m], "lp");
        leg->Draw();

        cv->SaveAs(TString(outDir) + "/zspec_vs_b_348.png");
        std::cout << "Saved " << outDir << "/zspec_vs_b_348.png" << std::endl;
    }

    // ---------- G4: <NpartA(b)> vs b (left) + Cb(b) (right) ----------
    {
        TCanvas* cv = new TCanvas("cNpartB", "cNpartB", 1800, 800);
        cv->Divide(2, 1);

        // left: <NpartA> vs b
        cv->cd(1);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.03);
        pNpartB[0]->SetTitle("#LT NpartA #GT vs b;b (fm);#LT NpartA #GT");
        pNpartB[0]->GetXaxis()->SetTitleSize(0.045);
        pNpartB[0]->GetYaxis()->SetTitleSize(0.045);
        pNpartB[0]->GetYaxis()->SetTitleOffset(1.15);
        pNpartB[0]->GetXaxis()->SetLabelSize(0.038);
        pNpartB[0]->GetYaxis()->SetLabelSize(0.038);
        pNpartB[0]->GetYaxis()->SetRangeUser(0, Ainit + 2);
        for (int m = 0; m < NTOT; ++m)
            pNpartB[m]->Draw(m == 0 ? "E1" : "E1 SAME");
        TLegend* leg = new TLegend(0.60, 0.65, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.85);
        leg->SetTextSize(0.03);
        for (int m = 0; m < NTOT; ++m) leg->AddEntry(pNpartB[m], labels[m], "lp");
        leg->Draw();

        // right: cumulative Cb(b) = integral_0^b P(b') db', P(b) normalized
        cv->cd(2);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.03);
        hCb[0]->SetTitle("C_{b}(b) = #int_{0}^{b} P(b') db';b (fm);C_{b}");
        hCb[0]->GetXaxis()->SetTitleSize(0.045);
        hCb[0]->GetYaxis()->SetTitleSize(0.045);
        hCb[0]->GetYaxis()->SetTitleOffset(1.15);
        hCb[0]->GetXaxis()->SetLabelSize(0.038);
        hCb[0]->GetYaxis()->SetLabelSize(0.038);
        hCb[0]->GetYaxis()->SetRangeUser(0, 1.05);
        for (int m = 0; m < NTOT; ++m)
            hCb[m]->Draw(m == 0 ? "HIST" : "HIST SAME");
        TLegend* leg2 = new TLegend(0.45, 0.15, 0.92, 0.40);
        leg2->SetBorderSize(0);
        leg2->SetFillColorAlpha(kWhite, 0.85);
        leg2->SetTextSize(0.03);
        for (int m = 0; m < NTOT; ++m) leg2->AddEntry(hCb[m], labels[m], "l");
        leg2->Draw();

        cv->SaveAs(TString(outDir) + "/npartA_vs_b_profile_348.png");
        std::cout << "Saved " << outDir << "/npartA_vs_b_profile_348.png" << std::endl;
    }

    // ---------- G5: <NpartA> vs Cb ----------
    {
        TCanvas* cv = new TCanvas("cNpartCb", "cNpartCb", 1000, 800);
        gPad->SetLeftMargin(0.11);
        pNpartCb[0]->SetTitle("#LT NpartA #GT vs C_{b};C_{b};#LT NpartA #GT");
        pNpartCb[0]->GetXaxis()->SetTitleSize(0.045);
        pNpartCb[0]->GetYaxis()->SetTitleSize(0.045);
        pNpartCb[0]->GetYaxis()->SetTitleOffset(1.15);
        pNpartCb[0]->GetXaxis()->SetLabelSize(0.038);
        pNpartCb[0]->GetYaxis()->SetLabelSize(0.038);
        pNpartCb[0]->GetYaxis()->SetRangeUser(0, Ainit + 2);
        for (int m = 0; m < NTOT; ++m)
            pNpartCb[m]->Draw(m == 0 ? "E1" : "E1 SAME");
        TLegend* leg = new TLegend(0.55, 0.65, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.85);
        leg->SetTextSize(0.03);
        for (int m = 0; m < NTOT; ++m) leg->AddEntry(pNpartCb[m], labels[m], "lp");
        leg->Draw();

        cv->SaveAs(TString(outDir) + "/npartA_vs_cb_348.png");
        std::cout << "Saved " << outDir << "/npartA_vs_cb_348.png" << std::endl;
    }
}
