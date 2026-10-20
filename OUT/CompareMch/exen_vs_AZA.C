// ============================================================
// exen_vs_AZA.C
//
// <E*/A> vs spectator residue characteristics: A_res, Z_res,
// alpha_res = A_res/A_init. Modes: ld=3 (ALADIN), ld=4 (Damped
// Ericson), ld=8 (Hybrid GS+ALADIN), plus DCM for comparison.
//
// AAMCC (side A): A_res = sum(Aa_cl), y = Ex_En_per_nucleon
// (MeV, event-level).
//
// DCM: projectile prefragment (A1 = Xe124) excitation BEFORE
// evaporation, from the .r12 file (fields A Z S E*[GeV] ...).
// E*/A = E*[GeV]*1000 / A_res (MeV per nucleon). Events selected
// with Mch>0 (charged midrapidity particle in the .r12 particle
// list), as in dcm_hist_mch.C. alpha_res = A_res/124 for DCM.
//
// Run:
//   root -l -q 'OUT/CompareMch/exen_vs_AZA.C'
// ============================================================
#include <TFile.h>
#include <TTree.h>
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

void exen_vs_AZA()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadLeftMargin(0.09);

    const int Ainit = 130;      // AAMCC Xe130
    const int AinitDCM = 124;   // DCM projectile Xe124
    const double yCut = 0.5;    // spectator rapidity cut (CM)

    struct Mode {
        const char* file;   // empty for DCM
        const char* label;
        Color_t color;
    };
    const int NM = 4;
    Mode modes[NM] = {
        {"OUT/scan_mode1_ld3/d25_v030/xexe130.root",   "ld=3 (ALADIN)",           kGreen+2},
        {"OUT/scan_mode1_ld4/d25_v030/xexe130.root",   "ld=4 (Damped Ericson)",   kMagenta+2},
        {"OUT/scan_mode1_hyb23/d25_v030/xexe130.root", "ld=8 (Hybrid GS+ALADIN)", kBlue+2},
        {"",                                           "DCM Xe124+Xe130",         kRed+2}
    };

    TProfile* pA[NM]  = {nullptr};
    TProfile* pZ[NM]  = {nullptr};
    TProfile* pAl[NM] = {nullptr};

    for (int m = 0; m < NM; ++m) {
        pA[m]  = new TProfile(Form("pA_%d", m),  "", 131, -0.5, 130.5);
        pZ[m]  = new TProfile(Form("pZ_%d", m),  "", 55, -0.5, 54.5);
        pAl[m] = new TProfile(Form("pAl_%d", m), "", 200, 0.0, 1.0);
    }

    // ---- AAMCC modes ----
    for (int m = 0; m < 3; ++m) {
        TFile* f = TFile::Open(modes[m].file);
        if (!f || f->IsZombie()) {
            std::cout << "Cannot open: " << modes[m].file << std::endl;
            continue;
        }
        TTree* g = (TTree*)f->Get("Glauber");
        TTree* c = (TTree*)f->Get("MST-Clusters");
        if (!g || !c) {
            std::cout << "Trees missing in " << modes[m].file << std::endl;
            f->Close();
            continue;
        }
        float exA = 0.f;
        std::vector<int>* Aa = nullptr;
        std::vector<int>* Za = nullptr;
        g->SetBranchAddress("Ex_En_per_nucleon", &exA);
        c->SetBranchAddress("Aa_cl", &Aa);
        c->SetBranchAddress("Za_cl", &Za);

        Long64_t n = g->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            g->GetEntry(i);
            c->GetEntry(i);
            if (!Aa || !Za || Aa->empty()) continue;
            int aRes = 0, zRes = 0;
            for (size_t k = 0; k < Aa->size(); ++k) {
                aRes += (*Aa)[k];
                zRes += (*Za)[k];
            }
            if (aRes < 2) continue;
            pA[m]->Fill(aRes, exA);
            pZ[m]->Fill(zRes, exA);
            pAl[m]->Fill(double(aRes) / double(Ainit), exA);
        }
        f->Close();
        std::cout << modes[m].label << ": " << n << " events" << std::endl;
    }

    // ---- DCM (from .r12: prefragment excitation BEFORE evaporation) ----
    // Format (hypcoa-b1.f, write 104/105):
    //   8-line text header
    //   per event:  header "NC12 npart MV12 b bx by"
    //               2 prefragment lines "A Z S E*[GeV] px py pz"
    //                 (1st = projectile A1=Xe124, 2nd = target)
    //               MV12 particle lines "id ie is ib x y z px py pz am"
    {
        const int m = 3;
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
                double pAres = 0, pZres = 0, pEexc = 0;   // A, Z, E*[GeV]
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

                // particle lines: Mch>0 selection (charged midrapidity)
                bool mch = false;
                for (int i = 0; i < mv; ++i) {
                    if (!std::getline(infile, line)) break;
                    int id, ie, is, ib;
                    double x, y, z, px, py, pz, am;
                    std::istringstream iss(line);
                    if (!(iss >> id >> ie >> is >> ib
                              >> x >> y >> z >> px >> py >> pz >> am)) continue;
                    if (ie == 0) continue;
                    double p2 = px * px + py * py + pz * pz;
                    double E = std::sqrt(p2 + am * am);
                    if (E <= std::fabs(pz) + 1e-9) continue;
                    double yRap = 0.5 * std::log((E + pz) / (E - pz));
                    if (yRap > -yCut && yRap < yCut) { mch = true; }
                }
                if (!mch) continue;
                if (!haveProj || pAres < 2) continue;
                ++nSel;

                double exA = pEexc * 1000.0 / pAres;   // MeV per nucleon
                pA[m]->Fill(pAres, exA);
                pZ[m]->Fill(pZres, exA);
                pAl[m]->Fill(pAres / double(AinitDCM), exA);
            }
            infile.close();
            std::cout << modes[m].label << ": " << nSel << "/" << nAll
                      << " events (Mch>0)" << std::endl;
        }
    }

    for (int m = 0; m < NM; ++m) {
        for (TProfile* p : {pA[m], pZ[m], pAl[m]}) {
            p->SetLineColor(modes[m].color);
            p->SetMarkerColor(modes[m].color);
            p->SetMarkerStyle(20 + m);
            p->SetMarkerSize(0.9);
            p->SetLineWidth(2);
        }
    }

    const char* outDir = "OUT/CompareMch";
    gSystem->mkdir(outDir, kTRUE);

    auto drawAll = [&](const char* canvasName, const char* title,
                       TProfile** prof, const char* xaxis,
                       const char* outFile) {
        TCanvas* cv = new TCanvas(canvasName, canvasName, 1400, 900);
        prof[0]->SetTitle(title);
        prof[0]->GetXaxis()->SetTitle(xaxis);
        prof[0]->GetYaxis()->SetTitle("#LT E*/A #GT (MeV)");
        prof[0]->GetXaxis()->SetTitleSize(0.04);
        prof[0]->GetYaxis()->SetTitleSize(0.04);
        prof[0]->GetXaxis()->SetTitleOffset(1.1);
        prof[0]->GetYaxis()->SetTitleOffset(1.2);
        prof[0]->GetXaxis()->SetLabelSize(0.035);
        prof[0]->GetYaxis()->SetLabelSize(0.035);

        double ymax = 0;
        double ymin = 1e30;
        for (int m = 0; m < NM; ++m) {
            double binmax = prof[m]->GetMaximum();
            if (binmax > ymax) ymax = binmax;
            for (int b = 1; b <= prof[m]->GetNbinsX(); ++b) {
                double v = prof[m]->GetBinContent(b);
                if (v > 0 && v < ymin) ymin = v;
            }
        }
        if (ymin >= 1e30 || ymin <= 0) ymin = 1e-3;
        cv->SetLogy();
        prof[0]->SetMinimum(ymin * 0.5);
        prof[0]->SetMaximum(ymax * 1.25);

        for (int m = 0; m < NM; ++m)
            prof[m]->Draw(m == 0 ? "E1" : "E1 SAME");

        TLegend* leg = new TLegend(0.62, 0.62, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.85);
        leg->SetTextSize(0.03);
        for (int m = 0; m < NM; ++m) leg->AddEntry(prof[m], modes[m].label, "lp");
        leg->Draw();

        cv->SaveAs(TString(outDir) + "/" + outFile);
        std::cout << "Saved " << outDir << "/" << outFile << std::endl;
    };

    drawAll("cA", Form("#LT E*/A #GT vs A_{res};A_{res};#LT E*/A #GT (MeV)"),
            pA, "A_{res}", "exen_vs_A_348.png");
    drawAll("cZ", Form("#LT E*/A #GT vs Z_{res};Z_{res};#LT E*/A #GT (MeV)"),
            pZ, "Z_{res}", "exen_vs_Z_348.png");
    drawAll("cal", Form("#LT E*/A #GT vs #alpha_{res};#alpha_{res} = A_{res}/A_{init};#LT E*/A #GT (MeV)"),
            pAl, "#alpha_{res}", "exen_vs_alpha_348.png");
}
