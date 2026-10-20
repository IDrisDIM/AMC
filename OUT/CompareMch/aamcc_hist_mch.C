#include "/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/include/UEvent.h"
#include "/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/include/UParticle.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

R__LOAD_LIBRARY(/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/build/libMcIniData.so)

// Spectator cut: fragments (Z>=2) are counted only if |y| > yCut (y in CM frame)
static const double yCut = 0.5;

static void GetAZ_pdg(int pdg, double &A, double &Z)
{
    int ap = std::abs(pdg);
    if (ap > 1000000000)
    {
        A = (ap / 10) % 1000;
        Z = (ap / 10000) % 1000;
        if (pdg < 0) Z = -Z;
        return;
    }
    TParticlePDG *pp = TDatabasePDG::Instance()->GetParticle(pdg);
    double q = pp ? pp->Charge() / 3.0 : 0.0;
    Z = (std::fabs(q) < 0.5) ? 0.0 : (q > 0 ? 1.0 : -1.0);
    if (ap == 2212 || ap == 2112) A = 1;
    else A = 0;
}

static bool IsHadron_pdg(int pdg)
{
    int ap = std::abs(pdg);
    if (ap > 1000000000) return true;
    return ap >= 100;
}

void aamcc_hist_mch(const char *filename = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/scan_ld2_dmin/d30_v030/xexe130_mcini_.root",
                    const char *outname = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch/AMC_mch.root")
{
    TFile fin(filename);
    if (!fin.IsOpen()) { std::cerr << "Cannot open " << filename << std::endl; return; }
    TTree *t = (TTree *)fin.Get("events");
    if (!t) { std::cerr << "No events tree" << std::endl; return; }
    UEvent *ev = nullptr;
    t->SetBranchAddress("event", &ev);

    TFile fout(outname, "RECREATE");

    TH1F *hA_frag = new TH1F("hA_frag", "Mass number (Z#geq2); A; Counts", 200, 0.5, 200.5);
    TH1F *hZ_frag = new TH1F("hZ_frag", "Charge (Z#geq2); Z; Counts", 100, 0.5, 100.5);
    TH1F *hY_Z2 = new TH1F("hY_Z2", "Rapidity (Z=2); y; Counts", 60, -3, 3);
    TH1F *hEta_Z2 = new TH1F("hEta_Z2", "Pseudo-rapidity (Z=2); #eta; Counts", 120, -5, 5);
    TH1F *hP_Z2 = new TH1F("hP_Z2", "Total momentum (Z=2); |p| (GeV/c); Counts", 100, 0, 50);
    TH1F *hPt_Z2 = new TH1F("hPt_Z2", "Transverse momentum (Z=2); p_{T} (GeV/c); Counts", 100, 0, 2);
    TH1F *hY_Z3 = new TH1F("hY_Z3", "Rapidity (Z=3); y; Counts", 60, -3, 3);
    TH1F *hEta_Z3 = new TH1F("hEta_Z3", "Pseudo-rapidity (Z=3); #eta; Counts", 120, -5, 5);
    TH1F *hP_Z3 = new TH1F("hP_Z3", "Total momentum (Z=3); |p| (GeV/c); Counts", 100, 0, 50);
    TH1F *hPt_Z3 = new TH1F("hPt_Z3", "Transverse momentum (Z=3); p_{T} (GeV/c); Counts", 100, 0, 2);
    TProfile *pZmax_vs_Zb2 = new TProfile("pZmax_vs_Zb2", "Z_{max} vs Zb2; Zb2 = #Sigma Z (Z#geq2); <Z_{max}>", 80, 0, 160);
    TProfile *pZmax_vs_Zb3 = new TProfile("pZmax_vs_Zb3", "Z_{max} vs Zb3; Zb3 = #Sigma Z (Z#geq3); <Z_{max}>", 80, 0, 160);
    TProfile *pN2_vs_Zb2 = new TProfile("pN2_vs_Zb2", "<N(Z=2)> vs Zb2; Zb2; <N(Z=2)>", 80, 0, 160);
    TProfile *pN2_vs_Zb3 = new TProfile("pN2_vs_Zb3", "<N(Z=2)> vs Zb3; Zb3; <N(Z=2)>", 80, 0, 160);
    TProfile *pN3_vs_Zb2 = new TProfile("pN3_vs_Zb2", "<N(Z=3)> vs Zb2; Zb2; <N(Z=3)>", 80, 0, 160);
    TProfile *pN3_vs_Zb3 = new TProfile("pN3_vs_Zb3", "<N(Z=3)> vs Zb3; Zb3; <N(Z=3)>", 80, 0, 160);
    TProfile *pZ2_vs_Zb2 = new TProfile("pZ2_vs_Zb2", "<Z^{2}> vs Zb2; Zb2; <Z^{2}>(Z#geq2)", 80, 0, 160);
    TProfile *pNfrag_vs_Zb2 = new TProfile("pNfrag_vs_Zb2", "<N_{frag}> vs Zb2; Zb2; <N_{frag}> (Z#geq2)", 80, 0, 160);
    TProfile *pA_vs_Zb2 = new TProfile("pA_vs_Zb2", "<A> vs Zb2; Zb2; <A> (Z#geq2)", 80, 0, 160);
    TProfile *pZ_vs_Zb2 = new TProfile("pZ_vs_Zb2", "<Z> vs Zb2; Zb2; <Z> (Z#geq2)", 80, 0, 160);

    Long64_t n = t->GetEntries();
    long nSel = 0;
    for (Long64_t j = 0; j < n; ++j)
    {
        t->GetEntry(j);
        TClonesArray *arr = ev->GetParticleList();
        int np = arr->GetEntries();

        bool mch = false;
        for (int i = 0; i < np; ++i)
        {
            UParticle *p = (UParticle *)arr->At(i);
            int pdg = p->GetPdg();
            if (!IsHadron_pdg(pdg)) continue;
            double A, Z;
            GetAZ_pdg(pdg, A, Z);
            if (Z == 0) continue;
            double pz = p->Pz(), E = p->E();
            if (E > std::fabs(pz) + 1e-9)
            {
                double y = 0.5 * std::log((E + pz) / (E - pz));
                if (y > -0.5 && y < 0.5) { mch = true; break; }
            }
        }
        if (!mch) continue;
        ++nSel;

        double psiRP = ev->GetPhi();
        double cosPhi = std::cos(psiRP);
        double sinPhi = std::sin(psiRP);

        int sumZb2 = 0, sumZb3 = 0, maxZ = 0, countZ2 = 0, countZ3 = 0, nFrag = 0;
        double sumZ2 = 0, sumA_frag = 0, sumZ_frag = 0;

        for (int i = 0; i < np; ++i)
        {
            UParticle *p = (UParticle *)arr->At(i);
            int pdg = p->GetPdg();
            double A, Z;
            GetAZ_pdg(pdg, A, Z);

            double px = p->Px(), py = p->Py(), pz = p->Pz(), E = p->E();
            double p2 = px * px + py * py + pz * pz;
            double pmod = std::sqrt(p2);
            double pxRP = px * cosPhi + py * sinPhi;
            double pyRP = -px * sinPhi + py * cosPhi;
            double pt = std::sqrt(pxRP * pxRP + pyRP * pyRP);

            int Zi_i = (int)Z;
            if (Zi_i == 2)
            {
                if (E > std::fabs(pz) + 1e-9)
                    hY_Z2->Fill(0.5 * std::log((E + pz) / (E - pz)));
                if (pmod > std::fabs(pz) + 1e-9)
                    hEta_Z2->Fill(0.5 * std::log((pmod + pz) / (pmod - pz)));
                hP_Z2->Fill(pmod);
                hPt_Z2->Fill(pt);
            }
            if (Zi_i == 3)
            {
                if (E > std::fabs(pz) + 1e-9)
                    hY_Z3->Fill(0.5 * std::log((E + pz) / (E - pz)));
                if (pmod > std::fabs(pz) + 1e-9)
                    hEta_Z3->Fill(0.5 * std::log((pmod + pz) / (pmod - pz)));
                hP_Z3->Fill(pmod);
                hPt_Z3->Fill(pt);
            }
            // Spectator selection: fragments (Z>=2) outside midrapidity, |y| > yCut
            bool isSpectator = false;
            if (Zi_i >= 2 && E > std::fabs(pz) + 1e-9)
            {
                double yf = 0.5 * std::log((E + pz) / (E - pz));
                if (std::fabs(yf) > yCut) isSpectator = true;
            }
            if (isSpectator)
            {
                sumZb2 += Zi_i;
                nFrag++;
                sumZ2 += Zi_i * Zi_i;
                sumA_frag += A;
                sumZ_frag += Zi_i;
                hA_frag->Fill(A);
                hZ_frag->Fill(Z);
                if (Zi_i > maxZ) maxZ = Zi_i;
                if (Zi_i >= 3) sumZb3 += Zi_i;
                if (Zi_i == 2) countZ2++;
                if (Zi_i == 3) countZ3++;
            }
        }

        pZmax_vs_Zb2->Fill(sumZb2, maxZ);
        pZmax_vs_Zb3->Fill(sumZb3, maxZ);
        pN2_vs_Zb2->Fill(sumZb2, countZ2);
        pN2_vs_Zb3->Fill(sumZb3, countZ2);
        pN3_vs_Zb2->Fill(sumZb2, countZ3);
        pN3_vs_Zb3->Fill(sumZb3, countZ3);
        pNfrag_vs_Zb2->Fill(sumZb2, nFrag);
        if (nFrag > 0)
        {
            pZ2_vs_Zb2->Fill(sumZb2, sumZ2 / nFrag);
            pA_vs_Zb2->Fill(sumZb2, sumA_frag / nFrag);
            pZ_vs_Zb2->Fill(sumZb2, sumZ_frag / nFrag);
        }
    }

    std::cout << "AAMCC: total=" << n << " selected(Mch>0)=" << nSel << std::endl;

    fout.cd();
    hA_frag->Write(); hZ_frag->Write();
    hY_Z2->Write(); hEta_Z2->Write(); hP_Z2->Write(); hPt_Z2->Write();
    hY_Z3->Write(); hEta_Z3->Write(); hP_Z3->Write(); hPt_Z3->Write();
    pZmax_vs_Zb2->Write(); pZmax_vs_Zb3->Write();
    pN2_vs_Zb2->Write(); pN2_vs_Zb3->Write();
    pN3_vs_Zb2->Write(); pN3_vs_Zb3->Write();
    pZ2_vs_Zb2->Write(); pNfrag_vs_Zb2->Write();
    pA_vs_Zb2->Write(); pZ_vs_Zb2->Write();
    fout.Close();
}
