#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"

// Spectator cut: fragments (Z>=2) are counted only if |y| > yCut (y in CM frame)
static const double yCut = 0.5;

struct Particle
{
    double A, Z, S, px, py, pz, mass, extra;
};
struct Event
{
    int num;
    double b, bx, by;
    std::vector<Particle> particles;
};

void dcm_hist_mch(const char *filename = "/home/oem/WORKLAST/DCMSMMmodel/runAuAu/XeXe/XeXe_event100.out",
                  const char *outname = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/CompareMch/DCM_mch.root")
{
    std::ifstream infile(filename);
    if (!infile.is_open()) { std::cerr << "Cannot open " << filename << std::endl; return; }

    std::vector<Event> events;
    std::string line;
    int eventNum; double b, bx, by;

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (!(iss >> eventNum >> b >> bx >> by)) continue;
        Event ev; ev.num = eventNum; ev.b = b; ev.bx = bx; ev.by = by;
        if (!std::getline(infile, line)) break;
        int npart;
        { std::istringstream iss2(line); if (!(iss2 >> npart)) continue; }
        for (int i = 0; i < npart; ++i)
        {
            if (!std::getline(infile, line)) break;
            std::istringstream iss3(line);
            Particle p;
            if (!(iss3 >> p.A >> p.Z >> p.S >> p.px >> p.py >> p.pz >> p.mass >> p.extra)) continue;
            ev.particles.push_back(p);
        }
        events.push_back(ev);
    }
    infile.close();

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

    long nSel = 0, nAll = 0;
    for (auto &ev : events)
    {
        ++nAll;
        bool mch = false;
        for (auto &p : ev.particles)
        {
            if (p.Z == 0) continue;
            double p2 = p.px * p.px + p.py * p.py + p.pz * p.pz;
            double E = std::sqrt(p2 + p.mass * p.mass);
            if (E > std::fabs(p.pz) + 1e-9)
            {
                double y = 0.5 * std::log((E + p.pz) / (E - p.pz));
                if (y > -0.5 && y < 0.5) { mch = true; break; }
            }
        }
        if (!mch) continue;
        ++nSel;

        double psiRP = std::atan2(ev.by, ev.bx);
        double cosPhi = std::cos(psiRP);
        double sinPhi = std::sin(psiRP);

        int sumZb2 = 0, sumZb3 = 0, maxZ = 0, countZ2 = 0, countZ3 = 0, nFrag = 0;
        double sumZ2 = 0, sumA_frag = 0, sumZ_frag = 0;

        for (auto &p : ev.particles)
        {
            double p2 = p.px * p.px + p.py * p.py + p.pz * p.pz;
            double pmod = std::sqrt(p2);
            double pxRP = p.px * cosPhi + p.py * sinPhi;
            double pyRP = -p.px * sinPhi + p.py * cosPhi;
            double pt = std::sqrt(pxRP * pxRP + pyRP * pyRP);
            double E = std::sqrt(p2 + p.mass * p.mass);

            int Zi_i = (int)p.Z;

            if (Zi_i == 2)
            {
                if (E > std::fabs(p.pz) + 1e-9)
                    hY_Z2->Fill(0.5 * std::log((E + p.pz) / (E - p.pz)));
                if (pmod > std::fabs(p.pz) + 1e-9)
                    hEta_Z2->Fill(0.5 * std::log((pmod + p.pz) / (pmod - p.pz)));
                hP_Z2->Fill(pmod);
                hPt_Z2->Fill(pt);
            }
            if (Zi_i == 3)
            {
                if (E > std::fabs(p.pz) + 1e-9)
                    hY_Z3->Fill(0.5 * std::log((E + p.pz) / (E - p.pz)));
                if (pmod > std::fabs(p.pz) + 1e-9)
                    hEta_Z3->Fill(0.5 * std::log((pmod + p.pz) / (pmod - p.pz)));
                hP_Z3->Fill(pmod);
                hPt_Z3->Fill(pt);
            }
            // Spectator selection: fragments (Z>=2) outside midrapidity, |y| > yCut
            bool isSpectator = false;
            if (Zi_i >= 2 && E > std::fabs(p.pz) + 1e-9)
            {
                double yf = 0.5 * std::log((E + p.pz) / (E - p.pz));
                if (std::fabs(yf) > yCut) isSpectator = true;
            }
            if (isSpectator)
            {
                sumZb2 += Zi_i;
                nFrag++;
                sumZ2 += Zi_i * Zi_i;
                sumA_frag += p.A;
                sumZ_frag += Zi_i;
                hA_frag->Fill(p.A);
                hZ_frag->Fill(p.Z);
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

    std::cout << "DCM: total=" << nAll << " selected(Mch>0)=" << nSel << std::endl;

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
