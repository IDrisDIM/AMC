#define readUnigen_cxx
#include "readUnigen.h"

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

// ======================= Классы центральности (BM@N) =======================
static Float_t minCentPercent[10] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
static Float_t maxCentPercent[10] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
static Float_t minB[10] = {0, 3.81907, 5.4932, 6.73342, 7.72698, 8.61946, 9.47531, 10.2768, 10.9616, 11.4984};
static Float_t maxB[10] = {3.81907, 5.4932, 6.73342, 7.72698, 8.61946, 9.47531, 10.2768, 10.9616, 11.4984, 12.0014};

static Int_t GetCentImpact_Unigen(Double_t imp)
{
    for (Int_t i = 0; i < 10; i++)
        if (imp >= minB[i] && imp < maxB[i])
            return i;
    return -1;
}

static void StyleLegend_Unigen(TLegend *leg)
{
    leg->SetBorderSize(0);
    leg->SetFillStyle(1001);
    leg->SetFillColorAlpha(kWhite, 0.75);
    leg->SetTextSize(0.026);
    leg->SetTextFont(42);
    leg->SetLineColor(kGray + 2);
    leg->SetLineWidth(0);
}

// ======================= BM@N (10-40%) =======================
static const int N_BMN_Y = 6;
static double bmn_y_v1_x[N_BMN_Y] = {0.0994, 0.2997, 0.5000, 0.6987, 0.8990, 1.0994};
static double bmn_y_v1_y[N_BMN_Y] = {0.0116, 0.0636, 0.1233, 0.2031, 0.2744, 0.3504};
static const int N_BMN_PT = 5;
static double bmn_pt_v1_x[N_BMN_PT] = {0.1965, 0.6018, 0.9992, 1.3995, 1.7972};
static double bmn_pt_v1_y[N_BMN_PT] = {0.1132, 0.2178, 0.3078, 0.3519, 0.3659};
static const int N_BMN_B = 6;
static double bmn_b_v1_x[N_BMN_B] = {2.7442, 5.0331, 6.5147, 7.7145, 8.7452, 9.6820};
static double bmn_b_v1_y[N_BMN_B] = {0.1886, 0.2479, 0.2605, 0.2508, 0.2294, 0.1993};
static double bmn_y_v2_x[N_BMN_Y] = {0.1037, 0.3040, 0.5028, 0.7032, 0.9024, 1.1014};
static double bmn_y_v2_y[N_BMN_Y] = {-0.0020, -0.0032, -0.0021, 0.0006, 0.0096, 0.0150};
static double bmn_pt_v2_x[N_BMN_PT] = {0.2069, 0.6074, 1.0048, 1.4047, 1.8047};
static double bmn_pt_v2_y[N_BMN_PT] = {-0.0061, 0.0024, -0.0005, -0.0027, -0.0023};
static double bmn_b_v2_x[N_BMN_B] = {2.7563, 5.0250, 6.5063, 7.7063, 8.7375, 9.6750};
static double bmn_b_v2_y[N_BMN_B] = {-0.0006, 0.0018, 0.0006, -0.0012, -0.0043, -0.0082};

static TGraphErrors *MakeBMnGraph_Unigen(int n, double *x, double *y, const char *name)
{
    TGraphErrors *g = new TGraphErrors(n, x, y);
    g->SetName(name);
    g->SetMarkerStyle(29);
    g->SetMarkerSize(1.8);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
    g->SetLineWidth(2);
    return g;
}

// ==================== Биннинг ====================
static const int N_Y_BINS = 6;
static const double Y_FLOW_MIN = 0.0;
static const double Y_FLOW_MAX = 1.2;
static const int N_PT_BINS = 5;
static const double PT_MIN = 0.0;
static const double PT_MAX = 2.0;

static void GetAZfromPDG(int pdg, double &A, double &Z)
{
    int apdg = std::abs(pdg);
    if (apdg >= 1000000000)
    {
        int A_ = (apdg / 10) % 1000;
        int Z_ = (apdg / 10000) % 1000;
        A = (double)A_;
        Z = (double)Z_;
        if (pdg < 0)
            Z = -Z;
        return;
    }
    switch (apdg)
    {
    case 2212:
        A = 1;
        Z = 1;
        break;
    case 2112:
        A = 1;
        Z = 0;
        break;
    case 211:
        A = 0;
        Z = 1;
        break;
    case 111:
        A = 0;
        Z = 0;
        break;
    case 321:
        A = 0;
        Z = 1;
        break;
    case 311:
        A = 0;
        Z = 0;
        break;
    case 130:
        A = 0;
        Z = 0;
        break;
    case 310:
        A = 0;
        Z = 0;
        break;
    case 3122:
        A = 0;
        Z = 0;
        break;
    case 3222:
        A = 0;
        Z = 1;
        break;
    case 3212:
        A = 0;
        Z = 0;
        break;
    case 3112:
        A = 0;
        Z = -1;
        break;
    case 3312:
        A = 0;
        Z = -1;
        break;
    case 3322:
        A = 0;
        Z = 0;
        break;
    case 3334:
        A = 0;
        Z = -1;
        break;
    case 22:
        A = 0;
        Z = 0;
        break;
    case 11:
        A = 0;
        Z = -1;
        break;
    case -11:
        A = 0;
        Z = 1;
        break;
    case 13:
        A = 0;
        Z = -1;
        break;
    case -13:
        A = 0;
        Z = 1;
        break;
    default:
        A = 0;
        Z = 0;
        break;
    }
    if (pdg < 0 && apdg != 11 && apdg != 13)
        Z = -Z;
}

static double GetMassFromPDG(int pdg, double A)
{
    auto particle = (TParticlePDG *)TDatabasePDG::Instance()->GetParticle(pdg);
    if (particle && particle->Mass() > 0.01)
        return particle->Mass();
    if (A >= 1.5)
        return A * 0.9315;
    return 0.0;
}

// ======================= Loop() =======================
// void readUnigen::Loop(const char *inputPattern = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/XeXe_densEn4dmin21/out_XeXe_dmi2p1_mcini_.root",
//                      const char *outputFileName = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/XeXe_densEn4dmin21/OutUnigen_XeXe_dmin2p1_mcini_.root")
// void readUnigen::Loop(const char *inputPattern = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/XeXe_densEn4dmin21_P300/XeXe_dmin21_p300_mcini_.root",
//                      const char *outputFileName = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/XeXe_densEn4dmin21_P300/Out_XeXe_dmin21_p300_mcini_.root")
void readUnigen::Loop(const char *inputPattern = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/out_XeXe_dmin2p5_Pmin300New2mode2newVCUT/xexe130_mcini_.root",
                      const char *outputFileName = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/out_XeXe_dmin2p5_Pmin300New2mode2newVCUT/Out_xexe130_mcini_.root")
{
    TChain *chain = new TChain("events");

    TString pattern(inputPattern);
    if (pattern.Contains("*") || pattern.Contains("?"))
    {
        TString cmd = TString::Format("ls %s 2>/dev/null", pattern.Data());
        FILE *fp = gSystem->OpenPipe(cmd.Data(), "r");
        if (!fp)
        {
            std::cerr << "Cannot open pipe for pattern " << inputPattern << std::endl;
            return;
        }
        char buf[4096];
        int nAdded = 0;
        while (fgets(buf, sizeof(buf), fp))
        {
            TString line(buf);
            line.ReplaceAll("\n", "");
            if (line.IsNull())
                continue;
            std::cout << "  add: " << line << std::endl;
            chain->Add(line);
            nAdded++;
        }
        gSystem->ClosePipe(fp);
        if (nAdded == 0)
        {
            std::cerr << "No files match pattern: " << inputPattern << std::endl;
            return;
        }
    }
    else
    {
        TObjArray *tokens = pattern.Tokenize(",");
        for (int i = 0; i < tokens->GetEntries(); ++i)
        {
            TString fname = ((TObjString *)tokens->At(i))->GetString();
            fname.ReplaceAll(" ", "");
            std::cout << "  add: " << fname << std::endl;
            chain->Add(fname);
        }
        delete tokens;
    }

    Init(chain);

    Long64_t nentries = fChain->GetEntriesFast();
    if (nentries == 0)
    {
        std::cerr << "No events in input!" << std::endl;
        return;
    }
    std::cout << "Total events: " << nentries << std::endl;

    // ==================== Директория для выходных файлов ====================
    TString outDir = gSystem->DirName(outputFileName);
    if (outDir.IsNull() || outDir == ".")
        outDir = gSystem->pwd();
    if (gSystem->AccessPathName(outDir, kWritePermission))
    {
        gSystem->mkdir(outDir, kTRUE);
    }
    std::cout << "PNG will be saved to: " << outDir << std::endl;

    // ==================== Гистограммы ====================
    TH1F *hB = new TH1F("hB", "Impact parameter; b (fm); Events", 100, 0, 20);
    TH1F *hY = new TH1F("hY", "Rapidity; y; Particles", 30, -3, 3);
    TH1F *hEta = new TH1F("hEta", "Pseudo-rapidity; #eta; Particles", 120, -5, 5);
    TH1F *hPt = new TH1F("hPt", "Transverse momentum; p_{T} (GeV/c); Particles", 25, 0, 2);
    TH1F *hP = new TH1F("hP", "Total momentum; |p| (GeV/c); Particles", 100, 0, 3);
    TH1F *hA = new TH1F("hA", "Mass number (all); A; Counts", 200, 0.5, 200.5);
    TH1F *hZ = new TH1F("hZ", "Charge (all); Z; Counts", 100, 0.5, 100.5);
    TH1F *hA_frag = new TH1F("hA_frag", "Mass number (Z#geq2); A; Counts", 200, 0.5, 200.5);
    TH1F *hZ_frag = new TH1F("hZ_frag", "Charge (Z#geq2); Z; Counts", 100, 0.5, 100.5);

    TH1F *hY_proton = new TH1F("hY_proton", "Rapidity (protons); y; Counts", 60, -3, 3);
    TH1F *hEta_proton = new TH1F("hEta_proton", "Pseudo-rapidity (protons); #eta; Counts", 120, -5, 5);
    TH1F *hP_proton = new TH1F("hP_proton", "Total momentum (protons); |p| (GeV/c); Counts", 100, 0, 3);
    TH1F *hPt_proton = new TH1F("hPt_proton", "Transverse momentum (protons); p_{T} (GeV/c); Counts", 100, 0, 2);

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

    TProfile *pPx_vs_y = new TProfile("pPx_vs_y", "<p_{x}> vs y (RP frame); y; <p_{x}> (GeV/c)", 30, -3, 3);
    TProfile *pPy_vs_y = new TProfile("pPy_vs_y", "<p_{y}> vs y (RP frame); y; <p_{y}> (GeV/c)", 30, -3, 3);
    TProfile *pPt_vs_b = new TProfile("pPt_vs_b", "<p_{T}> vs b; b (fm); <p_{T}> (GeV/c)", 20, 0, 20);

    const double V1_Y_MIN = 0.6, V1_Y_MAX = 1.2;
    const double V2_Y_MIN = 0.0, V2_Y_MAX = 1.2;
    const double V1_PT_MIN = 0.4, V1_PT_MAX = 2.0;
    const double V2_PT_MIN = 0.4, V2_PT_MAX = 2.0;

    TProfile *pV1_y_proton = new TProfile("pV1_y_proton", "v_{1}(y) protons, 10-40%; |y|; v_{1}", N_Y_BINS, Y_FLOW_MIN, Y_FLOW_MAX);
    TProfile *pV2_y_proton = new TProfile("pV2_y_proton", "v_{2}(y) protons, 10-40%; |y|; v_{2}", N_Y_BINS, Y_FLOW_MIN, Y_FLOW_MAX);
    TProfile *pV1_pt_proton = new TProfile("pV1_pt_proton", "v_{1}(p_{T}) protons, 10-40%; p_{T} (GeV/c); v_{1}", N_PT_BINS, PT_MIN, PT_MAX);
    TProfile *pV2_pt_proton = new TProfile("pV2_pt_proton", "v_{2}(p_{T}) protons, 10-40%; p_{T} (GeV/c); v_{2}", N_PT_BINS, PT_MIN, PT_MAX);
    TProfile *pV1_cent_proton = new TProfile("pV1_cent_proton", "v_{1} protons vs centrality; centrality (%); v_{1}", 10, 0, 100);
    TProfile *pV2_cent_proton = new TProfile("pV2_cent_proton", "v_{2} protons vs centrality; centrality (%); v_{2}", 10, 0, 100);

    Double_t bEdges[11];
    for (int i = 0; i < 10; ++i)
        bEdges[i] = minB[i];
    bEdges[10] = maxB[9];
    TProfile *pV1_b_proton = new TProfile("pV1_b_proton", "v_{1} protons vs <b>; <b> (fm); v_{1}", 10, bEdges);
    TProfile *pV2_b_proton = new TProfile("pV2_b_proton", "v_{2} protons vs <b>; <b> (fm); v_{2}", 10, bEdges);

    TProfile *pV1_y_proton_010 = new TProfile("pV1_y_proton_010", "v_{1}(y) protons, 0-10%; |y|; v_{1}", N_Y_BINS, Y_FLOW_MIN, Y_FLOW_MAX);
    TProfile *pV2_y_proton_010 = new TProfile("pV2_y_proton_010", "v_{2}(y) protons, 0-10%; |y|; v_{2}", N_Y_BINS, Y_FLOW_MIN, Y_FLOW_MAX);
    TProfile *pV1_pt_proton_010 = new TProfile("pV1_pt_proton_010", "v_{1}(p_{T}) protons, 0-10%; p_{T} (GeV/c); v_{1}", N_PT_BINS, PT_MIN, PT_MAX);
    TProfile *pV2_pt_proton_010 = new TProfile("pV2_pt_proton_010", "v_{2}(p_{T}) protons, 0-10%; p_{T} (GeV/c); v_{2}", N_PT_BINS, PT_MIN, PT_MAX);

    TProfile *hYMean_v1_1040 = new TProfile("hYMean_v1_1040", "", N_Y_BINS, Y_FLOW_MIN, Y_FLOW_MAX);
    TProfile *hYMean_v2_1040 = new TProfile("hYMean_v2_1040", "", N_Y_BINS, Y_FLOW_MIN, Y_FLOW_MAX);
    TProfile *hYMean_v1_010 = new TProfile("hYMean_v1_010", "", N_Y_BINS, Y_FLOW_MIN, Y_FLOW_MAX);
    TProfile *hYMean_v2_010 = new TProfile("hYMean_v2_010", "", N_Y_BINS, Y_FLOW_MIN, Y_FLOW_MAX);

    TProfile *hPtMean_v1_1040 = new TProfile("hPtMean_v1_1040", "", N_PT_BINS, PT_MIN, PT_MAX);
    TProfile *hPtMean_v2_1040 = new TProfile("hPtMean_v2_1040", "", N_PT_BINS, PT_MIN, PT_MAX);
    TProfile *hPtMean_v1_010 = new TProfile("hPtMean_v1_010", "", N_PT_BINS, PT_MIN, PT_MAX);
    TProfile *hPtMean_v2_010 = new TProfile("hPtMean_v2_010", "", N_PT_BINS, PT_MIN, PT_MAX);

    TProfile *hBMean_v1 = new TProfile("hBMean_v1", "", 10, bEdges);
    TProfile *hBMean_v2 = new TProfile("hBMean_v2", "", 10, bEdges);
    TProfile *hCentMean_v1 = new TProfile("hCentMean_v1", "", 10, 0, 100);
    TProfile *hCentMean_v2 = new TProfile("hCentMean_v2", "", 10, 0, 100);

    // ==================== Основной цикл ====================
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        fChain->GetEntry(jentry);

        hB->Fill(fB);
        Int_t centBin = GetCentImpact_Unigen(fB);
        if (centBin < 0)
            continue;
        double centPercent = 0.5 * (minCentPercent[centBin] + maxCentPercent[centBin]);

        bool is10_40 = (centBin == 1 || centBin == 2 || centBin == 3);
        bool is0_10 = (centBin == 0);

        double psiRP = fPhi;
        double cosPhi = std::cos(psiRP);
        double sinPhi = std::sin(psiRP);

        int sumZb2 = 0, sumZb3 = 0, maxZ = 0, countZ2 = 0, countZ3 = 0, nFrag = 0;
        double sumZ2 = 0, sumA_frag = 0, sumZ_frag = 0;

        for (int i = 0; i < fParticles_; i++)
        {
            int pdg = fParticles_fPdg[i];
            double A, Z;
            GetAZfromPDG(pdg, A, Z);

            double px = fParticles_fPx[i];
            double py = fParticles_fPy[i];
            double pz = fParticles_fPz[i];
            double E = fParticles_fE[i];

            double mass = GetMassFromPDG(pdg, A);
            if (mass <= 0)
            {
                double p2tmp = px * px + py * py + pz * pz;
                mass = std::sqrt(std::max(0.0, E * E - p2tmp));
            }

            hA->Fill(A);
            hZ->Fill(Z);
            double p2 = px * px + py * py + pz * pz;
            double pmod = std::sqrt(p2);
            hP->Fill(pmod);

            double pxRP = px * cosPhi + py * sinPhi;
            double pyRP = -px * sinPhi + py * cosPhi;
            double pt = std::sqrt(pxRP * pxRP + pyRP * pyRP);
            hPt->Fill(pt);
            pPt_vs_b->Fill(fB, pt);

            double aAbs = std::fabs(A);
            double zAbs = std::fabs(Z);
            bool isProton = (std::fabs(aAbs - 1.0) < 0.5) && (zAbs > 0.5);

            if (E > std::fabs(pz) + 1e-9)
            {
                double y = 0.5 * std::log((E + pz) / (E - pz));
                hY->Fill(y);
                pPx_vs_y->Fill(y, pxRP);
                pPy_vs_y->Fill(y, pyRP);

                if (pt > 1e-6)
                {
                    double v1_val = pxRP / pt;
                    double v2_val = (pxRP * pxRP - pyRP * pyRP) / (pt * pt);

                    if (isProton)
                    {
                        double absY = std::fabs(y);
                        double v1_fold = (y >= 0.0) ? v1_val : -v1_val;
                        double v2_fold = v2_val;

                        bool inV1Y = (absY >= V1_Y_MIN && absY <= V1_Y_MAX);
                        bool inV2Y = (absY >= V2_Y_MIN && absY <= V2_Y_MAX);
                        bool inV1P = (pt >= V1_PT_MIN && pt <= V1_PT_MAX);
                        bool inV2P = (pt >= V2_PT_MIN && pt <= V2_PT_MAX);

                        if (is10_40 && inV1P)
                        {
                            pV1_y_proton->Fill(absY, v1_fold);
                            hYMean_v1_1040->Fill(absY, absY);
                        }
                        if (is10_40 && inV2P)
                        {
                            pV2_y_proton->Fill(absY, v2_fold);
                            hYMean_v2_1040->Fill(absY, absY);
                        }
                        if (is10_40 && inV1Y)
                        {
                            pV1_pt_proton->Fill(pt, v1_fold);
                            hPtMean_v1_1040->Fill(pt, pt);
                        }
                        if (is10_40 && inV2Y)
                        {
                            pV2_pt_proton->Fill(pt, v2_fold);
                            hPtMean_v2_1040->Fill(pt, pt);
                        }
                        if (inV1Y && inV1P)
                        {
                            pV1_cent_proton->Fill(centPercent, v1_fold);
                            hCentMean_v1->Fill(centPercent, centPercent);
                            pV1_b_proton->Fill(fB, v1_fold);
                            hBMean_v1->Fill(fB, fB);
                        }
                        if (inV2Y && inV2P)
                        {
                            pV2_cent_proton->Fill(centPercent, v2_fold);
                            hCentMean_v2->Fill(centPercent, centPercent);
                            pV2_b_proton->Fill(fB, v2_fold);
                            hBMean_v2->Fill(fB, fB);
                        }

                        if (is0_10 && inV1P)
                        {
                            pV1_y_proton_010->Fill(absY, v1_fold);
                            hYMean_v1_010->Fill(absY, absY);
                        }
                        if (is0_10 && inV2P)
                        {
                            pV2_y_proton_010->Fill(absY, v2_fold);
                            hYMean_v2_010->Fill(absY, absY);
                        }
                        if (is0_10 && inV1Y)
                        {
                            pV1_pt_proton_010->Fill(pt, v1_fold);
                            hPtMean_v1_010->Fill(pt, pt);
                        }
                        if (is0_10 && inV2Y)
                        {
                            pV2_pt_proton_010->Fill(pt, v2_fold);
                            hPtMean_v2_010->Fill(pt, pt);
                        }
                    }
                }
            }

            if (pmod > std::fabs(pz) + 1e-9)
                hEta->Fill(0.5 * std::log((pmod + pz) / (pmod - pz)));

            int Zi_i = (int)Z;
            if (isProton)
            {
                if (E > std::fabs(pz) + 1e-9)
                    hY_proton->Fill(0.5 * std::log((E + pz) / (E - pz)));
                if (pmod > std::fabs(pz) + 1e-9)
                    hEta_proton->Fill(0.5 * std::log((pmod + pz) / (pmod - pz)));
                hP_proton->Fill(pmod);
                hPt_proton->Fill(pt);
            }
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

            if (Zi_i >= 2)
            {
                sumZb2 += Zi_i;
                nFrag++;
                sumZ2 += Zi_i * Zi_i;
                sumA_frag += A;
                sumZ_frag += Zi_i;
                hA_frag->Fill(A);
                hZ_frag->Fill(Z);
                if (Zi_i > maxZ)
                    maxZ = Zi_i;
            }
            if (Zi_i >= 3)
                sumZb3 += Zi_i;
            if (Zi_i == 2)
                countZ2++;
            if (Zi_i == 3)
                countZ3++;
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

    // ==================== TGraphErrors со средним в бине ====================
    auto MakeMeanGraph = [](TProfile *prof, TProfile *hMean, const char *name)
    {
        TGraphErrors *g = new TGraphErrors();
        g->SetName(name);
        int n = 0;
        for (int i = 1; i <= prof->GetNbinsX(); ++i)
        {
            if (hMean->GetBinEntries(i) < 1)
                continue;
            g->SetPoint(n, hMean->GetBinContent(i), prof->GetBinContent(i));
            g->SetPointError(n, 0, prof->GetBinError(i));
            n++;
        }
        return g;
    };

    TGraphErrors *gV1_y_1040 = MakeMeanGraph(pV1_y_proton, hYMean_v1_1040, "gV1_y_1040");
    TGraphErrors *gV2_y_1040 = MakeMeanGraph(pV2_y_proton, hYMean_v2_1040, "gV2_y_1040");
    TGraphErrors *gV1_pt_1040 = MakeMeanGraph(pV1_pt_proton, hPtMean_v1_1040, "gV1_pt_1040");
    TGraphErrors *gV2_pt_1040 = MakeMeanGraph(pV2_pt_proton, hPtMean_v2_1040, "gV2_pt_1040");
    TGraphErrors *gV1_b_1040 = MakeMeanGraph(pV1_b_proton, hBMean_v1, "gV1_b_1040");
    TGraphErrors *gV2_b_1040 = MakeMeanGraph(pV2_b_proton, hBMean_v2, "gV2_b_1040");
    TGraphErrors *gV1_cent_1040 = MakeMeanGraph(pV1_cent_proton, hCentMean_v1, "gV1_cent_1040");
    TGraphErrors *gV2_cent_1040 = MakeMeanGraph(pV2_cent_proton, hCentMean_v2, "gV2_cent_1040");

    TGraphErrors *gV1_y_010 = MakeMeanGraph(pV1_y_proton_010, hYMean_v1_010, "gV1_y_010");
    TGraphErrors *gV2_y_010 = MakeMeanGraph(pV2_y_proton_010, hYMean_v2_010, "gV2_y_010");
    TGraphErrors *gV1_pt_010 = MakeMeanGraph(pV1_pt_proton_010, hPtMean_v1_010, "gV1_pt_010");
    TGraphErrors *gV2_pt_010 = MakeMeanGraph(pV2_pt_proton_010, hPtMean_v2_010, "gV2_pt_010");

    auto StyleMean = [](TGraphErrors *g, int color)
    {
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->SetLineWidth(2);
    };
    StyleMean(gV1_y_1040, kRed + 2);
    StyleMean(gV2_y_1040, kRed + 2);
    StyleMean(gV1_pt_1040, kRed + 2);
    StyleMean(gV2_pt_1040, kRed + 2);
    StyleMean(gV1_b_1040, kRed + 2);
    StyleMean(gV2_b_1040, kRed + 2);
    StyleMean(gV1_cent_1040, kRed + 2);
    StyleMean(gV2_cent_1040, kRed + 2);
    StyleMean(gV1_y_010, kRed + 2);
    StyleMean(gV2_y_010, kRed + 2);
    StyleMean(gV1_pt_010, kRed + 2);
    StyleMean(gV2_pt_010, kRed + 2);

    // ==================== Сохранение ROOT-файла ====================
    TFile *fout = new TFile(outputFileName, "RECREATE");
    hB->Write();
    hY->Write();
    hEta->Write();
    hPt->Write();
    hP->Write();
    hA->Write();
    hZ->Write();
    hA_frag->Write();
    hZ_frag->Write();
    hY_proton->Write();
    hEta_proton->Write();
    hP_proton->Write();
    hPt_proton->Write();
    hY_Z2->Write();
    hEta_Z2->Write();
    hP_Z2->Write();
    hPt_Z2->Write();
    hY_Z3->Write();
    hEta_Z3->Write();
    hP_Z3->Write();
    hPt_Z3->Write();
    pZmax_vs_Zb2->Write();
    pZmax_vs_Zb3->Write();
    pN2_vs_Zb2->Write();
    pN2_vs_Zb3->Write();
    pN3_vs_Zb2->Write();
    pN3_vs_Zb3->Write();
    pZ2_vs_Zb2->Write();
    pNfrag_vs_Zb2->Write();
    pA_vs_Zb2->Write();
    pZ_vs_Zb2->Write();
    pPx_vs_y->Write();
    pPy_vs_y->Write();
    pPt_vs_b->Write();
    pV1_y_proton->Write();
    pV2_y_proton->Write();
    pV1_pt_proton->Write();
    pV2_pt_proton->Write();
    pV1_cent_proton->Write();
    pV2_cent_proton->Write();
    pV1_b_proton->Write();
    pV2_b_proton->Write();
    pV1_y_proton_010->Write();
    pV2_y_proton_010->Write();
    pV1_pt_proton_010->Write();
    pV2_pt_proton_010->Write();
    hYMean_v1_1040->Write();
    hYMean_v2_1040->Write();
    hYMean_v1_010->Write();
    hYMean_v2_010->Write();
    hPtMean_v1_1040->Write();
    hPtMean_v2_1040->Write();
    hPtMean_v1_010->Write();
    hPtMean_v2_010->Write();
    hBMean_v1->Write();
    hBMean_v2->Write();
    hCentMean_v1->Write();
    hCentMean_v2->Write();
    gV1_y_1040->Write();
    gV2_y_1040->Write();
    gV1_pt_1040->Write();
    gV2_pt_1040->Write();
    gV1_b_1040->Write();
    gV2_b_1040->Write();
    gV1_cent_1040->Write();
    gV2_cent_1040->Write();
    gV1_y_010->Write();
    gV2_y_010->Write();
    gV1_pt_010->Write();
    gV2_pt_010->Write();
    fout->Close();

    std::cout << "Done. File saved: " << outputFileName << std::endl;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    // ==================== Канвас 1 ====================
    TCanvas *c1 = new TCanvas("c1", "Kinematics & impact parameter", 1800, 900);
    c1->Divide(4, 2);
    c1->cd(1);
    hB->Draw();
    c1->cd(2);
    hY->Draw();
    c1->cd(3);
    hEta->Draw();
    c1->cd(4);
    hPt->Draw();
    c1->cd(5);
    hP->Draw();
    c1->cd(6);
    pPx_vs_y->SetMarkerStyle(21);
    pPx_vs_y->SetMarkerColor(kRed + 2);
    pPx_vs_y->SetLineColor(kRed + 2);
    pPx_vs_y->Draw();
    {
        TLine *l = new TLine(-3, 0, 3, 0);
        l->SetLineStyle(2);
        l->Draw("same");
    }
    c1->cd(7);
    pPy_vs_y->SetMarkerStyle(22);
    pPy_vs_y->SetMarkerColor(kBlue + 2);
    pPy_vs_y->SetLineColor(kBlue + 2);
    pPy_vs_y->Draw();
    {
        TLine *l = new TLine(-3, 0, 3, 0);
        l->SetLineStyle(2);
        l->Draw("same");
    }
    c1->cd(8);
    pPt_vs_b->SetMarkerStyle(20);
    pPt_vs_b->SetMarkerColor(kBlack);
    pPt_vs_b->SetLineColor(kBlack);
    pPt_vs_b->Draw();
    c1->SaveAs(TString::Format("%s/kinematics.png", outDir.Data()));

    // ==================== Канвас 2 ====================
    TCanvas *c2 = new TCanvas("c2", "Fragmentation", 1600, 900);
    c2->Divide(3, 2);
    c2->cd(1);
    gPad->SetLogy();
    hA->SetLineColor(kBlack);
    hA_frag->SetLineColor(kRed + 2);
    hA_frag->SetLineWidth(2);
    hA->Draw();
    hA_frag->Draw("same");
    c2->cd(2);
    gPad->SetLogy();
    hZ->SetLineColor(kBlack);
    hZ_frag->SetLineColor(kRed + 2);
    hZ_frag->SetLineWidth(2);
    hZ->Draw();
    hZ_frag->Draw("same");
    c2->cd(3);
    pZmax_vs_Zb2->SetMarkerStyle(20);
    pZmax_vs_Zb2->SetMarkerColor(kRed + 2);
    pZmax_vs_Zb2->Draw();
    c2->cd(4);
    pZmax_vs_Zb3->SetMarkerStyle(21);
    pZmax_vs_Zb3->SetMarkerColor(kBlue + 2);
    pZmax_vs_Zb3->Draw();
    c2->cd(5);
    pN2_vs_Zb2->SetMarkerStyle(22);
    pN2_vs_Zb2->SetMarkerColor(kGreen + 2);
    pN2_vs_Zb2->Draw();
    c2->cd(6);
    pN2_vs_Zb3->SetMarkerStyle(23);
    pN2_vs_Zb3->SetMarkerColor(kMagenta);
    pN2_vs_Zb3->Draw();
    c2->SaveAs(TString::Format("%s/fragmentation.png", outDir.Data()));

    // ==================== Канвас 3 ====================
    TCanvas *c3 = new TCanvas("c3", "Fragmentation moments", 1500, 900);
    c3->Divide(3, 2);
    c3->cd(1);
    pN3_vs_Zb2->SetMarkerStyle(20);
    pN3_vs_Zb2->Draw();
    c3->cd(2);
    pN3_vs_Zb3->SetMarkerStyle(21);
    pN3_vs_Zb3->Draw();
    c3->cd(3);
    pNfrag_vs_Zb2->SetMarkerStyle(22);
    pNfrag_vs_Zb2->Draw();
    c3->cd(4);
    pZ2_vs_Zb2->SetMarkerStyle(23);
    pZ2_vs_Zb2->Draw();
    c3->cd(5);
    pA_vs_Zb2->SetMarkerStyle(24);
    pA_vs_Zb2->Draw();
    c3->cd(6);
    pZ_vs_Zb2->SetMarkerStyle(25);
    pZ_vs_Zb2->Draw();
    c3->SaveAs(TString::Format("%s/fragmentation_moments.png", outDir.Data()));

    // ==================== Канвас 4: Proton flow ====================
    TCanvas *c6 = new TCanvas("c6", "Proton flow", 1600, 1800);
    c6->Divide(2, 3);

    c6->cd(1);
    {
        gV1_y_010->Draw("APE1");
        gV2_y_010->Draw("PE1 SAME");
        gV1_y_010->SetTitle("v_{1}, v_{2}(|y|) protons, 0-10%; |y|; v_{1}, v_{2}");
        gV1_y_010->GetXaxis()->SetLimits(0, 1.2);
        TLine *l = new TLine(0, 0, 1.2, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.13, 0.15, 0.45, 0.35);
        leg->AddEntry(gV1_y_010, "v_{1}, 0-10%", "lp");
        leg->AddEntry(gV2_y_010, "v_{2}, 0-10%", "lp");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    c6->cd(2);
    {
        gV1_pt_010->Draw("APE1");
        gV2_pt_010->Draw("PE1 SAME");
        gV1_pt_010->SetTitle("v_{1}, v_{2}(p_{T}) protons, 0-10%; p_{T} (GeV/c); v_{1}, v_{2}");
        gV1_pt_010->GetXaxis()->SetLimits(0, 2.0);
        TLine *l = new TLine(0, 0, 2.0, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.55, 0.72, 0.92, 0.9);
        leg->AddEntry(gV1_pt_010, "v_{1}, 0-10%", "lp");
        leg->AddEntry(gV2_pt_010, "v_{2}, 0-10%", "lp");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    c6->cd(3);
    {
        gV1_y_1040->Draw("APE1");
        gV2_y_1040->Draw("PE1 SAME");
        gV1_y_1040->SetTitle("v_{1}, v_{2}(|y|) protons, 10-40%; |y|; v_{1}, v_{2}");
        gV1_y_1040->GetXaxis()->SetLimits(0, 1.2);
        TLine *l = new TLine(0, 0, 1.2, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.13, 0.15, 0.45, 0.35);
        leg->AddEntry(gV1_y_1040, "v_{1}, 10-40%", "lp");
        leg->AddEntry(gV2_y_1040, "v_{2}, 10-40%", "lp");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    c6->cd(4);
    {
        gV1_pt_1040->Draw("APE1");
        gV2_pt_1040->Draw("PE1 SAME");
        gV1_pt_1040->SetTitle("v_{1}, v_{2}(p_{T}) protons, 10-40%; p_{T} (GeV/c); v_{1}, v_{2}");
        gV1_pt_1040->GetXaxis()->SetLimits(0, 2.0);
        TLine *l = new TLine(0, 0, 2.0, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.55, 0.72, 0.92, 0.9);
        leg->AddEntry(gV1_pt_1040, "v_{1}, 10-40%", "lp");
        leg->AddEntry(gV2_pt_1040, "v_{2}, 10-40%", "lp");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    c6->cd(5);
    {
        gV1_cent_1040->Draw("APE1");
        gV2_cent_1040->Draw("PE1 SAME");
        gV1_cent_1040->SetTitle("v_{1}, v_{2} protons vs centrality; centrality (%); v_{1}, v_{2}");
        TLine *l = new TLine(0, 0, 100, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.55, 0.72, 0.92, 0.9);
        leg->AddEntry(gV1_cent_1040, "v_{1}", "lp");
        leg->AddEntry(gV2_cent_1040, "v_{2}", "lp");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    c6->cd(6);
    {
        gV1_b_1040->Draw("APE1");
        gV2_b_1040->Draw("PE1 SAME");
        gV1_b_1040->SetTitle("v_{1}, v_{2} protons vs <b>; <b> (fm); v_{1}, v_{2}");
        TLine *l = new TLine(0, 0, 12, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.55, 0.72, 0.92, 0.9);
        leg->AddEntry(gV1_b_1040, "v_{1}", "lp");
        leg->AddEntry(gV2_b_1040, "v_{2}", "lp");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    c6->SaveAs(TString::Format("%s/proton_flow.png", outDir.Data()));

    // ==================== Flow vs BM@N (10-40%) ====================
    TCanvas *cExp = new TCanvas("cExp", "Proton flow vs BM@N (10-40%)", 1800, 1200);
    cExp->Divide(3, 2);

    cExp->cd(1);
    {
        TGraphErrors *gBMn = MakeBMnGraph_Unigen(N_BMN_Y, bmn_y_v1_x, bmn_y_v1_y, "gBMn_v1_y");
        gV1_y_1040->SetTitle("v_{1}(y) protons, 10-40%; |y|; v_{1}");
        gV1_y_1040->SetMinimum(-0.05);
        gV1_y_1040->SetMaximum(0.45);
        gV1_y_1040->GetXaxis()->SetLimits(0, 1.2);
        gV1_y_1040->Draw("APE1");
        gBMn->Draw("P SAME");
        TLine *l = new TLine(0, 0, 1.2, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.13, 0.7, 0.55, 0.9);
        leg->AddEntry(gV1_y_1040, "UrQMD-AMC", "lp");
        leg->AddEntry(gBMn, "BM@N 10-40%", "p");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    cExp->cd(2);
    {
        TGraphErrors *gBMn = MakeBMnGraph_Unigen(N_BMN_PT, bmn_pt_v1_x, bmn_pt_v1_y, "gBMn_v1_pt");
        gV1_pt_1040->SetTitle("v_{1}(p_{T}) protons, 10-40%; p_{T} (GeV/c); v_{1}");
        gV1_pt_1040->SetMinimum(0.0);
        gV1_pt_1040->SetMaximum(0.45);
        gV1_pt_1040->GetXaxis()->SetLimits(0, 2.0);
        gV1_pt_1040->Draw("APE1");
        gBMn->Draw("P SAME");
        TLine *l = new TLine(0, 0, 2.0, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.55, 0.15, 0.92, 0.35);
        leg->AddEntry(gV1_pt_1040, "UrQMD-AMC", "lp");
        leg->AddEntry(gBMn, "BM@N 10-40%", "p");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    cExp->cd(3);
    {
        TGraphErrors *gBMn = MakeBMnGraph_Unigen(N_BMN_B, bmn_b_v1_x, bmn_b_v1_y, "gBMn_v1_b");
        gV1_b_1040->SetTitle("v_{1} protons vs <b>; <b> (fm); v_{1}");
        gV1_b_1040->SetMinimum(0.1);
        gV1_b_1040->SetMaximum(0.35);
        gV1_b_1040->Draw("APE1");
        gBMn->Draw("P SAME");
        TLine *l = new TLine(0, 0, 12, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.55, 0.7, 0.92, 0.9);
        leg->AddEntry(gV1_b_1040, "UrQMD-AMC", "lp");
        leg->AddEntry(gBMn, "BM@N 10-40%", "p");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    cExp->cd(4);
    {
        TGraphErrors *gBMn = MakeBMnGraph_Unigen(N_BMN_Y, bmn_y_v2_x, bmn_y_v2_y, "gBMn_v2_y");
        gV2_y_1040->SetMarkerStyle(24);
        gV2_y_1040->SetTitle("v_{2}(y) protons, 10-40%; |y|; v_{2}");
        gV2_y_1040->SetMinimum(-0.05);
        gV2_y_1040->SetMaximum(0.05);
        gV2_y_1040->GetXaxis()->SetLimits(0, 1.2);
        gV2_y_1040->Draw("APE1");
        gBMn->Draw("P SAME");
        TLine *l = new TLine(0, 0, 1.2, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.13, 0.15, 0.55, 0.35);
        leg->AddEntry(gV2_y_1040, "UrQMD-AMC", "lp");
        leg->AddEntry(gBMn, "BM@N 10-40%", "p");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    cExp->cd(5);
    {
        TGraphErrors *gBMn = MakeBMnGraph_Unigen(N_BMN_PT, bmn_pt_v2_x, bmn_pt_v2_y, "gBMn_v2_pt");
        gV2_pt_1040->SetMarkerStyle(24);
        gV2_pt_1040->SetTitle("v_{2}(p_{T}) protons, 10-40%; p_{T} (GeV/c); v_{2}");
        gV2_pt_1040->SetMinimum(-0.05);
        gV2_pt_1040->SetMaximum(0.05);
        gV2_pt_1040->GetXaxis()->SetLimits(0, 2.0);
        gV2_pt_1040->Draw("APE1");
        gBMn->Draw("P SAME");
        TLine *l = new TLine(0, 0, 2.0, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.55, 0.15, 0.92, 0.35);
        leg->AddEntry(gV2_pt_1040, "UrQMD-AMC", "lp");
        leg->AddEntry(gBMn, "BM@N 10-40%", "p");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    cExp->cd(6);
    {
        TGraphErrors *gBMn = MakeBMnGraph_Unigen(N_BMN_B, bmn_b_v2_x, bmn_b_v2_y, "gBMn_v2_b");
        gV2_b_1040->SetMarkerStyle(24);
        gV2_b_1040->SetTitle("v_{2} protons vs <b>; <b> (fm); v_{2}");
        gV2_b_1040->SetMinimum(-0.05);
        gV2_b_1040->SetMaximum(0.05);
        gV2_b_1040->Draw("APE1");
        gBMn->Draw("P SAME");
        TLine *l = new TLine(0, 0, 12, 0);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 1);
        l->Draw("same");
        TLegend *leg = new TLegend(0.55, 0.7, 0.92, 0.9);
        leg->AddEntry(gV2_b_1040, "UrQMD-AMC", "lp");
        leg->AddEntry(gBMn, "BM@N 10-40%", "p");
        StyleLegend_Unigen(leg);
        leg->Draw();
    }
    cExp->SaveAs(TString::Format("%s/flow_vs_BMn.png", outDir.Data()));

    // ==================== Канвас 5: Протоны ====================
    TCanvas *cProt = new TCanvas("cProt", "Protons: y, #eta, |p|, p_{T}", 1200, 900);
    cProt->Divide(2, 2);
    cProt->cd(1);
    hY_proton->Draw();
    cProt->cd(2);
    hEta_proton->Draw();
    cProt->cd(3);
    hP_proton->Draw();
    cProt->cd(4);
    hPt_proton->Draw();
    cProt->SaveAs(TString::Format("%s/protons_kinematics.png", outDir.Data()));

    // ==================== Канвас 6: Z = 2 ====================
    TCanvas *cZ2 = new TCanvas("cZ2", "Z=2: y, #eta, |p|, p_{T}", 1200, 900);
    cZ2->Divide(2, 2);
    cZ2->cd(1);
    hY_Z2->Draw();
    cZ2->cd(2);
    hEta_Z2->Draw();
    cZ2->cd(3);
    hP_Z2->Draw();
    cZ2->cd(4);
    hPt_Z2->Draw();
    cZ2->SaveAs(TString::Format("%s/Z2_kinematics.png", outDir.Data()));

    // ==================== Канвас 7: Z = 3 ====================
    TCanvas *cZ3 = new TCanvas("cZ3", "Z=3: y, #eta, |p|, p_{T}", 1200, 900);
    cZ3->Divide(2, 2);
    cZ3->cd(1);
    hY_Z3->Draw();
    cZ3->cd(2);
    hEta_Z3->Draw();
    cZ3->cd(3);
    hP_Z3->Draw();
    cZ3->cd(4);
    hPt_Z3->Draw();
    cZ3->SaveAs(TString::Format("%s/Z3_kinematics.png", outDir.Data()));
}