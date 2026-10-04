// ============================================================
// excitation_energy_comparison.C
//
// Строит зависимости энергии возбуждения префрагмента от его
// массового числа A для всех 7 параметризаций из ExcitationEnergy.cc:
//   1 - Ericson
//   2 - Gaimard-Schmidt
//   3 - ALADIN
//   4 - Damped Ericson (Hybrid of 1 and 3)
//   5 - Parabolic Approximation
//   6 - Corrected ALADIN
//   7 - Hybrid Fit
//
// Запуск:
//   root -l -q 'excitation_energy_comparison.C()'
//   root -l -q 'excitation_energy_comparison.C(124, 500)'
// ============================================================
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

static TRandom3* gRng = nullptr;
static double Gauss() { return gRng->Gaus(0, 1); }

// ============================================================
// 1) Ericson
// ============================================================
static double EricsonDist(double E, double Ebound, int a, int A) {
    const double g0 = 16.0;
    int R = A - a;
    if (R <= 0) return 0.0;
    if (E > Ebound * R) return 0.0;
    if (R < 2)  return (R == 1) ? g0 : 0.0;
    double s = g0;
    for (int k = 2; k <= R; ++k) s *= (g0 * E) / ((double)k * (k - 1));
    return s;
}

static double SampleEricson(int A, int initA, double Ebound = 40.0) {
    if (A >= initA) return 0.0;
    const int N = 10000;
    double Up = 100.0 * initA;
    double dx = Up / N;
    TH1F* h = new TH1F("h_tmp_er", "", N, 0, Up);
    for (int n = 0; n < N; ++n) {
        double E = (n + 0.5) * dx;
        h->SetBinContent(n + 1, EricsonDist(E, Ebound, A, initA));
    }
    double e = h->GetRandom();
    delete h;
    return e;
}

// ============================================================
// 2) Gaimard-Schmidt
// ============================================================
static double GaimardSchmidtDist(double E, double EvapEn, int a, int A) {
    const double g0 = 16.0, g1 = 0.7;
    int R = A - a;
    if (R < 2) return 0.0;
    if (E > 0.5 * EvapEn * R) return 0.0;

    double term = g0;
    for (int k = 2; k < R - 1; ++k) term *= (g0 * E) / ((double)k * (k - 1));
    double res = term;

    for (int m = 1; m < R; ++m) {
        term *= (g1 * E / g0) * (-(double)(R - m + 1)) / ((double)m * (R + m - 1));
        res += term;
    }
    if (res != res) res = 0.0;
    if (res < 0)    res = 0.0;
    return res;
}

static double SampleGaimardSchmidt(int A, int initA, double Ebound = 40.0) {
    if (A >= initA) return 0.0;
    const int N = 10000;
    double Up = 100.0 * initA;
    double dx = Up / N;
    TH1F* h = new TH1F("h_tmp_gs", "", N, 0, Up);
    for (int n = 0; n < N; ++n) {
        double E = (n + 0.5) * dx;
        h->SetBinContent(n + 1, GaimardSchmidtDist(E, Ebound, A, initA));
    }
    double e = h->GetRandom();
    delete h;
    return e;
}

// ============================================================
// 3) ALADIN
// ============================================================
static double SampleALADIN(int A, int initA, double e0, double sigma0, double b0) {
    if (A >= initA) return 0.0;
    double alpha = (double)A / initA;
    double energy = -1;
    int safety = 0;
    while ((energy < 0 || energy != energy) && safety++ < 10000) {
        double sigmaA = Gauss() * sigma0 * (1.0 + b0 * (1.0 - alpha));
        double arg = 1.0 - alpha - sigmaA;
        if (arg < 0) continue;
        energy = e0 * A * std::sqrt(arg);
    }
    return (energy < 0) ? 0.0 : energy;
}

// ============================================================
// alphaSwitch (используется в 4 и 7)
// ============================================================
static double ComputeAlphaSwitch(double e0, double Ebound, int initA) {
    double a = (e0 / Ebound) * (e0 / Ebound);
    for (int k = 0; k < 100; ++k) {
        double inner = (e0 / Ebound) *
                       (a * (1.0 - 1.0 / initA) + 1.0 / initA - a * a);
        if (inner < 0) inner = 0;
        a = std::pow(inner, 2.0 / 3.0);
    }
    return a;
}

// ============================================================
// 4) Damped Ericson (Hybrid of 1 and 3)
// ============================================================
static double SampleDampedEricson(int A, int initA, double alphaSwitch,
                                  double e0, double sigma0, double b0,
                                  double Ebound = 40.0) {
    if (initA - A < alphaSwitch * initA) return SampleEricson(A, initA, Ebound);
    else                                 return SampleALADIN(A, initA, e0, sigma0, b0);
}

// ============================================================
// 5) Parabolic Approximation
// ============================================================
static double SampleParabolic(int A, int initA, double Pe, double Pm,
                              double sigmaP, double bP0, double bP1) {
    if (A >= initA) return 0.0;
    double alpha = (double)A / initA;
    double sigmaE = Gauss() * sigmaP *
                    (1.0 + bP0 * (1.0 - alpha) + bP1 * (1.0 - alpha) * (1.0 - alpha));
    double e = Pe * A * (1.0 - alpha) * (alpha + Pm) + sigmaE;
    return (e < 0) ? 0.0 : e;
}

// ============================================================
// 6) Corrected ALADIN
// ============================================================
static double SampleCorrectedALADIN(int A, int initA,
                                    double e0, double c0, double sigma0,
                                    double b0, double b1) {
    if (A >= initA) return 0.0;
    double alpha = (double)A / initA;
    double energy = -1;
    int safety = 0;
    while ((energy < 0 || energy != energy) && safety++ < 10000) {
        double sigmaE = Gauss() * sigma0 *
                        (1.0 + b0 * (1.0 - alpha) + b1 * (1.0 - alpha) * (1.0 - alpha));
        double inner = std::sqrt(e0 * e0 + c0 * (1.0 - alpha)) - e0;
        if (inner < 0) inner = 0;
        energy = A * std::sqrt(inner) + A * sigmaE;
    }
    return (energy < 0) ? 0.0 : energy;
}

// ============================================================
// 7) Hybrid Fit
// ============================================================
static double SampleHybridFit(int A, int initA, double alphaSwitch,
                              double a0, double a1, double a2, double a3,
                              double a4, double a5, double a6,
                              double sigma1, double sigma2, double sigma3) {
    if (A >= initA) return 0.0;
    double alpha = (double)A / initA;
    double energy = -1;
    int safety = 0;
    while (energy < 0 && safety++ < 10000) {
        double r = Gauss();
        double sigmaE;
        if (alpha > alphaSwitch) {
            sigmaE = r * ((sigma3 - 1.0) / (1.0 - alphaSwitch) * alpha +
                          (1.0 - sigma3 * alphaSwitch) / (1.0 - alphaSwitch))
                     * sigma1 * (1.0 + sigma2 * (1.0 - alpha));
        } else {
            sigmaE = r * sigma1 * (1.0 + sigma2 * (1.0 - alpha));
        }
        double poly = a0 + a1 * alpha + a2 * std::pow(alpha, 2) +
                      a3 * std::pow(alpha, 3) + a4 * std::pow(alpha, 4) +
                      a5 * std::pow(alpha, 5) + a6 * std::pow(alpha, 6);
        energy = A * (poly + sigmaE);
    }
    return (energy < 0) ? 0.0 : energy;
}

// ============================================================
// Основная функция
// ============================================================
void excitation_energy_comparison(int initA = 124, int Nsamples = 300)
{
    if (!gRng) gRng = new TRandom3(42);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    // --- Параметры моделей (значения по умолчанию для тяжёлых ядер) ---
    const double e0_aladin = 11.5;      // MeV
    const double sigma0    = 0.005;
    const double b0_aladin = 2.0;
    const double c0_aladin = 0.1;
    const double Pe_par    = 24.0;      // MeV
    const double Pm_par    = 0.2;
    const double Ebound    = 40.0;

    // --- Hybrid Fit ---
    const double a0 = 11.46648905, a1 = -1.84830078, a2 = -58.53674677,
                 a3 = 284.66431513, a4 = -637.51406293, a5 = 652.80324427,
                 a6 = -251.28205381;
    const double sigma1 = 0.4, sigma2 = 0.5, sigma3 = 0.2;

    // --- Corrected ALADIN (значения из CorrectedALADINParameters.dat;
    //     подставьте свои, если знаете реальные) ---
    const double e0_corr = 10.0, c0_corr = 0.5;
    const double sigma0_corr = 0.05, b0_corr = 1.0, b1_corr = 0.1;

    double alphaSw = ComputeAlphaSwitch(e0_aladin, Ebound, initA);
    std::cout << "initA = " << initA
              << ", alphaSwitch = " << alphaSw
              << "  (Ericson при A > " << initA * (1.0 - alphaSw) << ")"
              << std::endl;

    const int NP = 7;
    const char* labels[NP] = {
        "1: Ericson",
        "2: Gaimard-Schmidt",
        "3: ALADIN",
        "4: Damped Ericson (Hybrid)",
        "5: Parabolic Approx",
        "6: Corrected ALADIN",
        "7: Hybrid Fit"
    };
    Color_t colors[NP] = {kRed+2, kBlue+2, kGreen+2, kMagenta+2,
                          kOrange+2, kCyan+2, kViolet+2};

    std::vector<int> Avec;
    for (int A = 2; A < initA; ++A) Avec.push_back(A);

    std::vector<TGraphErrors*> g_mean(NP), g_rms(NP);

    for (int p = 0; p < NP; ++p) {
        g_mean[p] = new TGraphErrors();
        g_rms[p]  = new TGraphErrors();
        g_mean[p]->SetName(TString::Format("g_mean_%d", p + 1));
        g_rms[p]->SetName(TString::Format("g_rms_%d",  p + 1));

        int n = 0;
        for (size_t iA = 0; iA < Avec.size(); ++iA) {
            int A = Avec[iA];
            double sum = 0, sum2 = 0;
            int nok = 0;
            for (int s = 0; s < Nsamples; ++s) {
                double e = 0;
                switch (p) {
                    case 0: e = SampleEricson(A, initA, Ebound); break;
                    case 1: e = SampleGaimardSchmidt(A, initA, Ebound); break;
                    case 2: e = SampleALADIN(A, initA, e0_aladin, sigma0, b0_aladin); break;
                    case 3: e = SampleDampedEricson(A, initA, alphaSw,
                                                    e0_aladin, sigma0, b0_aladin, Ebound); break;
                    case 4: e = SampleParabolic(A, initA, Pe_par, Pm_par,
                                                sigma0, c0_aladin, 0.01); break;
                    case 5: e = SampleCorrectedALADIN(A, initA,
                                                      e0_corr, c0_corr,
                                                      sigma0_corr, b0_corr, b1_corr); break;
                    case 6: e = SampleHybridFit(A, initA, alphaSw,
                                                a0, a1, a2, a3, a4, a5, a6,
                                                sigma1, sigma2, sigma3); break;
                }
                if (e != e || e < 0) continue;
                sum  += e;
                sum2 += e * e;
                nok++;
            }
            if (nok < 5) continue;

            double mean = sum / nok;
            double var  = sum2 / nok - mean * mean;
            if (var < 0) var = 0;
            double rms = std::sqrt(var);

            g_mean[p]->SetPoint(n, A, mean / A);
            g_mean[p]->SetPointError(n, 0, rms / std::sqrt((double)nok) / A);
            g_rms[p] ->SetPoint(n, A, rms / A);
            g_rms[p] ->SetPointError(n, 0, 0);
            ++n;
        }

        g_mean[p]->SetLineColor(colors[p]);
        g_mean[p]->SetMarkerColor(colors[p]);
        g_mean[p]->SetMarkerStyle(20 + p);
        g_mean[p]->SetMarkerSize(1.0);
        g_mean[p]->SetLineWidth(2);

        g_rms[p]->SetLineColor(colors[p]);
        g_rms[p]->SetMarkerColor(colors[p]);
        g_rms[p]->SetMarkerStyle(20 + p);
        g_rms[p]->SetMarkerSize(1.0);
        g_rms[p]->SetLineWidth(2);
    }

    // ==================== Canvas 1: <E*/A>(A) ====================
    {
        TCanvas* c1 = new TCanvas("c1", "<E*/A> vs A", 1200, 800);
        g_mean[0]->SetTitle(TString::Format(
            "<E*/A> vs A for ^{%d}Xe; A_{pf}; <E*/A> (MeV)", initA));
        g_mean[0]->GetXaxis()->SetTitleSize(0.045);
        g_mean[0]->GetYaxis()->SetTitleSize(0.045);
        g_mean[0]->GetXaxis()->SetTitleOffset(1.1);
        g_mean[0]->GetYaxis()->SetTitleOffset(1.3);

        double ymax = 0;
        for (int p = 0; p < NP; ++p)
            for (int i = 0; i < g_mean[p]->GetN(); ++i) {
                double x, y; g_mean[p]->GetPoint(i, x, y);
                if (y > ymax) ymax = y;
            }
        g_mean[0]->SetMinimum(0);
        g_mean[0]->SetMaximum(ymax * 1.15);

        g_mean[0]->Draw("APE1");
        for (int p = 1; p < NP; ++p) g_mean[p]->Draw("PE1 SAME");

        TLegend* leg = new TLegend(0.5, 0.5, 0.92, 0.92);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.8);
        leg->SetTextSize(0.032);
        for (int p = 0; p < NP; ++p) leg->AddEntry(g_mean[p], labels[p], "lp");
        leg->Draw();

        c1->SaveAs("excitation_energy_mean_vs_A.png");
    }

    // ==================== Canvas 2: RMS(E*/A)(A) ====================
    {
        TCanvas* c2 = new TCanvas("c2", "RMS(E*/A) vs A", 1200, 800);
        g_rms[0]->SetTitle(TString::Format(
            "RMS(E*/A) vs A for ^{%d}Xe; A_{pf}; RMS(E*/A) (MeV)", initA));
        g_rms[0]->GetXaxis()->SetTitleSize(0.045);
        g_rms[0]->GetYaxis()->SetTitleSize(0.045);
        g_rms[0]->GetXaxis()->SetTitleOffset(1.1);
        g_rms[0]->GetYaxis()->SetTitleOffset(1.3);

        double ymax = 0;
        for (int p = 0; p < NP; ++p)
            for (int i = 0; i < g_rms[p]->GetN(); ++i) {
                double x, y; g_rms[p]->GetPoint(i, x, y);
                if (y > ymax) ymax = y;
            }
        g_rms[0]->SetMinimum(0);
        g_rms[0]->SetMaximum(ymax * 1.15);

        g_rms[0]->Draw("APE1");
        for (int p = 1; p < NP; ++p) g_rms[p]->Draw("PE1 SAME");

        TLegend* leg = new TLegend(0.5, 0.5, 0.92, 0.92);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.8);
        leg->SetTextSize(0.032);
        for (int p = 0; p < NP; ++p) leg->AddEntry(g_rms[p], labels[p], "lp");
        leg->Draw();

        c2->SaveAs("excitation_energy_rms_vs_A.png");
    }

    // ==================== Canvas 3: <E*>(A) ====================
    {
        TCanvas* c3 = new TCanvas("c3", "<E*> vs A", 1200, 800);
        TGraphErrors* g_tot[NP];
        for (int p = 0; p < NP; ++p) {
            g_tot[p] = new TGraphErrors();
            g_tot[p]->SetName(TString::Format("g_tot_%d", p + 1));
            for (int i = 0; i < g_mean[p]->GetN(); ++i) {
                double x, y; g_mean[p]->GetPoint(i, x, y);
                g_tot[p]->SetPoint(i, x, y * x);
            }
            g_tot[p]->SetLineColor(colors[p]);
            g_tot[p]->SetMarkerColor(colors[p]);
            g_tot[p]->SetMarkerStyle(20 + p);
            g_tot[p]->SetMarkerSize(1.0);
            g_tot[p]->SetLineWidth(2);
        }
        g_tot[0]->SetTitle(TString::Format(
            "<E*> vs A for ^{%d}Xe; A_{pf}; <E*> (MeV)", initA));
        g_tot[0]->GetXaxis()->SetTitleSize(0.045);
        g_tot[0]->GetYaxis()->SetTitleSize(0.045);
        g_tot[0]->GetXaxis()->SetTitleOffset(1.1);
        g_tot[0]->GetYaxis()->SetTitleOffset(1.3);

        double ymax = 0;
        for (int p = 0; p < NP; ++p)
            for (int i = 0; i < g_tot[p]->GetN(); ++i) {
                double x, y; g_tot[p]->GetPoint(i, x, y);
                if (y > ymax) ymax = y;
            }
        g_tot[0]->SetMinimum(0);
        g_tot[0]->SetMaximum(ymax * 1.15);

        g_tot[0]->Draw("APE1");
        for (int p = 1; p < NP; ++p) g_tot[p]->Draw("PE1 SAME");

        TLegend* leg = new TLegend(0.5, 0.15, 0.92, 0.55);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.8);
        leg->SetTextSize(0.032);
        for (int p = 0; p < NP; ++p) leg->AddEntry(g_tot[p], labels[p], "lp");
        leg->Draw();

        c3->SaveAs("excitation_energy_total_vs_A.png");
    }

    std::cout << "\nСохранены файлы:" << std::endl;
    std::cout << "  excitation_energy_mean_vs_A.png"  << std::endl;
    std::cout << "  excitation_energy_rms_vs_A.png"   << std::endl;
    std::cout << "  excitation_energy_total_vs_A.png" << std::endl;
}

// ============================================================
// Comp() — сравнение моделей 1 (Ericson), 3 (ALADIN), 4 (Hybrid)
// в координатах alpha = A_pf / A_init.
//
// Строит 3 канваса:
//   excitation_alpha_mean.png           <E*/A> vs alpha
//   excitation_alpha_rms.png            RMS(E*/A) vs alpha
//   excitation_alpha_distributions.png  распределения E*/A при трёх alpha
// ============================================================
void Comp(int initA = 124, int Nsamples = 2000, int nAlphaBins = 40)
{
    if (!gRng) gRng = new TRandom3(42);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    // ---- Параметры моделей (те же, что в excitation_energy_comparison) ----
    const double e0_aladin = 11.5;      // MeV
    const double sigma0    = 0.005;
    const double b0_aladin = 2.0;
    const double Ebound    = 40.0;

    double alphaSwitch = ComputeAlphaSwitch(e0_aladin, Ebound, initA);
    double alphaSwitchLine = 1.0 - alphaSwitch;   // в терминах alpha = A/initA

    std::cout << "initA = " << initA
              << ", alphaSwitch (доля удалённых) = " << alphaSwitch
              << ", alpha переключения = " << alphaSwitchLine << std::endl;

    // ---- Три модели, которые сравниваем ----
    const int NM = 3;
    const char* modelNames[NM] = {"Ericson (1)", "ALADIN (3)", "Hybrid (4)"};
    Color_t modelColors[NM]   = {kRed+2, kGreen+2, kMagenta+2};
    Style_t modelMarkers[NM]  = {20, 21, 22};

    // ---- Сетка по alpha ----
    std::vector<double> alphas;
    for (int i = 0; i <= nAlphaBins; ++i) {
        double a = 0.10 + 0.88 * i / nAlphaBins;   // от 0.10 до 0.98
        alphas.push_back(a);
    }

    TGraphErrors* gMean[NM];
    TGraphErrors* gRms[NM];
    for (int m = 0; m < NM; ++m) {
        gMean[m] = new TGraphErrors();
        gRms[m]  = new TGraphErrors();
        gMean[m]->SetName(TString::Format("gMean_m%d", m));
        gRms[m] ->SetName(TString::Format("gRms_m%d",  m));
    }

    // ---- Заполнение графиков ----
    for (size_t ia = 0; ia < alphas.size(); ++ia) {
        double alpha = alphas[ia];
        int A = (int)std::round(alpha * initA);
        if (A < 2)      A = 2;
        if (A >= initA) A = initA - 1;

        double sum [NM] = {0, 0, 0};
        double sum2[NM] = {0, 0, 0};
        int    nok [NM] = {0, 0, 0};

        for (int s = 0; s < Nsamples; ++s) {
            double e[NM];
            e[0] = SampleEricson(A, initA, Ebound);
            e[1] = SampleALADIN(A, initA, e0_aladin, sigma0, b0_aladin);
            e[2] = SampleDampedEricson(A, initA, alphaSwitch,
                                       e0_aladin, sigma0, b0_aladin, Ebound);
            for (int m = 0; m < NM; ++m) {
                if (e[m] == e[m] && e[m] >= 0) {
                    sum[m]  += e[m];
                    sum2[m] += e[m] * e[m];
                    nok[m]++;
                }
            }
        }

        for (int m = 0; m < NM; ++m) {
            if (nok[m] < 5) continue;
            double mean = sum[m] / nok[m];
            double var  = sum2[m] / nok[m] - mean * mean;
            if (var < 0) var = 0;
            double rms = std::sqrt(var);

            int n = gMean[m]->GetN();
            gMean[m]->SetPoint(n, alpha, mean / A);
            gMean[m]->SetPointError(n, 0, rms / std::sqrt((double)nok[m]) / A);
            gRms[m]->SetPoint(n, alpha, rms / A);
            gRms[m]->SetPointError(n, 0, 0);
        }
    }

    // ---- Стиль ----
    for (int m = 0; m < NM; ++m) {
        for (TGraphErrors* g : {gMean[m], gRms[m]}) {
            g->SetLineColor(modelColors[m]);
            g->SetMarkerColor(modelColors[m]);
            g->SetMarkerStyle(modelMarkers[m]);
            g->SetMarkerSize(1.1);
            g->SetLineWidth(2);
        }
    }

    // ============================================================
    // CANVAS 1: <E*/A> vs alpha
    // ============================================================
    {
        TCanvas* c = new TCanvas("cMeanAlpha",
                                 Form("<E*/A> vs #alpha for ^{%d}Xe", initA),
                                 1200, 800);

        TGraphErrors* base = gMean[0];
        base->SetTitle(Form(
            "<E*/A> vs #alpha for ^{%d}Xe; #alpha = A_{pf}/A_{init}; <E*/A> (MeV)",
            initA));
        base->GetXaxis()->SetTitleSize(0.045);
        base->GetYaxis()->SetTitleSize(0.045);
        base->GetXaxis()->SetTitleOffset(1.1);
        base->GetYaxis()->SetTitleOffset(1.3);

        double ymax = 0;
        for (int m = 0; m < NM; ++m)
            for (int i = 0; i < gMean[m]->GetN(); ++i) {
                double x, y; gMean[m]->GetPoint(i, x, y);
                if (y > ymax) ymax = y;
            }

        base->GetXaxis()->SetLimits(0, 1);
        base->SetMinimum(0);
        base->SetMaximum(ymax * 1.15);

        base->Draw("APE1");
        for (int m = 1; m < NM; ++m) gMean[m]->Draw("PE1 SAME");

        TLine* lsw = new TLine(alphaSwitchLine, 0, alphaSwitchLine, ymax * 1.15);
        lsw->SetLineStyle(2);
        lsw->SetLineColor(kBlack);
        lsw->SetLineWidth(2);
        lsw->Draw("same");

        TLegend* leg = new TLegend(0.15, 0.68, 0.55, 0.9);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.8);
        leg->SetTextSize(0.032);
        for (int m = 0; m < NM; ++m)
            leg->AddEntry(gMean[m], modelNames[m], "lp");
        leg->AddEntry(lsw, Form("#alpha_{sw} = %.3f", alphaSwitchLine), "l");
        leg->Draw();

        c->SaveAs("excitation_alpha_mean.png");
    }

    // ============================================================
    // CANVAS 2: RMS(E*/A) vs alpha
    // ============================================================
    {
        TCanvas* c = new TCanvas("cRmsAlpha",
                                 Form("RMS(E*/A) vs #alpha for ^{%d}Xe", initA),
                                 1200, 800);

        TGraphErrors* base = gRms[0];
        base->SetTitle(Form(
            "RMS(E*/A) vs #alpha for ^{%d}Xe; #alpha = A_{pf}/A_{init}; RMS(E*/A) (MeV)",
            initA));
        base->GetXaxis()->SetTitleSize(0.045);
        base->GetYaxis()->SetTitleSize(0.045);
        base->GetXaxis()->SetTitleOffset(1.1);
        base->GetYaxis()->SetTitleOffset(1.3);

        double ymax = 0;
        for (int m = 0; m < NM; ++m)
            for (int i = 0; i < gRms[m]->GetN(); ++i) {
                double x, y; gRms[m]->GetPoint(i, x, y);
                if (y > ymax) ymax = y;
            }

        base->GetXaxis()->SetLimits(0, 1);
        base->SetMinimum(0);
        base->SetMaximum(ymax * 1.15);

        base->Draw("APE1");
        for (int m = 1; m < NM; ++m) gRms[m]->Draw("PE1 SAME");

        TLine* lsw = new TLine(alphaSwitchLine, 0, alphaSwitchLine, ymax * 1.15);
        lsw->SetLineStyle(2);
        lsw->SetLineColor(kBlack);
        lsw->SetLineWidth(2);
        lsw->Draw("same");

        TLegend* leg = new TLegend(0.15, 0.68, 0.55, 0.9);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(kWhite, 0.8);
        leg->SetTextSize(0.032);
        for (int m = 0; m < NM; ++m)
            leg->AddEntry(gRms[m], modelNames[m], "lp");
        leg->AddEntry(lsw, Form("#alpha_{sw} = %.3f", alphaSwitchLine), "l");
        leg->Draw();

        c->SaveAs("excitation_alpha_rms.png");
    }

    // ============================================================
    // CANVAS 3: Распределения E*/A при трёх значениях alpha
    // ============================================================
    {
        TCanvas* c = new TCanvas("cDist", "E*/A distributions", 1800, 650);
        c->Divide(3, 1);

        double alphasShow[3] = {0.95, alphaSwitchLine, 0.5};
        const char* padTitles[3] = {
            "Peripheral (#alpha = 0.95)",
            Form("Transition (#alpha_{sw} = %.3f)", alphaSwitchLine),
            "Central (#alpha = 0.5)"
        };

        const int    NBINS = 80;
        const double EMAX  = 15.0;
        const int    NEV   = 20000;

        for (int pad = 0; pad < 3; ++pad) {
            c->cd(pad + 1);
            double alpha = alphasShow[pad];
            int A = (int)std::round(alpha * initA);
            if (A < 2)      A = 2;
            if (A >= initA) A = initA - 1;

            TH1F* hists[NM];
            for (int m = 0; m < NM; ++m) {
                hists[m] = new TH1F(
                    Form("hDist_p%d_m%d", pad, m),
                    Form("%s; E*/A (MeV); Probability / bin", padTitles[pad]),
                    NBINS, 0, EMAX);
                hists[m]->SetStats(0);
            }

            for (int s = 0; s < NEV; ++s) {
                double e[NM];
                e[0] = SampleEricson(A, initA, Ebound);
                e[1] = SampleALADIN(A, initA, e0_aladin, sigma0, b0_aladin);
                e[2] = SampleDampedEricson(A, initA, alphaSwitch,
                                           e0_aladin, sigma0, b0_aladin, Ebound);
                for (int m = 0; m < NM; ++m)
                    if (e[m] == e[m] && e[m] >= 0)
                        hists[m]->Fill(e[m] / A);
            }

            double ymax = 0;
            for (int m = 0; m < NM; ++m) {
                if (hists[m]->Integral() > 0)
                    hists[m]->Scale(1.0 / hists[m]->Integral());
                hists[m]->SetLineColor(modelColors[m]);
                hists[m]->SetLineWidth(2);
                if (hists[m]->GetMaximum() > ymax)
                    ymax = hists[m]->GetMaximum();
            }

            for (int m = 0; m < NM; ++m) {
                hists[m]->SetMaximum(ymax * 1.15);
                if (m == 0) hists[m]->Draw("HIST");
                else        hists[m]->Draw("HIST SAME");
            }

            TLegend* leg = new TLegend(0.5, 0.65, 0.92, 0.9);
            leg->SetBorderSize(0);
            leg->SetFillColorAlpha(kWhite, 0.8);
            leg->SetTextSize(0.04);
            for (int m = 0; m < NM; ++m)
                leg->AddEntry(hists[m], modelNames[m], "l");
            leg->Draw();

            TLatex* lat = new TLatex();
            lat->SetNDC();
            lat->SetTextSize(0.045);
            lat->DrawLatex(0.15, 0.85, Form("A = %d, #alpha = %.3f", A, alpha));
        }
        c->SaveAs("excitation_alpha_distributions.png");
    }

    std::cout << "\nComp() сохранил:" << std::endl;
    std::cout << "  excitation_alpha_mean.png"          << std::endl;
    std::cout << "  excitation_alpha_rms.png"           << std::endl;
    std::cout << "  excitation_alpha_distributions.png" << std::endl;
}