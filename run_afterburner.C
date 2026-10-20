#include "TString.h"
#include "TSystem.h"
#include <fstream>
#include <cstdio>
#include <cstdlib>

void run_afterburner(const char* inputFile = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/Input/xexe130.root",
                     double dmin = 2.6,
                     double V0 = 0.3,
                     int metricMode = 1,
                     
                     const char* outputBase = "/home/oem/WORKLAST/AAMCC-newDataHandle/OUT/out_XeXe_dmin2p6_mode1/xexe130",
                     int levelDensity = 4)
{
    if (!inputFile || !outputBase) {
        printf("Использование: run_afterburner(\"входной_файл\", dmin, V0, metricMode, \"базовое_имя_выхода_без_.root\")\n");
        printf("Пример (DCM-style): run_afterburner(\"/path/urqmd.root\", 2.0, 0.3, 1, \"XeXe_dcm\")\n");
        printf("Пример (pair-CM):   run_afterburner(\"/path/urqmd.root\", 2.0, 0.3, 2, \"XeXe_paircm\")\n");
        printf("Пример (spatial):   run_afterburner(\"/path/urqmd.root\", 2.1, 0,   0, 1, \"XeXe_spatial\")\n");
        printf("  dmin       — критическое расстояние (фм)\n");
        printf("  V0         — безразмерный порог скорости |dV| (V=P/E; <=0 → только пространственный MST)\n");
        printf("  metricMode — 0: spatial, 1: pair-CM (1 pass), 2: pair-CM (2 passes, constant V0)\n");
        printf("  В режиме 2 два прохода при одном V0; во 2-м dmin_eff = d0*(A_i+A_j)^(1/3)\n");
        return;
    }

    const char* GRATE_PATH = "/home/oem/WORKLAST/AAMCC-newDataHandle/build/GRATE";

    // ==== Создаём каталог для выходных файлов, если его нет ====
    // ROOT TFile::Open("RECREATE") не создаёт промежуточные каталоги,
    // поэтому выходной путь вида ".../OUT/подпапка/имя" без существующего
    // каталога приводил к падению. Создаём каталог заранее.
    {
        TString outPath(outputBase);
        Ssiz_t slashPos = outPath.Last('/');
        if (slashPos != kNPOS) {
            TString outDir = outPath(0, slashPos);
            if (!outDir.IsNull()) {
                gSystem->mkdir(outDir, kTRUE);
                if (gSystem->AccessPathName(outDir)) {
                    printf("Ошибка: не могу создать каталог %s\n", outDir.Data());
                    return;
                }
            }
        }
    }

    // ==== 8 ответов для режима afterburner ====
    // 1) читать Abrasion из файла
    // 2) генератор URQMD
    // 3) путь к входному файлу
    // 4) критическое расстояние (dmin) в фм
    // 5) MST V0 (безразмерный порог |dV|, <=0 → только пространственный MST)
    // 6) режим метрики MST (0 - spatial, 1 - pair-CM 1 pass, 2 - pair-CM 2 passes, constant V0)
    // 7) функция плотности уровней: 4 = Hybrid
    // 8) базовое имя выходного файла (без .root)
    TString answers;
    answers += "1\n";
    answers += "1\n";
    answers += TString(inputFile) + "\n";
    answers += TString::Format("%.4f\n", dmin);
    answers += TString::Format("%.4f\n", V0);
    answers += TString::Format("%d\n", metricMode);
    answers += TString::Format("%d\n", levelDensity);
    answers += TString(outputBase) + "\n";

    // Пишем ответы во временный файл
    TString tmpFile = "grate_answers.tmp";
    std::ofstream out(tmpFile.Data());
    if (!out.is_open()) {
        printf("Ошибка: не могу создать временный файл %s\n", tmpFile.Data());
        return;
    }
    out << answers;
    out.close();

    // Запускаем GRATE с абсолютным путём
    TString cmd = TString::Format("%s < %s", GRATE_PATH, tmpFile.Data());
    int ret = gSystem->Exec(cmd.Data());

    gSystem->Unlink(tmpFile.Data());

    if (ret == 0) {
        printf("Успешно. dmin = %.4f fm, V0 = %.4f (|dV|), metricMode = %d\n",
               dmin, V0, metricMode);
        printf("Созданы файлы: %s.root и %s_mcini_.root\n", outputBase, outputBase);
    } else {
        printf("Ошибка при выполнении GRATE. Код возврата: %d\n", ret);
    }
}
