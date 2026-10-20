// check_Z.C
R__ADD_INCLUDE_PATH(/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/include)
R__LOAD_LIBRARY(/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/build/libMcIniData.so)

#include "TFile.h"
#include "TTree.h"
#include "UEvent.h"
#include "UParticle.h"
#include <iostream>
#include <cmath>

void check_Z(const char* fn = "input.root", int nShow = 5)
{
    TFile* f = TFile::Open(fn);
    TTree* t = nullptr;
    f->GetObject("events", t);
    if (!t) { std::cout << "no events tree" << std::endl; return; }

    UEvent* ev = nullptr;
    t->SetBranchAddress("event", &ev);
    t->GetEntry(0);

    std::cout << "Npa = " << ev->GetNpa() << std::endl;

    double minx=1e9,maxx=-1e9, miny=1e9,maxy=-1e9, minz=1e9,maxz=-1e9;
    int n = 0;
    auto arr = ev->GetParticleList();
    for (int k = 0; k < ev->GetNpa(); ++k) {
        auto* p = (UParticle*)arr->At(k);
        if (!p) continue;
        if (p->GetStatus() != 0 || p->GetParent() != 0) continue;
        if (p->GetPdg() != 2212 && p->GetPdg() != 2112) continue;

        double x = p->X(), y = p->Y(), z = p->Z();
        minx = std::min(minx,x); maxx = std::max(maxx,x);
        miny = std::min(miny,y); maxy = std::max(maxy,y);
        minz = std::min(minz,z); maxz = std::max(maxz,z);
        if (n < nShow) {
            std::cout << "pdg=" << p->GetPdg()
                      << "  px=" << p->Px()
                      << "  py=" << p->Py()
                      << "  pz=" << p->Pz()
                      << "  x=" << x << "  y=" << y << "  z=" << z
                      << std::endl;
        }
        ++n;
    }
    std::cout << "nucleons = " << n << std::endl;
    std::cout << "X range: " << minx << " .. " << maxx << std::endl;
    std::cout << "Y range: " << miny << " .. " << maxy << std::endl;
    std::cout << "Z range: " << minz << " .. " << maxz << std::endl;
}