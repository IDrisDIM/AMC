#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include "AAMCConstants.hh"

#include "URun.h"

// AAMCCrun getTheRunData(const TString&);
AAMCCrun getTheRunData(const TString &inputFileName = "test.root", double critDist = 1.8);