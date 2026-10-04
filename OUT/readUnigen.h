//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec 28 13:54:42 2025 by ROOT version 6.30/04
// from TTree events/AAMCC
// found on file: /home/oem/WORKLAST/unigen/urqmd_aamcc_xecs_3.26gev_cascade_12273365_109.unigen.root
//////////////////////////////////////////////////////////

#ifndef readUnigen_h
#define readUnigen_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

class readUnigen
{
public:
   TTree *fChain;  //! pointer to the analyzed TTree or TChain
   Int_t fCurrent; //! current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxfParticles = 410;
   static constexpr Int_t kMaxnucleons = 1;

   // Declaration of leaf types
   // UEvent          *event;
   UInt_t fUniqueID;
   UInt_t fBits;
   Int_t fEventNr;
   Double_t fB;
   Double_t fPhi;
   Int_t fNes;
   Int_t fStepNr;
   Double_t fStepT;
   Int_t fNpa;
   Int_t fParticles_;
   UInt_t fParticles_fUniqueID[kMaxfParticles];   //[fParticles_]
   UInt_t fParticles_fBits[kMaxfParticles];       //[fParticles_]
   Int_t fParticles_fIndex[kMaxfParticles];       //[fParticles_]
   Int_t fParticles_fPdg[kMaxfParticles];         //[fParticles_]
   Int_t fParticles_fStatus[kMaxfParticles];      //[fParticles_]
   Int_t fParticles_fParent[kMaxfParticles];      //[fParticles_]
   Int_t fParticles_fParentDecay[kMaxfParticles]; //[fParticles_]
   Int_t fParticles_fMate[kMaxfParticles];        //[fParticles_]
   Int_t fParticles_fDecay[kMaxfParticles];       //[fParticles_]
   Int_t fParticles_fChild[kMaxfParticles][2];    //[fParticles_]
   Double32_t fParticles_fPx[kMaxfParticles];     //[fParticles_]
   Double32_t fParticles_fPy[kMaxfParticles];     //[fParticles_]
   Double32_t fParticles_fPz[kMaxfParticles];     //[fParticles_]
   Double32_t fParticles_fE[kMaxfParticles];      //[fParticles_]
   Double32_t fParticles_fX[kMaxfParticles];      //[fParticles_]
   Double32_t fParticles_fY[kMaxfParticles];      //[fParticles_]
   Double32_t fParticles_fZ[kMaxfParticles];      //[fParticles_]
   Double32_t fParticles_fT[kMaxfParticles];      //[fParticles_]
   Double32_t fParticles_fWeight[kMaxfParticles]; //[fParticles_]
                                                  // EventInitialState *iniState;

   UShort_t id;
   UInt_t nColl;
   UInt_t nPart;
   Int_t nucleons_;
   UInt_t nucleons_fUniqueID[kMaxnucleons]; //[nucleons_]
   UInt_t nucleons_fBits[kMaxnucleons];     //[nucleons_]
   UShort_t nucleons_id[kMaxnucleons];      //[nucleons_]
   Int_t nucleons_pdgId[kMaxnucleons];      //[nucleons_]
   TLorentzVector nucleons_momentum[kMaxnucleons];
   TLorentzVector nucleons_position[kMaxnucleons];
   UShort_t nucleons_collisionType[kMaxnucleons]; //[nucleons_]
   vector<unsigned short> nucleons_collidedNucleonIndices[kMaxnucleons];

   // List of branches
   TBranch *b_event_fUniqueID;                 //!
   TBranch *b_event_fBits;                     //!
   TBranch *b_event_fEventNr;                  //!
   TBranch *b_event_fB;                        //!
   TBranch *b_event_fPhi;                      //!
   TBranch *b_event_fNes;                      //!
   TBranch *b_event_fStepNr;                   //!
   TBranch *b_event_fStepT;                    //!
   TBranch *b_event_fNpa;                      //!
   TBranch *b_event_fParticles_;               //!
   TBranch *b_fParticles_fUniqueID;            //!
   TBranch *b_fParticles_fBits;                //!
   TBranch *b_fParticles_fIndex;               //!
   TBranch *b_fParticles_fPdg;                 //!
   TBranch *b_fParticles_fStatus;              //!
   TBranch *b_fParticles_fParent;              //!
   TBranch *b_fParticles_fParentDecay;         //!
   TBranch *b_fParticles_fMate;                //!
   TBranch *b_fParticles_fDecay;               //!
   TBranch *b_fParticles_fChild;               //!
   TBranch *b_fParticles_fPx;                  //!
   TBranch *b_fParticles_fPy;                  //!
   TBranch *b_fParticles_fPz;                  //!
   TBranch *b_fParticles_fE;                   //!
   TBranch *b_fParticles_fX;                   //!
   TBranch *b_fParticles_fY;                   //!
   TBranch *b_fParticles_fZ;                   //!
   TBranch *b_fParticles_fT;                   //!
   TBranch *b_fParticles_fWeight;              //!
   TBranch *b_iniState_fUniqueID;              //!
   TBranch *b_iniState_fBits;                  //!
   TBranch *b_iniState_id;                     //!
   TBranch *b_iniState_nColl;                  //!
   TBranch *b_iniState_nPart;                  //!
   TBranch *b_iniState_nucleons_;              //!
   TBranch *b_nucleons_fUniqueID;              //!
   TBranch *b_nucleons_fBits;                  //!
   TBranch *b_nucleons_id;                     //!
   TBranch *b_nucleons_pdgId;                  //!
   TBranch *b_nucleons_momentum;               //!
   TBranch *b_nucleons_position;               //!
   TBranch *b_nucleons_collisionType;          //!
   TBranch *b_nucleons_collidedNucleonIndices; //!

   readUnigen(TTree *tree = 0);
   virtual ~readUnigen();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop(const char *inputFileName, const char *outputFileName);
   virtual Bool_t Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef readUnigen_cxx
readUnigen::readUnigen(TTree *tree) : fChain(0)
{
}

readUnigen::~readUnigen()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t readUnigen::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t readUnigen::LoadTree(Long64_t entry)
{
   // Set the environment to read one entry
   if (!fChain)
      return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0)
      return centry;
   if (fChain->GetTreeNumber() != fCurrent)
   {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void readUnigen::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   cout << endl
        << " init tree=" << tree;
   // Set branch addresses and branch pointers
   if (!tree)
      return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   cout << endl
        << " SetBranchAdress";
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_event_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_event_fBits);
   fChain->SetBranchAddress("fEventNr", &fEventNr, &b_event_fEventNr);
   fChain->SetBranchAddress("fB", &fB, &b_event_fB);
   fChain->SetBranchAddress("fPhi", &fPhi, &b_event_fPhi);
   fChain->SetBranchAddress("fNes", &fNes, &b_event_fNes);
   fChain->SetBranchAddress("fStepNr", &fStepNr, &b_event_fStepNr);
   fChain->SetBranchAddress("fStepT", &fStepT, &b_event_fStepT);
   fChain->SetBranchAddress("fNpa", &fNpa, &b_event_fNpa);
   fChain->SetBranchAddress("fParticles", &fParticles_, &b_event_fParticles_);
   fChain->SetBranchAddress("fParticles.fUniqueID", fParticles_fUniqueID, &b_fParticles_fUniqueID);
   fChain->SetBranchAddress("fParticles.fBits", fParticles_fBits, &b_fParticles_fBits);
   fChain->SetBranchAddress("fParticles.fIndex", fParticles_fIndex, &b_fParticles_fIndex);
   fChain->SetBranchAddress("fParticles.fPdg", fParticles_fPdg, &b_fParticles_fPdg);
   fChain->SetBranchAddress("fParticles.fStatus", fParticles_fStatus, &b_fParticles_fStatus);
   fChain->SetBranchAddress("fParticles.fParent", fParticles_fParent, &b_fParticles_fParent);
   fChain->SetBranchAddress("fParticles.fParentDecay", fParticles_fParentDecay, &b_fParticles_fParentDecay);
   fChain->SetBranchAddress("fParticles.fMate", fParticles_fMate, &b_fParticles_fMate);
   fChain->SetBranchAddress("fParticles.fDecay", fParticles_fDecay, &b_fParticles_fDecay);
   fChain->SetBranchAddress("fParticles.fChild[2]", fParticles_fChild, &b_fParticles_fChild);
   fChain->SetBranchAddress("fParticles.fPx", fParticles_fPx, &b_fParticles_fPx);
   fChain->SetBranchAddress("fParticles.fPy", fParticles_fPy, &b_fParticles_fPy);
   fChain->SetBranchAddress("fParticles.fPz", fParticles_fPz, &b_fParticles_fPz);
   fChain->SetBranchAddress("fParticles.fE", fParticles_fE, &b_fParticles_fE);
   fChain->SetBranchAddress("fParticles.fX", fParticles_fX, &b_fParticles_fX);
   fChain->SetBranchAddress("fParticles.fY", fParticles_fY, &b_fParticles_fY);
   fChain->SetBranchAddress("fParticles.fZ", fParticles_fZ, &b_fParticles_fZ);
   fChain->SetBranchAddress("fParticles.fT", fParticles_fT, &b_fParticles_fT);
   fChain->SetBranchAddress("fParticles.fWeight", fParticles_fWeight, &b_fParticles_fWeight);
   //    fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_iniState_fUniqueID);
   //    fChain->SetBranchAddress("fBits", &fBits, &b_iniState_fBits);
   fChain->SetBranchAddress("id", &id, &b_iniState_id);
   fChain->SetBranchAddress("nColl", &nColl, &b_iniState_nColl);
   fChain->SetBranchAddress("nPart", &nPart, &b_iniState_nPart);
   fChain->SetBranchAddress("nucleons", &nucleons_, &b_iniState_nucleons_);
   fChain->SetBranchAddress("nucleons.fUniqueID", &nucleons_fUniqueID, &b_nucleons_fUniqueID);
   fChain->SetBranchAddress("nucleons.fBits", &nucleons_fBits, &b_nucleons_fBits);
   fChain->SetBranchAddress("nucleons.id", &nucleons_id, &b_nucleons_id);
   fChain->SetBranchAddress("nucleons.pdgId", &nucleons_pdgId, &b_nucleons_pdgId);
   fChain->SetBranchAddress("nucleons.momentum", &nucleons_momentum, &b_nucleons_momentum);
   fChain->SetBranchAddress("nucleons.position", &nucleons_position, &b_nucleons_position);
   fChain->SetBranchAddress("nucleons.collisionType", &nucleons_collisionType, &b_nucleons_collisionType);
   fChain->SetBranchAddress("nucleons.collidedNucleonIndices", &nucleons_collidedNucleonIndices, &b_nucleons_collidedNucleonIndices);
   Notify();
}

Bool_t readUnigen::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void readUnigen::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t readUnigen::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}
#endif // #ifdef readUnigen_cxx
