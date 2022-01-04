// -------------------------------------------------------------------------
// -----                MpdSmashGenerator source file                  -----
// -----                Adopted  by I.Altsybeev 17/07/2020             -----
// -----                SMASH v.>=1.8 with ROOT output is required     -----
// -------------------------------------------------------------------------
#include "MpdSmashGenerator.h"

#include "FairPrimaryGenerator.h"
#include "FairMCEventHeader.h"

#include "TRandom.h"
#include "TMath.h"
#include "TString.h"
//#include "TTree.h"
//#include "TChain.h"
#include "TFile.h"

#include "TMCProcess.h"
#include "TObjArray.h"
#include "TPDGCode.h"
#include "TParticle.h"

#include "TVirtualMCStack.h"
#include "TLorentzVector.h"
//#include "TDatabasePDG.h"
//#include "TParticlePDG.h"

#include <iostream>
#include <cstring>

#include <stdio.h>

using namespace std;
using namespace TMath;

// -----   Default constructor   ------------------------------------------

MpdSmashGenerator::MpdSmashGenerator()
   : FairGenerator(), fInputFile(NULL), fFileName(""), fImpPar(0), fNpart(0), fEmptyEv(false)

     ,
     fIsRandomRP(false), fPsiRP(0), fRandom(0)
{
}
// ------------------------------------------------------------------------

// -----   Standard constructor   -----------------------------------------

MpdSmashGenerator::MpdSmashGenerator(TString fileName)
   : FairGenerator(), fInputFile(NULL), fFileName(fileName)

     ,
     fImpPar(0), fNpart(0), fEmptyEv(false)

     ,
     fIsRandomRP(false), fPsiRP(0), fRandom(0)
{
   cout << "-I MpdSmashGenerator: Opening input file " << fFileName << endl;
   fInputFile = new TFile(fFileName.Data());
   if (!fInputFile) {
      Fatal("MpdSmashGenerator", "Cannot open input file.");
      exit(1);
   }

   fInputTree = (TTree *)fInputFile->Get("particles");
   //    fInputTree = new TChain("out");
   //    fInputTree->Add(fFileName);

   cout << "-I MpdSmashGenerator: Number of entries in tree: " << fInputTree->GetEntries() << endl;

   fInputTree->SetBranchAddress("impact_b", &fImpPar);
   fInputTree->SetBranchAddress("npart", &fNpart);
   fInputTree->SetBranchAddress("empty_event", &fEmptyEv);
   fInputTree->SetBranchAddress("px", fPx);
   fInputTree->SetBranchAddress("py", fPy);
   fInputTree->SetBranchAddress("pz", fPz);
   fInputTree->SetBranchAddress("pdgcode", fPID);
   //    fInputTree->SetBranchAddress( "charge", fCharge );

   fEventNumber = 0;
   fPsiRP       = 0.;
   fIsRandomRP  = kTRUE; // by default RP is random

   fRandom = new TRandom3();
   fRandom->SetSeed(0);
   //    fPNCorr = 0.; // by default no correction
}
// ------------------------------------------------------------------------

// -----   Destructor   ---------------------------------------------------

MpdSmashGenerator::~MpdSmashGenerator()
{
   delete fInputFile;
   delete fRandom;
}
// ------------------------------------------------------------------------

// -----   Public method ReadEvent   --------------------------------------

Bool_t MpdSmashGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{
   // Check for input file
   if (!fInputFile) {
      cout << "-E MpdSmashGenerator: Input file not open! " << endl;
      return kFALSE;
   }

   // Check for primary generator
   if (!primGen) {
      cout << "-E- MpdSmashGenerator::ReadEvent: "
           << "No PrimaryGenerator!" << endl;
      return kFALSE;
   }

   // Skip empty SMASH events:
   //    bool hasInteraction = false;
   //    while( !hasInteraction )
   //    {
   //        fInputTree->GetEntry(fEventNumber);
   //        if ( fEmptyEv )
   //        {
   //            fEventNumber++;
   //            continue;
   //        }
   //        hasInteraction = true;
   //    }

   fInputTree->GetEntry(fEventNumber);

   cout << "-I MpdSmashGenerator: Event " << fEventNumber << ",  b = " << fImpPar << " fm,  multiplicity " << fNpart
        << endl; // << ", Elab: " << energ << endl;

   if (fNpart > MAX_N_PART) {
      cout << "-E- MpdSmashGenerator: ERROR! fNpart>MAX_N_PART (fNpart=" << fNpart << ")" << endl;
      exit(1);
   }

   /* set random Reaction Plane angle */
   fPsiRP = 0;
   if (fIsRandomRP) fPsiRP = fRandom->Uniform(2.0 * TMath::Pi());

   // Set event id and impact parameter in MCEvent if not yet done
   FairMCEventHeader *event = primGen->GetEvent();
   if (event && (!event->IsSet())) {
      // event->SetEventID(evnr);
      event->SetB(fImpPar);
      event->SetRotZ(fPsiRP);
      event->MarkSet(kTRUE);
   }

   // Loop over tracks in the current event
   for (int itrack = 0; itrack < fNpart; itrack++) {
      Double_t px  = fPx[itrack];
      Double_t py  = fPy[itrack];
      Double_t pz  = fPz[itrack];
      Int_t    pid = fPID[itrack];

      /* rotate RP angle */
      if (fPsiRP != 0.) {
         Double_t cosPsi = TMath::Cos(fPsiRP);
         Double_t sinPsi = TMath::Sin(fPsiRP);
         Double_t px1    = px * cosPsi - py * sinPsi;
         Double_t py1    = px * sinPsi + py * cosPsi;
         px              = px1;
         py              = py1;
      }

      // Give track to PrimaryGenerator
      primGen->AddTrack(pid, px, py, pz, 0.0, 0.0, 0.0); // pdg, mom [GeV/c], vertex [cm]
   }

   fEventNumber++;

   return kTRUE;
}
// ------------------------------------------------------------------------

ClassImp(MpdSmashGenerator);
