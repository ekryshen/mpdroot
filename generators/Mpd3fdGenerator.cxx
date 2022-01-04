// -------------------------------------------------------------------------
// -----                Mpd3fdGenerator source file                  -----
// -----                Created 24/06/04  by V. Friese                 -----
// -------------------------------------------------------------------------
#include "Mpd3fdGenerator.h"

#include "FairPrimaryGenerator.h"
#include "FairMCEventHeader.h"

#include "TMCProcess.h"
#include "TObjArray.h"
#include "TPDGCode.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TString.h"
#include "TVirtualMCStack.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TMath.h"
#include "TString.h"
#include "TRegexp.h"

#include <iostream>
#include <cstring>

#include <stdio.h>

using namespace std;
using namespace TMath;

// -----   Default constructor   ------------------------------------------

Mpd3fdGenerator::Mpd3fdGenerator() : FairGenerator(), fInputFile(NULL), fFileName("") {}
// ------------------------------------------------------------------------

// -----   Standard constructor   -----------------------------------------

Mpd3fdGenerator::Mpd3fdGenerator(TString fileName) : FairGenerator(), fInputFile(NULL), fFileName(fileName)
{
   //  fFileName = fileName;
   cout << "-I Mpd3fdGenerator: Opening input file " << fFileName << endl;
   fInputFile = new TFile(fFileName.Data());
   if (!fInputFile) {
      Fatal("Mpd3fdGenerator", "Cannot open input file.");
      exit(1);
   }
   fDstTree = new TChain("out");
   fDstTree->Add(fFileName);

   // Deactivate All
   fDstTree->SetBranchStatus("*", 0);
   // Activate only selected branches
   fDstTree->SetBranchStatus("npart", 1);
   fDstTree->SetBranchStatus("px", 1);
   fDstTree->SetBranchStatus("py", 1);
   fDstTree->SetBranchStatus("pz", 1);
   // fDstTree->SetBranchStatus("x",1);
   // fDstTree->SetBranchStatus("y",1);
   // fDstTree->SetBranchStatus("z",1);
   fDstTree->SetBranchStatus("E", 1);
   fDstTree->SetBranchStatus("id", 1);

   fDstTree->SetBranchAddress("px", fPx);
   fDstTree->SetBranchAddress("py", fPy);
   fDstTree->SetBranchAddress("pz", fPz);
   // fDstTree->SetBranchAddress("x", fX);  // [fm]
   // fDstTree->SetBranchAddress("y", fY);  // [fm]
   // fDstTree->SetBranchAddress("z", fZ);  // [fm]
   fDstTree->SetBranchAddress("E", fE);
   fDstTree->SetBranchAddress("npart", &fNpart);
   fDstTree->SetBranchAddress("id", fPID);

   fEventNumber = 0;
   fPsiRP       = 0.;
   fisRP        = kTRUE; // by default RP is random
   frandom      = new TRandom2();
   frandom->SetSeed(0);
   fPNCorr = 0.; // by default no correction
}
// ------------------------------------------------------------------------

// -----   Destructor   ---------------------------------------------------

Mpd3fdGenerator::~Mpd3fdGenerator()
{
   delete fInputFile;
   delete fDstTree;
   delete frandom;
}
// ------------------------------------------------------------------------

// -----   Public method ReadEvent   --------------------------------------

Bool_t Mpd3fdGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{

   // ---> Check for input file
   if (!fInputFile) {
      cout << "-E Mpd3fdGenerator: Input file not open! " << endl;
      return kFALSE;
   }

   // ---> Check for primary generator
   if (!primGen) {
      cout << "-E- Mpd3fdGenerator::ReadEvent: "
           << "No PrimaryGenerator!" << endl;
      return kFALSE;
   }

   fDstTree->GetEntry(fEventNumber);

   // Int_t events = fDstTree->GetEntries();

   // ---> Define event variables to be read from file
   Int_t energ, impact;

   const TRegexp elb("_elb[0-9]*gev_");
   TSubString    subelb = fFileName(elb); // select substring by regexp "elab"
   sscanf(subelb.Data(), "_elb%dgev_", &energ);
   const TRegexp imp("_[0-9]*fm_");
   TSubString    subimp = fFileName(imp); // select substring by regexp "impact"
   sscanf(subelb.Data(), "_%dfm_", &impact);

   Float_t b = (Float_t)impact;

   /*
   // ---> Calculate beta and gamma for Lorentztransformation
   TDatabasePDG* pdgDB = TDatabasePDG::Instance();
   TParticlePDG* kProton = pdgDB->GetParticle(2212);
   Double_t kProtonMass = kProton->Mass();

   Double_t eBeam = energ + kProtonMass;
   Double_t pBeam = TMath::Sqrt(eBeam * eBeam - kProtonMass * kProtonMass);
   Double_t betaCM = pBeam / (eBeam + kProtonMass);
   Double_t gammaCM = TMath::Sqrt(1. / (1. - betaCM * betaCM));
   */

   cout << "-I Mpd3fdGenerator: Event " << fEventNumber << ",  b = " << b << " fm,  multiplicity " << fNpart
        << ", Elab: " << energ << endl;

   if (fNpart > kBatyukConst) {
      cout << "-E- 3fdGenerator SELFCHECK ERROR" << endl;
      exit(1);
   }

   /* set random Reaction Plane angle */
   if (fisRP) {
      fPsiRP = frandom->Uniform(2.0 * TMath::Pi());
   }

   // Set event id and impact parameter in MCEvent if not yet done
   FairMCEventHeader *event = primGen->GetEvent();
   if (event && (!event->IsSet())) {
      // event->SetEventID(evnr);
      event->SetB(b);
      event->SetRotZ(fPsiRP);
      event->MarkSet(kTRUE);
   }

   // ---> Loop over tracks in the current event
   for (int itrack = 0; itrack < fNpart; itrack++) {
      Double_t px = fPx[itrack];
      Double_t py = fPy[itrack];
      Double_t pz = fPz[itrack];
      /*
      // Lorentztransformation to lab
      Double_t e = fE[itrack];
      Double_t mass = Sqrt(e * e - px * px - py * py - pz * pz);
      if (gCoordinateSystem == sysLaboratory) { // this variable/enum is not even declared,
                                                // consider to remove this commented part, as it is dead code
        pz = gammaCM * (pz + betaCM * e);
      }
      Double_t ee = sqrt(mass * mass + px * px + py * py + pz * pz);
      */

      /* in Theseus 2018-03-17-bc2a06d <neutrons>=<protons> */
      Int_t pid = fPID[itrack];
      if (fPNCorr != 0. && (pid == 2212 || pid == 2112)) {
         Float_t xx = frandom->Rndm();
         if (xx > fPNCorr) {
            pid = 2112;
         } else {
            pid = 2212;
         };
      }
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

ClassImp(Mpd3fdGenerator);
