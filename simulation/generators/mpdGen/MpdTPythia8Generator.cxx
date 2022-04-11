/**
 * MpdTPythia8Generator.cxx
 * *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 * The MpdTPythia8Generator produces particles from TPythia8 output tree.
 **/

#include "MpdTPythia8Generator.h"
#include "FairMCEventHeader.h"
#include "FairPrimaryGenerator.h"

#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
//#include "TMCParticle.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVirtualMC.h"

#include <fstream>
#include <iostream>
#include <iomanip>
// #define debug_hsd

using namespace std;

// ------------------------------------------------------------------------
MpdTPythia8Generator::MpdTPythia8Generator()
   : FairGenerator(), fInputFile(NULL), fTree(NULL), fParticles(NULL), fNevents(0), fEdit(0)
{
   // Constructor
}

// ------------------------------------------------------------------------
MpdTPythia8Generator::MpdTPythia8Generator(TString fileName, TString treeName, TString branchName)
   : FairGenerator(), fParticles(new TClonesArray("TParticle", 100)),
     // fParticles(new TClonesArray("TMCParticle",100)),
     fNevents(0), fEdit(0)
{
   // Constructor

   // gSystem->Setenv("PYTHIA8DATA", "/home/zinchenk/pythia8/pythia8215/share/Pythia8/xmldoc");
   // TPythia8 pyt8; // to update PDG data base
   AddParticlesToPdgDataBase();

   fInputFile = new TFile(fileName);
   fTree      = (TTree *)fInputFile->Get(treeName.Data());
   fTree->SetBranchAddress(branchName.Data(), &fParticles);

   cout << "\n ***** MpdPythia8Generator ***** " << endl;
   cout << " Number of events in input file: " << fileName << " " << fTree->GetEntries() << endl;
   cout << " *****  ***** \n" << endl;
}

// ------------------------------------------------------------------------
MpdTPythia8Generator::~MpdTPythia8Generator()
{
   // Destructor

   fParticles->Delete();
   fInputFile->Close();
   delete fInputFile;
}

//-------------------------------------------------------------------------
Bool_t MpdTPythia8Generator::ReadEvent(FairPrimaryGenerator *primGen)
{
   // Generate event

   static Int_t first = 1;

   if (first) {
      first = 0;
      // cout << " histo " << fHist << " " << fHist->GetNbinsX() << " " << fHist->GetNbinsY() << endl;
      // if (TDatabasePDG::Instance()->GetParticle(fPdgCode) == NULL) {
      //  Define particle with fPdgCode (take it from DCM-QGSM input file)
      //  To be done ...
      //}
   }

   //
   fTree->GetEntry(fNevents++);

   Int_t nPart = fParticles->GetEntriesFast(), iFirstIndx = 0;
   cout << " Number of particles in event " << fNevents - 1 << ": " << nPart << endl;

   // TMCParticle *part = NULL;
   TParticle *part  = NULL;
   Bool_t     track = kTRUE;

   for (Int_t ipart = 0; ipart < nPart; ++ipart) {
      // part = (TMCParticle*) fParticles->UncheckedAt(ipart);
      part = (TParticle *)fParticles->UncheckedAt(ipart);
      // if (fEdit == 0 && part->GetKS() == 21) { ++iFirstIndx; continue; }

      // track = (part->GetKS() == 1) ? kTRUE : kFALSE;
      track = (part->GetStatusCode() == 1) ? kTRUE : kFALSE;

      // FairMCEventHeader* eventHeader = primGen->GetEvent();
      /* Set event impact parameter in MCEvent if not yet done */
      /*
        if (eventHeader && (!eventHeader->IsSet()) ) {
        eventHeader->SetB(fB);
        eventHeader->MarkSet(kTRUE);
        }
      */

      // Int_t parent = 0; //TMath::Max (part->GetParent() - 1 - iFirstIndx, -1);
      // Int_t parent = TMath::Max (part->GetFirstMother() - 1 - iFirstIndx, -1);
      Int_t parent = part->GetFirstMother();
      // cout << ipart << " " << part->GetStatusCode() << " " << part->GetPdgCode() << " " << part->GetFirstMother() <<
      // " " << parent << endl;

      // if (parent == -1) continue; // do not add beam protons - 5.12.2017
      primGen->AddTrack(part->GetPdgCode(), part->Px(), part->Py(), part->Pz(), part->Vx() * 0.1, part->Vy() * 0.1,
                        part->Vz() * 0.1, parent, track);
      // part->Vx(), part->Vy(), part->Vz(), -TMath::Abs(parent)-1, track); //AZ 5.12.2017
   }

   return kTRUE;
}
// ------------------------------------------------------------------------

void MpdTPythia8Generator::AddParticlesToPdgDataBase()
{
   // Add some pythia specific particle code to the data base - from TPythia8

   TDatabasePDG *pdgDB = TDatabasePDG::Instance();
   pdgDB->AddParticle("psi(3770)", "psi(3770)", 3.77315, kTRUE, 0, 0, "psi(3770)", 30443);
   pdgDB->AddParticle("f_0(980)", "f_0(980)", 1.00000, kTRUE, 0, 3, "f_0 meson", 9010221);
   pdgDB->AddParticle("Pomeron", "Pomeron", 0.00000, kTRUE, 0, 0, "Pomeron", 990);
   pdgDB->AddParticle("p_diffr+", "p_diffr+", 0.0, kTRUE, 0, 0, "p_diffr+", 9902210);
   pdgDB->AddParticle("J/psi[3S1(8)]", "J/psi[3S1(8)]", 3.29692, kTRUE, 0, 0, "Ccbar multiplet", 9940003);
   pdgDB->AddParticle("chi_2c[3S1(8)]", "chi_2c[3S1(8)]", 3.75620, kTRUE, 0, 0, "Ccbar multiplet", 9940005);
   pdgDB->AddParticle("chi_0c[3S1(8)]", "chi_0c[3S1(8)]", 3.61475, kTRUE, 0, 0, "Ccbar multiplet", 9940011);
   pdgDB->AddParticle("chi_1c[3S1(8)]", "chi_1c[3S1(8)]", 3.71066, kTRUE, 0, 0, "Ccbar multiplet", 9940023);
   pdgDB->AddParticle("psi(2S)[3S1(8)]", "psi(2S)[3S1(8)]", 3.88611, kTRUE, 0, 0, "QCD diffr. state", 9940103);
   pdgDB->AddParticle("J/psi[1S0(8)]", "J/psi[1S0(8)]", 3.29692, kTRUE, 0, 0, "Ccbar multiplet", 9941003);
   pdgDB->AddParticle("psi(2S)[1S0(8)]", "psi(2S)[1S0(8)]", 3.88611, kTRUE, 0, 0, "Ccbar multiplet", 9941103);
   pdgDB->AddParticle("J/psi[3PJ(8)]", "J/psi[3PJ(8)]", 3.29692, kTRUE, 0, 0, "Ccbar multiplet", 9942003);
   pdgDB->AddParticle("psi(3770)[3PJ(8)]", "psi(3770)[3PJ(8)]", 3.97315, kTRUE, 0, 0, "Ccbar multiplet", 9942033);
   pdgDB->AddParticle("psi(2S)[3PJ(8)]", "psi(2S)[3PJ(8)]", 3.88611, kTRUE, 0, 0, "Ccbar multiplet", 9942103);
   pdgDB->AddParticle("Upsilon[3S1(8)]", "Upsilon[3S1(8)]", 9.6603, kTRUE, 0, 0, "Bbbar multiplet", 9950003);
   pdgDB->AddParticle("chi_2b[3S1(8)]", "chi_2b[3S1(8)]", 10.1122, kTRUE, 0, 0, "Bbbar multiplet", 9950005);
   pdgDB->AddParticle("Upsilon[1S0(8)]", "Upsilon[1S0(8)]", 9.6603, kTRUE, 0, 0, "Bbbar multiplet", 9951003);
   pdgDB->AddParticle("Upsilon(2S)[1S0(8)]", "Upsilon(2S)[1S0(8)]", 10.22326, kTRUE, 0, 0, "Bbbar multiplet", 9951103);
   pdgDB->AddParticle("Upsilon(2S)[3PJ(8)]", "Upsilon(2S)[3PJ(8)]", 10.22326, kTRUE, 0, 0, "Bbbar multiplet", 9952103);
}

// ------------------------------------------------------------------------
