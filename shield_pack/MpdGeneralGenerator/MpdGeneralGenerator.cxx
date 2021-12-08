// -------------------------------------------------------------------------
// -----          MpdGeneralGenerator source file                      -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

#include "MpdGeneralGenerator.h"

// // ------------------------------------------------------------------------
MpdGeneralGenerator::MpdGeneralGenerator(TGenerator *gen) : FairGenerator(), fGenerator(gen), fDebug(kFALSE){};
MpdGeneralGenerator::MpdGeneralGenerator(TGenerator &gen) : FairGenerator(), fGenerator(&gen), fDebug(kFALSE){};
MpdGeneralGenerator::MpdGeneralGenerator() : FairGenerator(), fGenerator(NULL), fDebug(kFALSE) {}
//     template <class T> MpdGeneralGenerator::MpdGeneralGenerator():
//     FairGenerator(), fDebug(kFALSE){
//         fGenerator = new T();
//     }
// ------------------------------------------------------------------------
Bool_t MpdGeneralGenerator::Init()
{
   // Initialize generator

   return kTRUE;
}

// ------------------------------------------------------------------------
Bool_t MpdGeneralGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{
   printf("MpdGeneralGenerator::ReadEvent\n");
   //     TClonesArray *p = new TClonesArray("TParticles");
   fGenerator->GenerateEvent();
   //     fGenerator->ImportParticles(p, "");
   TObjArray *p = fGenerator->ImportParticles();

   TParticle *par        = 0;
   int        numEntries = p->GetEntriesFast();
   for (Int_t k = 0; k < numEntries; k++) {
      par = (TParticle *)(p->UncheckedAt(k));
      if (fDebug)
         printf("GeneralGen: kf=%d, p=(%.2f, %.2f, %.2f) GeV, x=(%.1f, %.1f, %.1f) cm\n", par->GetPdgCode(), par->Px(),
                par->Py(), par->Pz(), par->Vx(), par->Vz(), par->Vz());
      if (par->GetPDG() == NULL) continue;
      primGen->AddTrack(par->GetPdgCode(), par->Px(), par->Py(), par->Pz(), par->Vx(), par->Vz(), par->Vz());
   }
   return kTRUE;
}
// ------------------------------------------------------------------------

ClassImp(MpdGeneralGenerator)
