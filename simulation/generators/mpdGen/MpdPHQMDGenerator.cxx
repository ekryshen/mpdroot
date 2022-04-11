/**
 *@class MpdPHQMDGenerator
 *@author V.Kireyeu, V. Voronyuk
 **/

#include <stdio.h>
#include <stdlib.h>

#include <TMath.h>
#include <TMCParticleType.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TLorentzVector.h>

#include "FairMCEventHeader.h"
#include "MpdMCEventHeader.h"
#include "MpdPHQMDGenerator.h"

// ---------------------------------------------------------------------
MpdPHQMDGenerator::MpdPHQMDGenerator(){
   /* It is better to leave empty */
};
// ---------------------------------------------------------------------
MpdPHQMDGenerator::MpdPHQMDGenerator(const char *filename = "phsd.dat.gz", const char *f79name = "fort.79.gz")
{
   fgzFile   = gzopen(filename, "rb"); // zlib
   fgzFile79 = gzopen(f79name, "rb");  // zlib
   if (!fgzFile) {
      printf("-E- MpdPHQMDGenerator: can not open file: %s\n", filename);
      exit(1);
   }
   printf("-I- MpdPHQMDGenerator: open %s\n", filename);

   if (!fgzFile79) {
      printf("-E- MpdPHQMDGenerator: can not open file: %s\n", f79name);
      exit(1);
   }
   printf("-I- MpdPHQMDGenerator: open %s\n", f79name);

   fPsiRP  = 0.;
   fisRP   = kTRUE; // by default RP is random
   frandom = new TRandom2();
   frandom->SetSeed(0);

   TDatabasePDG *fPDG     = TDatabasePDG::Instance();
   TParticlePDG *particle = 0;
   kProtonMass            = fPDG->GetParticle(2212)->Mass();
   kNeutronMass           = fPDG->GetParticle(2112)->Mass();
   kLambdaMass            = fPDG->GetParticle(3122)->Mass();

   // Light nuclei
   particle = fPDG->AddParticle("d", "d", 1.8761199, 0, 1.6916047e-33, 1 * 3, "nucleus", 1000010020, 0, kPTHadron);
   particle = fPDG->AddParticle("t", "t", 2.8094298, 0, 1.6916047e-33, 1 * 3, "nucleus", 1000010030, 0, kPTHadron);
   particle = fPDG->AddParticle("He3", "He3", 2.8094101, 0, 1.6916047e-33, 2 * 3, "nucleus", 1000020030, 0, kPTHadron);
   particle = fPDG->AddParticle("He4", "He4", 3.7283999, 0, 1.6916047e-33, 2 * 3, "nucleus", 1000020040, 0, kPTHadron);
   particle = fPDG->AddParticle("He5", "He5", 5.01222363, 0, 1.6916047e-33, 2 * 3, "nucleus", 1000020050, 0, kPTHadron);
   particle = fPDG->AddParticle("He6", "He6", 5.6056, 0, 1.6916047e-33, 2 * 3, "nucleus", 1000020060, 0, kPTHadron);
   particle = fPDG->AddParticle("Li6", "Li6", 5.6016, 0, 1.6916047e-33, 3 * 3, "nucleus", 1000030060, 0, kPTHadron);
   particle = fPDG->AddParticle("Li7", "Li7", 6.5339, 0, 1.6916047e-33, 3 * 3, "nucleus", 1000030070, 0, kPTHadron);
   particle = fPDG->AddParticle("Li8", "Li8", 7.4714, 0, 1.6916047e-33, 3 * 3, "nucleus", 1000030080, 0, kPTHadron);
   particle = fPDG->AddParticle("Be7", "Be7", 6.5342, 0, 1.6916047e-33, 4 * 3, "nucleus", 1000040070, 0, kPTHadron);
   particle = fPDG->AddParticle("Be8", "Be8", 7.4569, 0, 1.6916047e-33, 4 * 3, "nucleus", 1000040080, 0, kPTHadron);
   particle = fPDG->AddParticle("Be9", "Be9", 8.3959, 0, 1.6916047e-33, 4 * 3, "nucleus", 1000040090, 0, kPTHadron);
   particle = fPDG->AddParticle("B8", "B8", 7.4724, 0, 1.6916047e-33, 5 * 3, "nucleus", 1000050080, 0, kPTHadron);
   particle = fPDG->AddParticle("B9", "B9", 8.3959, 0, 1.6916047e-33, 5 * 3, "nucleus", 1000050090, 0, kPTHadron);

   // Hypernuclei desribed in UserDecay.C
};

// ---------------------------------------------------------------------
MpdPHQMDGenerator::~MpdPHQMDGenerator()
{
   if (fgzFile) {
      gzclose(fgzFile);
      fgzFile = NULL;
   }
   if (fgzFile79) {
      gzclose(fgzFile79);
      fgzFile79 = NULL;
   }
   delete frandom;
};
// ---------------------------------------------------------------------
Bool_t MpdPHQMDGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{
   if (!ReadHeader()) {
      printf("-I- MpdPHQMDGenerator: End of file reached.\n");
      return kFALSE; // end of file
   }

   /* set random Reaction Plane angle */
   if (fisRP) {
      fPsiRP = frandom->Uniform(2.0 * TMath::Pi());
   }
   // printf("-I- MpdPHQMDGenerator: RP angle = %e\n",fPsiRP);

   /* Set event impact parameter in MCEventHeader if not yet done */
   FairMCEventHeader *eventHeader = primGen->GetEvent();
   if (eventHeader && (!eventHeader->IsSet())) {
      eventHeader->SetB(fb);
      eventHeader->MarkSet(kTRUE);
      eventHeader->SetRotZ(fPsiRP);
   }

   Int_t   ipdg;
   Float_t px, py, pz;
   Float_t mass, amass;
   Int_t   nZ, nN, nL, nS;

   /* read tracks: fort.79 - baryons + clusters */
   for (Int_t i = 1; i <= fntr79; i++) {
      if (gzgets(fgzFile79, fbuffer79, 256) == Z_NULL) {
         printf("-E- MpdPHQMDGenerator: unexpected end of baryons file %d/%d\n", i, fntr79);
         exit(1);
      }
      int res79 = sscanf(fbuffer79, "%*d %d %d %e %e %e %*e %*e %*e %*e %e %e %d %d", &nZ, &nN, &px, &py, &pz, &mass,
                         &amass, &nL, &nS);
      if (res79 != 9) {
         printf("-E- MpdPHQMDGenerator: selftest error in track, scan %d of 9\n", res79);
         exit(1);
      }
      ipdg = BaryonPDG(nN, nZ, nL, nS, amass);

      Float_t        pAbs = 0.0, tmpZ, st, phi;
      TVector3       velocity;
      TLorentzVector pMother, tmpMom;
      if (nZ > 1 && nN == 0 && nL == 0 && nS == 0) { // bounded protons case
         pAbs = 0.5 * sqrt(4 * pow(kProtonMass, 2 * nZ)) / amass;
         pMother.SetPxPyPzE(px, py, pz, sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2) + pow(amass, 2)));

         for (int bound = 0; bound < nZ; ++bound) {
            tmpZ = 1. - 2. * (gRandom->Rndm());
            st   = sqrt((1. - tmpZ) * (1. + tmpZ)) * pAbs;
            phi  = 2. * TMath::Pi() * (gRandom->Rndm());

            tmpMom.SetXYZM(st * cos(phi), st * sin(phi), tmpZ * pAbs, kProtonMass);
            tmpMom.Boost(velocity);
            primGen->AddTrack(2212, tmpMom.Px(), tmpMom.Py(), tmpMom.Pz(), 0., 0., 0.);
         }

      } else if (nZ == 0 && nN > 1 && nL == 0 && nS == 0) { // bounded neutrons case
         pAbs = 0.5 * sqrt(4 * pow(kNeutronMass, 2 * nN)) / amass;
         pMother.SetPxPyPzE(px, py, pz, sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2) + pow(amass, 2)));

         for (int bound = 0; bound < nN; ++bound) {
            tmpZ = 1. - 2. * (gRandom->Rndm());
            st   = sqrt((1. - tmpZ) * (1. + tmpZ)) * pAbs;
            phi  = 2. * TMath::Pi() * (gRandom->Rndm());

            tmpMom.SetXYZM(st * cos(phi), st * sin(phi), tmpZ * pAbs, kNeutronMass);
            tmpMom.Boost(velocity);
            primGen->AddTrack(2112, tmpMom.Px(), tmpMom.Py(), tmpMom.Pz(), 0., 0., 0.);
         }
      } else if (nZ == 0 && nN > 1 && nL >= 1 && nS == 0) { // bounded lambdas + neutrons case
         pAbs = 0.5 * sqrt(4 * pow(kNeutronMass, 2 * nN) * pow(kLambdaMass, 2 * nL)) / amass;
         pMother.SetPxPyPzE(px, py, pz, sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2) + pow(amass, 2)));

         for (int bound = 0; bound < nN; ++bound) { // free neutrons
            tmpZ = 1. - 2. * (gRandom->Rndm());
            st   = sqrt((1. - tmpZ) * (1. + tmpZ)) * pAbs;
            phi  = 2. * TMath::Pi() * (gRandom->Rndm());

            tmpMom.SetXYZM(st * cos(phi), st * sin(phi), tmpZ * pAbs, kNeutronMass);
            tmpMom.Boost(velocity);
            primGen->AddTrack(2112, tmpMom.Px(), tmpMom.Py(), tmpMom.Pz(), 0., 0., 0.);
         }
         for (int bound = 0; bound < nL; ++bound) { // free lambdas
            tmpZ = 1. - 2. * (gRandom->Rndm());
            st   = sqrt((1. - tmpZ) * (1. + tmpZ)) * pAbs;
            phi  = 2. * TMath::Pi() * (gRandom->Rndm());

            tmpMom.SetXYZM(st * cos(phi), st * sin(phi), tmpZ * pAbs, kLambdaMass);
            tmpMom.Boost(velocity);
            primGen->AddTrack(3122, tmpMom.Px(), tmpMom.Py(), tmpMom.Pz(), 0., 0., 0.);
         }
      } else { // normal case
         primGen->AddTrack(ipdg, px, py, pz, 0., 0., 0.);
      }
   }

   /* read tracks: phsd.dat */
   for (Int_t i = 1; i <= fntr; i++) {
      /* read track */
      if (gzgets(fgzFile, fbuffer, 256) == Z_NULL) {
         printf("-E- MpdPHQMDGenerator: unexpected end of file %d/%d\n", i, fntr);
         exit(1);
      }
      /* scan values */
      int res = sscanf(fbuffer, "%d %*d %e %e %e", &ipdg, &px, &py, &pz);
      if (res != 4) {
         printf("-E- MpdPHQMDGenerator: selftest error in track, scan %d of 4\n", res);
         exit(1);
      }
      if (ipdg == 0) {
         printf("-W- MpdPHQMDGenerator: particle with pdg=0\n");
         continue;
      }
      if (ipdg > 999) continue; // only mesons should be readed

      // replacement of heavy particles for Geant4
      // this decays will be improved in the next versions of PHQMD
      if (ipdg == +413) ipdg = +421;   // D*+     -> D0
      if (ipdg == -413) ipdg = -421;   // D*-     -> D0_bar
      if (ipdg == +423) ipdg = +421;   // D*0     -> D0
      if (ipdg == -423) ipdg = -421;   // D*0_bar -> D0_bar
      if (ipdg == +4214) ipdg = +4122; // Sigma*_c+  -> Lambda_c+
      if (ipdg == -4214) ipdg = -4122; // Sigma*_c-  -> Lambda_c-
      if (ipdg == +4114) ipdg = +4122; // Sigma*_c0  -> Lambda_c+
      if (ipdg == -4114) ipdg = -4122; // Sigma*_c0  -> Lambda_c-
      if (ipdg == +4224) ipdg = +4122; // Sigma*_c++ -> Lambda_c+
      if (ipdg == -4224) ipdg = -4122; // Sigma*_c-- -> Lambda_c-

      /* rotate RP angle */
      if (fPsiRP != 0.) {
         Double_t cosPsi = TMath::Cos(fPsiRP);
         Double_t sinPsi = TMath::Sin(fPsiRP);
         Double_t px1    = px * cosPsi - py * sinPsi;
         Double_t py1    = px * sinPsi + py * cosPsi;
         px              = px1;
         py              = py1;
      }

      /* add track to simulation */
      primGen->AddTrack(ipdg, px, py, pz, 0., 0., 0.);
   }

   return kTRUE;
};
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
Int_t MpdPHQMDGenerator::BaryonPDG(Int_t N, Int_t Z, Int_t L, Int_t S, Float_t mass)
{
   Int_t I   = 0;
   Int_t A   = N + Z + L + S;
   Int_t LS  = L + S;
   Int_t PDG = 999; // some default value

   if (N == 1 && Z == 0 && L == 0 && S == 0)
      PDG = 2112; // neutron
   else if (N == 0 && Z == 1 && L == 0 && S == 0)
      PDG = 2212; // proton
   else if (N == 0 && Z == 0 && L == 1 && S == 0)
      PDG = 3122; // lambda0
   else if (N == 0 && Z == 1 && L == 0 && S == 1 && abs(mass - kSigmaPMass) < 0.0001)
      PDG = 3222; // sigma+
   else if (N == 0 && Z == -1 && L == 0 && S == 1 && abs(mass - kSigmaPMass) < 0.0001)
      PDG = -3222; // sigma+

   else if (N == 0 && Z == 1 && L == 0 && S == 1 && abs(mass - kSigmaNMass) < 0.0001)
      PDG = 3212; // sigma0
   else if (N == 0 && Z == -1 && L == 0 && S == 1 && abs(mass - kSigmaNMass) < 0.0001)
      PDG = -3212; // sigma0

   else if (N == 0 && Z == 1 && L == 0 && S == 1 && abs(mass - kSigmaMMass) < 0.0001)
      PDG = 3112; // sigma-
   else if (N == 0 && Z == -1 && L == 0 && S == 1 && abs(mass - kSigmaMMass) < 0.0001)
      PDG = 3112; // sigma-
   else
      PDG = 10 * 10E7 + LS * 10E6 + Z * 10E3 + A * 10 + I;

   return PDG;
};
// ---------------------------------------------------------------------

Bool_t MpdPHQMDGenerator::ReadHeader()
{
   /* read header */
   gzgets(fgzFile, fbuffer, 256);     // 1st line
   if (gzeof(fgzFile)) return kFALSE; // end of file reached
   int res = sscanf(fbuffer, "%d %*d %*d %e", &fntr, &fb);
   if (res != 2) {
      printf("-E- MpdPHQMDGenerator: selftest error in header, scan %d of 2\n", res);
      exit(1);
   }
   gzgets(fgzFile, fbuffer, 256); // 2nd line
   if (gzeof(fgzFile)) {
      printf("-E- MpdPHQMDGenerator: unexpected end of file\n");
      exit(1);
   }

   gzgets(fgzFile79, fbuffer79, 256);   // 1st line - unused
   gzgets(fgzFile79, fbuffer79, 256);   // 2nd line - unused
   gzgets(fgzFile79, fbuffer79, 256);   // event header
   if (gzeof(fgzFile79)) return kFALSE; // end of file reached
   int res79 = sscanf(fbuffer79, "%d", &fntr79);
   if (res79 != 1) {
      printf("-E- MpdPHQMDGenerator: selftest error in header, scan %d of 1\n", res79);
      exit(1);
   }
   if (gzeof(fgzFile79)) {
      printf("-E- MpdPHQMDGenerator: unexpected end of baryons file\n");
      exit(1);
   }
   return kTRUE;
};
// ---------------------------------------------------------------------
void MpdPHQMDGenerator::SkipTrack()
{
   gzgets(fgzFile, fbuffer, 256);
   if (gzeof(fgzFile)) {
      printf("-E- MpdPHQMDGenerator: unexpected end of file\n");
      exit(1);
   }
};
// ---------------------------------------------------------------------
Bool_t MpdPHQMDGenerator::SkipEvents(int n)
{
   for (Int_t k = 1; k <= n; ++k) {
      if (!ReadHeader()) return kFALSE; // end of file
      for (Int_t i = 1; i <= fntr; ++i) SkipTrack();
   }
   return kTRUE;
};
// ---------------------------------------------------------------------

ClassImp(MpdPHQMDGenerator);
