/**
 *@class MpdEPOSGenerator
 *@author K.Shtejer <kshtejer@jinr.ru>
 * MpdEPOSGenerator reads output of EPOS transport model in HepMC format.
 * (crmc_epos199.hepmc / crmc_epos199.hepmc.gz)
 * This interfase assumes that the two header lines and the footer line
 * starting with "HepMC::" have been removed.
 *
 * Last updates: June 8, 2018
 *
 **/

#include <stdio.h>
#include <stdlib.h>

#include <TMath.h>

#include "FairMCEventHeader.h"
#include "MpdMCEventHeader.h"
#include "MpdEPOSGenerator.h"

// ---------------------------------------------------------------------
MpdEPOSGenerator::MpdEPOSGenerator(){
   // It is better to leave empty
};
// ---------------------------------------------------------------------
MpdEPOSGenerator::MpdEPOSGenerator(const char *filename = "crmc_epos199.hepmc.gz")
{
   fgzFile = gzopen(filename, "rb"); // zlib
   if (!fgzFile) {
      printf("-E- MpdEPOSGenerator: can not open file: %s\n", filename);
      exit(1);
   }
   printf("-I- MpdEPOSGenerator: open %s\n", filename);

   fPsiRP  = 0.;
   fisRP   = kTRUE; // by default RP is random
   frandom = new TRandom2();
   frandom->SetSeed(0);
};
// ---------------------------------------------------------------------
MpdEPOSGenerator::~MpdEPOSGenerator()
{
   if (fgzFile) {
      gzclose(fgzFile);
      fgzFile = NULL;
   }
   delete frandom;
};
// ---------------------------------------------------------------------
Bool_t MpdEPOSGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{

   gzgets(fgzFile, fbuffer, 256);     // 1st line "E": event
   if (gzeof(fgzFile)) return kFALSE; // end of file reached
   res = sscanf(fbuffer, "%s %*d %*d %*e %*e %*e %*d %*d %d %*[^\n]", charId, &numVtx);
   if (strcmp(charId, "E") == 0) {
      printf(" => %s \n", charId);
   }
   if (res != 2) {
      printf("-E- MpdEPOSGenerator: selftest error in header, scan %d of 2\n", res);
      exit(1);
   }

   gzgets(fgzFile, fbuffer, 256); // 2st line "U": unities
   if (gzeof(fgzFile)) return kFALSE;
   res = sscanf(fbuffer, "%s %*[^\n]", charId);
   if (strcmp(charId, "U") == 0) {
      printf(" => %s \n", charId);
   }
   if (res != 1) {
      printf("-E- MpdEPOSGenerator: selftest error in header, scan %d of 1\n", res);
      exit(1);
   }

   gzgets(fgzFile, fbuffer, 256); // 3st line "C": cross section
   if (gzeof(fgzFile)) return kFALSE;
   res = sscanf(fbuffer, "%s %*[^\n]", charId);
   if (strcmp(charId, "C") == 0) {
      printf(" => %s \n", charId);
   }
   if (res != 1) {
      printf("-E- MpdEPOSGenerator: selftest error in header, scan %d of 1\n", res);
      exit(1);
   }

   gzgets(fgzFile, fbuffer, 256); // 4st line "H": heavy ion info
   if (gzeof(fgzFile)) return kFALSE;
   res = sscanf(fbuffer, "%s %*d %*d %*d %*d %*d %*d %*d %*d %*d %e %*[^\n]", charId, &fb);
   if (strcmp(charId, "H") == 0) {
      printf(" => %s \n", charId);
   }
   if (res != 2) {
      printf("-E- MpdEPOSGenerator: selftest error in header, scan %d of 2\n", res);
      exit(1);
   }

   gzgets(fgzFile, fbuffer, 256); // 5st line "F"
   if (gzeof(fgzFile)) return kFALSE;
   res = sscanf(fbuffer, "%s %*d %*[^\n]", charId);
   if (strcmp(charId, "F") == 0) {
      printf(" => %s \n", charId);
   }
   if (res != 1) {
      printf("-E- MpdEPOSGenerator: selftest error in header, scan %d of 1\n", res);
      exit(1);
   }

   printf("-I- Info from EPOS-> Event-> numVtx, ImpParam(fm) : %d, %e \n", numVtx, fb);

   /* set random Reaction Plane angle */
   if (fisRP) {
      fPsiRP = frandom->Uniform(2.0 * TMath::Pi());
   }
   // printf("-I- MpdPHSDGenerator: RP angle = %e\n",fPsiRP);

   /* Set event impact parameter in MCEventHeader if not yet done */
   FairMCEventHeader *eventHeader = primGen->GetEvent();
   if (eventHeader && (!eventHeader->IsSet())) {
      eventHeader->SetB(fb);
      eventHeader->MarkSet(kTRUE);
      eventHeader->SetRotZ(fPsiRP);
   }

   // read vertexes //
   for (Int_t j = 1; j <= numVtx; j++) {
      gzgets(fgzFile, fbuffer, 256); // line "V": Vertex
      res  = sscanf(fbuffer, "%s %*d %d %*e %*e %*e %*e %d %d %*[^\n]", charId, &ivtx, &inpart, &outpart);
      fntr = inpart + outpart;

      // read particles //
      for (Int_t i = 1; i <= fntr; i++) {
         gzgets(fgzFile, fbuffer, 256); // line "P": Particle
         res = sscanf(fbuffer, "%s %d %d %e %e %e %*e %e %d %*[^\n]", charId, &itrk, &ipdg, &px, &py, &pz, &imass,
                      &istatus);
         if (istatus == 1) {
            // add track to simulation //
            primGen->AddTrack(ipdg, px, py, pz, 0., 0., 0.);
            printf("\n -I- Info from EPOS-> Track -> itrk, istatus, ipdg, px, py, pz : %d %d %d %e %e %e", itrk,
                   istatus, ipdg, px, py, pz);
         } else {
            printf("\n -I- Info from EPOS-> Skipped Track -> itrk, istatus : %d %d", itrk, istatus);
         }
      }

      if (gzeof(fgzFile)) return kFALSE;
   }

   return kTRUE;
};

// ---------------------------------------------------------------------
ClassImp(MpdEPOSGenerator);
