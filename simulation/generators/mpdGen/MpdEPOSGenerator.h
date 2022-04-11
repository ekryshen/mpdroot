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

#ifndef Mpd_EPOSGENERATOR_H
#define Mpd_EPOSGENERATOR_H

#ifndef __CINT__
#include <zlib.h>
#endif

#include <TRandom2.h>

#include "FairGenerator.h"
#include "FairPrimaryGenerator.h"

class MpdEPOSGenerator : public FairGenerator {
public:
   MpdEPOSGenerator();
   /** Default constructor **/
   MpdEPOSGenerator(const char *filename);
   /** Destructor **/
   virtual ~MpdEPOSGenerator();
   /** Read event **/
   Bool_t ReadEvent(FairPrimaryGenerator *primGen);
   /** Skip n events, return kTRUE if successful **/
   // Bool_t SkipEvents(Int_t n);

   /** Set fixed reaction plane angle **/
   void SetPsiRP(Double_t PsiRP)
   {
      fPsiRP = PsiRP;
      fisRP  = kFALSE;
   };

private:
   Int_t   res;        //! number of scans in lines
   char    charId[2];  //! First character identifying each line
   Float_t fb;         //! impact parameter in given event
   Int_t   numVtx;     //! number of vertexes in given event
   Int_t   ivtx;       //! id of vertex (Vtx defined by incoming & outgoing particles)
   Int_t   inpart;     //! number of incoming particles "orphan"
   Int_t   outpart;    //! number of outgoing particles
   Int_t   fntr;       //! number of tracks in given vertex
   Int_t   itrk;       //! id of track (particle)
   Int_t   ipdg;       //! PDG code
   Float_t px, py, pz; //! Momemntum components
   Float_t imass;      //! Generated mass for this particle
   Int_t   istatus;    //! Code status (decayed or not)

#ifndef __CINT__
   gzFile fgzFile; //! file
#endif
   char fbuffer[256]; //! reading buffer

   TRandom2 *frandom; //!
   Double_t  fPsiRP;  //! reaction plane angle
   Bool_t    fisRP;   //! random/fixed reaction plane

   ClassDef(MpdEPOSGenerator, 0);
};
#endif
