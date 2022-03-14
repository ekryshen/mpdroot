/**
 *@class MpdPHQMDGenerator
 *@author V.Kireyeu, V. Voronyuk
 **/

#ifndef Mpd_PHQMDGENERATOR_H
#define Mpd_PHQMDGENERATOR_H

#ifndef __CINT__
#include <zlib.h>
#endif

#include <TRandom2.h>

#include "FairGenerator.h"
#include "FairPrimaryGenerator.h"

class MpdPHQMDGenerator : public FairGenerator {
public:
   MpdPHQMDGenerator();
   /** Default constructor **/
   MpdPHQMDGenerator(const char *filename, const char *f79name);
   /** Destructor **/
   virtual ~MpdPHQMDGenerator();
   /** Read event **/
   Bool_t ReadEvent(FairPrimaryGenerator *primGen);
   /** Skip n events, return kTRUE if successful **/
   Bool_t SkipEvents(Int_t n);

   /** Set fixed reaction plane angle **/
   void SetPsiRP(Double_t PsiRP)
   {
      fPsiRP = PsiRP;
      fisRP  = kFALSE;
   };

   /** Calculate nuclei or baryon PDG **/
   Int_t BaryonPDG(Int_t N, Int_t Z, Int_t L, Int_t S, Float_t mass);

private:
   Int_t   fntr;   //! number of tracks in given event
   Int_t   fntr79; //! number of tracks in given event
   Float_t fb;     //! impact parameter in given event
#ifndef __CINT__
   gzFile fgzFile;   //! file
   gzFile fgzFile79; //! fort.79 output file
#endif
   char fbuffer[256];   //! reading buffer
   char fbuffer79[256]; //! reading buffer for fort.79 file

   TRandom2 *frandom; //!
   Double_t  fPsiRP;  //! reaction plane angle
   Bool_t    fisRP;   //! random/fixed reaction plane

   Bool_t ReadHeader(); //! read event header, return kTRUE if successful
   void   SkipTrack();  //! skip one track

   const Float_t kSigmaPMass = 0.1189E+01;
   const Float_t kSigmaNMass = 0.1192E+01;
   const Float_t kSigmaMMass = 0.1197E+01;
   Float_t       kProtonMass;
   Float_t       kNeutronMass;
   Float_t       kLambdaMass;

   ClassDef(MpdPHQMDGenerator, 0);
};
#endif
