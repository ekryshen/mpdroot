/*
 * File:   Mpd3fdGenerator.h
 * Author: Sergei Merts
 *
 * Created on June 22, 2015, 2:06 PM
 */

#ifndef MPD3FDGENERATOR_H
#define MPD3FDGENERATOR_H

#include "FairGenerator.h"
#include "TFile.h"
#include "TChain.h"
#include "TRandom2.h"
#include <fstream>

const UInt_t kBatyukConst = 1e5;

class TVirtualMCStack;
class FairPrimaryGenerator;

class Mpd3fdGenerator : public FairGenerator {
public:
   /** Default constructor without arguments should not be used. **/
   Mpd3fdGenerator();

   /** Standard constructor.
    * @param fileName The input file name
    **/
   Mpd3fdGenerator(TString fileName);

   /** Destructor. **/
   ~Mpd3fdGenerator();

   /** Set fixed reaction plane angle **/
   void SetPsiRP(Double_t PsiRP)
   {
      fPsiRP = PsiRP;
      fisRP  = kFALSE;
   };

   /** Correction factor = Z/A for Theseus 2018-03-17-bc2a06d **/
   void SetProtonNumberCorrection(Float_t c) { fPNCorr = c; };

   Bool_t ReadEvent(FairPrimaryGenerator *primGen);

   void SkipEvents(Int_t ev) { fEventNumber = ev; }

private:
   TFile   *fInputFile;        //!  Input file
   TString  fFileName;         //!  Input file name
   TChain  *fDstTree;          //!
   Double_t fPx[kBatyukConst]; //!
   Double_t fPy[kBatyukConst]; //!
   Double_t fPz[kBatyukConst]; //!
   // Double_t fX[kBatyukConst]; //!  [fm] freezout
   // Double_t fY[kBatyukConst]; //!  [fm] freezout
   // Double_t fZ[kBatyukConst]; //!  [fm] freezout
   Double_t  fE[kBatyukConst];   //!
   Int_t     fPID[kBatyukConst]; //!  pdg
   UInt_t    fNpart;             //!
   Int_t     fEventNumber;       //!
   Float_t   fPNCorr;            //! correction factor for proton multiplicity
   TRandom2 *frandom;            //!
   Double_t  fPsiRP;             //! reaction plane angle
   Bool_t    fisRP;              //! random/fixed reaction plane

   Mpd3fdGenerator(const Mpd3fdGenerator &);
   Mpd3fdGenerator &operator=(const Mpd3fdGenerator &);

   ClassDef(Mpd3fdGenerator, 1);
};

#endif /* MPD3FDGENERATOR_H */
