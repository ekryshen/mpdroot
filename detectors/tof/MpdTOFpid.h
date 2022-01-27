//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      class for particle identification
//
//
// Environment:
//      Software developed for MPD at NICA.
//
// Author List:
//      Gyulnara Eyyubova
//
//-----------------------------------------------------------

#ifndef MpdTOFpid_HH
#define MpdTOFpid_HH

// Base Class Headers ----------------
#include "TSystem.h"
#include <TGraph.h>
#include <TFile.h>

class MpdTOFpid {
public:
   // Constructors/Destructors ---------
   MpdTOFpid();
   virtual ~MpdTOFpid();

   Int_t   GetTofProbs(Float_t P, Float_t beta, Float_t &Ppi, Float_t &Pk, Float_t &Pp, Float_t &Pe, Int_t method);
   Float_t HyperbolicFunction(Float_t x, Float_t *p);

   Int_t BayesFunction(Float_t *measProb, Float_t *aprioriProb, Float_t *bayesProb, Int_t N);
   void  ReadTOFResponse();
   void  SetNSigmaBeta(Float_t n) { fnBeta = n; }

private:
   Float_t fnBeta; // number of sigmas for n-sigma method
   Float_t ProtonParHyper[2];
   Float_t PionParHyper[2];
   Float_t KaonParHyper[2];
   Float_t ElectronParHyper[2];

   Float_t sigmasTofPi[72];
   Float_t sigmasTofPr[72];
   Float_t sigmasTofKa[72];
   Float_t sigmasTofEl[72];

public:
   ClassDef(MpdTOFpid, 1);
};
#endif
