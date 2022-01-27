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

#ifndef MpdTPCpid_HH
#define MpdTPCpid_HH

// Base Class Headers ----------------
#include "TSystem.h"
#include <TGraph.h>
#include <TFile.h>

class MpdTPCpid {
public:
   // Constructors/Destructors ---------
   MpdTPCpid();
   virtual ~MpdTPCpid();

   Int_t   GetTpcProbs(Float_t P, Float_t dedx, Int_t nHits, Float_t &Ppi, Float_t &PK, Float_t &Pp, Float_t &Pe,
                       Int_t method);
   Float_t BetheBlochFunction(Float_t x, Float_t *p);
   Int_t   BayesFunction(Float_t *measProb, Float_t *aprioriProb, Float_t *bayesProb, Int_t N);
   void    ReadTPCResponse();
   void    SetNSigmaDedx(Float_t n) { fn = n; }

private:
   Float_t fn; // number of sigmas for n-sigma method

   Float_t ProtonPar[6];
   Float_t PionPar[6];
   Float_t KaonPar[6];
   Float_t ElectronPar[6];

   Float_t sigmasPi[72];
   Float_t sigmasPr[72];
   Float_t sigmasKa[72];
   Float_t sigmasEl[72];

public:
   ClassDef(MpdTPCpid, 1)
};

#endif
