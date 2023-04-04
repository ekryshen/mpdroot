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
//      Sergey Merts
//
//-----------------------------------------------------------

#ifndef MPDPARTICLEIDENTIFICATION_HH
#define MPDPARTICLEIDENTIFICATION_HH
// Base Class Headers ----------------
#include "TSystem.h"
#include "MpdTPCpid.h"
#include "MpdTOFpid.h"

class MpdParticleIdentification : public MpdTPCpid, public MpdTOFpid {
public:
   // Constructors/Destructors ---------
   MpdParticleIdentification();
   virtual ~MpdParticleIdentification();

   // Int_t GetTpcProbs(Float_t P, Float_t dedx, Int_t nHits, Float_t& Ppi, Float_t& PK, Float_t& Pp, Float_t& Pe, Int_t
   // method); The function is realized now in the MpdTPCpid-class with corrected probability coefficients Int_t
   // GetTofProbs(Float_t P, Float_t beta, Float_t& Ppi, Float_t& PK, Float_t& Pp, Float_t& Pe, Int_t method);
   // The function is realized now in the detectors/tof/MpdTOFpid-class with corrected probability coefficients
   Int_t GetCombinedProbs(Float_t *tofProbs, Float_t *tpcProbs, Float_t *resultProbs, Int_t N);

   // Float_t HyperbolicFunction(Float_t x, Float_t *p);
   // Float_t BetheBlochFunction(Float_t x, Float_t *p);
   Float_t ParabolicFunction(Float_t x, Float_t *p);
   // Int_t BayesFunction(Float_t *measProb, Float_t *aprioriProb, Float_t *bayesProb, Int_t N);

private:
public:
   ClassDef(MpdParticleIdentification, 1);
};
#endif // #ifndef MPDPARTICLEIDENTIFICATION_HH
