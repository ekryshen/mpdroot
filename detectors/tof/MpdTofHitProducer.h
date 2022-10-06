//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_HIT_PRODUCER_H
#define __MPD_TOF_HIT_PRODUCER_H

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofHitProducer
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <iostream>

#include <TString.h>

#include "MpdTofHitProducerIdeal.h"
//------------------------------------------------------------------------------------------------------------------------
class TRandom2;
class TEfficiency;
class TH1D;
class TH2D;

class FairMCEventHeader;
//------------------------------------------------------------------------------------------------------------------------
class MpdTofHitProducer : public MpdTofHitProducerIdeal {
   Double_t  fTimeSigma = 0.080;                    // Uncertainties of time, gaus sigma [ns]
   Double_t  fStripLength;                          // [cm] strip length, default: 65 cm
   Double_t  fErrX = 0.5, fErrZ = 1.25 / sqrt(12.); // Uncertainties of coordinates, gaus sigma [cm]
   TRandom2 *pRandom = nullptr;

   // value - distance to the strip edge [cm], gap - 1,3 = outer strip gap, 2 = middle strip gap
   Bool_t IsHitCreated(Double_t value, Int_t gap);
   Bool_t IsCrossHitCreated(Double_t value, Int_t gap);
//V
   FairMCEventHeader *fEventHeaderB  = nullptr;
   Float_t b_mc, t0res_mc;

public:
   MpdTofHitProducer(const char *name = "TOF HitProducer", Bool_t useMCdata = true, Int_t verbose = 1,
                     Bool_t DoTest = false, const char *flnm = "QA.MpdTofHitProducer.root", double StripLength = 65.);
   virtual ~MpdTofHitProducer();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *option);
   virtual void       Finish();

   void SetTimeResolution(Double_t sigma) { fTimeSigma = sigma; };
   void SetAlongStripzResolution(Double_t Xerr) { fErrX = Xerr; };
   void SetSeed(UInt_t seed = 0);

   virtual void AddParameters(TString &buf) const;

   ClassDef(MpdTofHitProducer, 3);
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_HIT_PRODUCER_H
