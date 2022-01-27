//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_HIT_PRODUCER_IDEAL_H
#define __MPD_TOF_HIT_PRODUCER_IDEAL_H

#include <TVector3.h>
#include <TList.h>

#include "FairTask.h"
//------------------------------------------------------------------------------------------------------------------------
class TClonesArray;
class MpdTofHitProducerQA;
//------------------------------------------------------------------------------------------------------------------------
class MpdTofHitProducerIdeal : public FairTask {
protected:
   TClonesArray *aMcPoints  = nullptr; //! <--- MC input
   TClonesArray *aMcTracks  = nullptr; //! <--- MC input
   TClonesArray *aExpDigits = nullptr; //! <--- Exp input
   TClonesArray *aTofHits   = nullptr; //! ---> output

   Bool_t fDoTest;
   Bool_t fDoMergeHits;
   Bool_t fUseMCData;
   Bool_t fOnlyPrimary = false;

   MpdTofHitProducerQA *pQA = nullptr; //!

   void  AddHit(Int_t detUID, const TVector3 &posHit, const TVector3 &posHitErr, Int_t mcPointIndex, Int_t mcTrackIndex,
                Double_t time, Int_t flag);
   void  AddHit(Int_t detUID, const TVector3 &posHit, const TVector3 &posHitErr, Int_t expDigitIndex, Double_t time,
                Int_t flag);
   Int_t CompressHits();
   size_t MergeHitsOnStrip(); // save only the fastest hit in the strip

public:
   MpdTofHitProducerIdeal(const char *name = "TOF Ideal HitProducer", Bool_t useMCdata = true, Int_t verbose = 1,
                          Bool_t DoTest = false, Bool_t DoMergeHits = false,
                          const char *flnm = "QA.MpdTofHitProducerIdeal.root", bool IsEndcap = false);
   virtual ~MpdTofHitProducerIdeal();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *option);
   virtual void       Finish();

   void SetOnlyPrimary(Bool_t opt = true) { fOnlyPrimary = opt; }
   void Dump(const char *comment = nullptr, std::ostream &out = std::cout) const;

   virtual void AddParameters(TString &buf) const;

   ClassDef(MpdTofHitProducerIdeal, 3);
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_HIT_PRODUCER_IDEAL_H
