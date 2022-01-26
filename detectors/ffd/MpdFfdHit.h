//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_FFD_HIT_H
#define __MPD_FFD_HIT_H

#include "FairHit.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdFfdHit : public FairHit {
protected:
   Double_t fTime = 0.; // Time since event start [ns]
   Int_t    fFlag = 0;  // Flag for general purposes [TDC, event tagging...]
   size_t   fNpe  = 0;  // weight: number of pe

public:
   MpdFfdHit();
   MpdFfdHit(Int_t suid, TVector3 pos, TVector3 dpos, Int_t refIndex, Double_t time, size_t npe, Int_t flag = 0);
   virtual ~MpdFfdHit();

   void  Print(const Option_t *opt = 0) const;
   Int_t GetFlag() const { return fFlag; };
   void  SetFlag(Int_t flag) { fFlag = flag; };

   ClassDef(MpdFfdHit, 2)
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_FFD_HIT_H
