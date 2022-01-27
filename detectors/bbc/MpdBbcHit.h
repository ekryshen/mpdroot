#ifndef MPDBBCHIT_H
#define MPDBBCHIT_H

#include "FairHit.h"

class MpdBbcHit : public FairHit {

public:
   /** Default constructor **/
   MpdBbcHit();

   /** Constructor with hit parameters (1)**/
   MpdBbcHit(Int_t detectorID, TVector3 pos, TVector3 dpos, Int_t refIndex, Int_t flag);

   /** Constructor with hit parameters (2) [not the flag]**/
   MpdBbcHit(Int_t detectorID, TVector3 pos, TVector3 dpos, Int_t refIndex);

   /** Destructor **/
   virtual ~MpdBbcHit();

   void Print(const Option_t *opt = 0) const;

   /** Accessors **/
   Int_t GetFlag() const { return fFlag; };

   /** Modifiers **/
   void SetFlag(Int_t flag) { fFlag = flag; };

protected:
   Int_t fTrackID;
   Int_t fFlag; // Flag for general purposes [TDC, event tagging...]

   ClassDef(MpdBbcHit, 1);
};

#endif
