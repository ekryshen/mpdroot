#ifndef TPCLASERGRIDPREPROCESS_H
#define TPCLASERGRIDPREPROCESS_H

#include "FairTask.h"
#include "TClonesArray.h"
#include "TVector3.h"

class TpcLaserGridPreprocess : public FairTask {
public:
   TpcLaserGridPreprocess();
   ~TpcLaserGridPreprocess();

   void ClearSecondaries(Bool_t opt = kTRUE) { clrSecondaries = opt; }
   void SetMakeQA(Bool_t opt = kFALSE) { makeQA = opt; }
   void SetBeamsCount(Int_t bc) { beamsCount = bc; }
   void SetMCPointsSmearing(Bool_t opt = kTRUE) { mcPointsSmearing = opt; }

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   virtual void       Finish();

private:
   TString       inputBranchName;
   TClonesArray *mcPointArray;
   TClonesArray *mcTrackArray;

   Int_t          beamsCount          = 224;
   const Double_t centerBeamMaxLength = 130.;
   const Double_t tubeInR             = 40.;
   const Double_t beamRadius          = 0.05;

   Bool_t clrSecondaries   = kTRUE;
   Bool_t mcPointsSmearing = kTRUE;
   Bool_t makeQA           = kFALSE;

   Double_t distance_to_line(TVector3 pt, TVector3 l1, TVector3 l2) const;

   ClassDef(TpcLaserGridPreprocess, 0)
};

#endif
