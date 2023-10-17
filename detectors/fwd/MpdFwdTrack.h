//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdPoint
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------

#ifndef __MPD_FWD_TRACK_H
#define __MPD_FWD_TRACK_H

#include "TObject.h"
#include "vector"

class MpdFwdTrack : public TObject {
public:
   MpdFwdTrack();
   virtual ~MpdFwdTrack();
   void AddHitIndex(Int_t hitIndex);
   void Print();
   Int_t GetNIndices() { return fHitIndices.size(); }
   Int_t GetHitIndex(Int_t i) { return fHitIndices[i]; }
   Double_t GetPtMC() {return fPtMC;}
   void SetPtMC(Double_t ptMC) { fPtMC = ptMC; }
   Double_t GetPMC() {return fPMC;}
   void SetPMC(Double_t pMC) { fPMC = pMC; }
   Double_t GetPdgCode() {return fPdgCode;}
   void SetPdgCode(Double_t pdgcode) { fPdgCode = pdgcode; }
private:
   std::vector<Int_t> fHitIndices{};
   Double_t fPtMC;
   Double_t fPMC;
   Int_t fPdgCode;
   ClassDef(MpdFwdTrack, 5)
};
#endif // #ifndef __MPD_FWD_HIT_H
