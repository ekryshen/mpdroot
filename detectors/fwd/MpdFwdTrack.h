//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdTrack
/// \brief Reconstructed track 
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
   void SetStateVector(Double_t* v) { for (Int_t i=0;i<6;i++) fV[i] = v[i]; }
   Double_t GetZ()   { return fV[0]; }
   Double_t GetX()   { return fV[1]; }
   Double_t GetY()   { return fV[2]; }
   Double_t GetTy()  { return fV[3]; }
   Double_t GetTz()  { return fV[4]; }
   Double_t GetQpt() { return fV[5]; }
private:
   std::vector<Int_t> fHitIndices{};
   Double_t fV[6]; // state vector: z, x, y, ty = py/px, tz = pt/pz, q/pt
   ClassDef(MpdFwdTrack, 5)
};
#endif // #ifndef __MPD_FWD_HIT_H
