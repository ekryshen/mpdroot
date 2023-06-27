//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_HIT_H
#define __MPD_TOF_HIT_H

#include <vector>

#include "MpdTofUtils.h"
#include "FairHit.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdTofHit : public FairHit {
protected:
   Double_t fTime; // Time since event start [ns]
   Int_t    fFlag; // Flag for general purposes [TDC, event tagging...]
   Double_t fAmplitude = 0.;

public:
   MpdTofHit();
   MpdTofHit(Int_t suid, TVector3 pos, TVector3 dpos, Int_t refIndex, Double_t tof, Int_t flag = 0, Double_t amplitude = 0.);

   void Print(const char *comment = nullptr, std::ostream &os = std::cout) const;

   Double_t GetTime() const { return fTime; };
   Int_t    GetFlag() const { return fFlag; };
   Double_t GetAmplitude() const { return fAmplitude; };

   void SetTime(Double_t time) { fTime = time; };
   void SetFlag(Int_t flag) { fFlag = flag; };
   void AddFlag(Int_t flag) { fFlag = fFlag | flag; };

   bool CheckSuid(Int_t suid) const;
   bool CheckTid(Int_t tid) const;
   bool IsSameTid(const MpdTofHit &) const;
   void getLinks(const MpdTofUtils::k_LinkType type, std::vector<Int_t> &) const;

   inline void GetXYZ(Double_t *carray) const
   {
      carray[0] = fX;
      carray[1] = fY;
      carray[2] = fZ;
   }

   ClassDef(MpdTofHit, 2);
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_HIT_H
