//------------------------------------------------------------------------------------------------------------------------

#include <iostream>

#include "MpdTofPoint.h"

using namespace std;

ClassImp(MpdTofPoint);
//------------------------------------------------------------------------------------------------------------------------
MpdTofPoint::MpdTofPoint() : FairMCPoint() {}
//------------------------------------------------------------------------------------------------------------------------
MpdTofPoint::MpdTofPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t tof, Double_t length,
                         Double_t eLoss)
   : FairMCPoint(trackID, detID, pos, mom, tof, length, eLoss)
{
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofPoint::Print(const Option_t *opt) const
{
   cout << "\n-I- MpdTofPoint: tid=" << fTrackID << ", duid=" << fDetectorID << ", Position(" << fX << ", " << fY
        << ", " << fZ << ") cm"
        << ", Momentum(" << fPx << ", " << fPy << ", " << fPz << ") GeV/c"
        << ", time=" << fTime << " ns, tr. length " << fLength << " cm,  Eloss=" << fELoss * 1.e06 << " keV.";
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofPoint::PrintSuid(Int_t suid, const char *comment, ostream &os)
{
   if (comment != nullptr) os << comment;

   Int_t sector, detector, gap, strip;
   ParseSuid(suid, sector, detector, gap, strip);

   os << " suid=0x" << hex << suid << dec << " : sector=" << sector << ", detector=" << detector << ", gap=" << gap
      << ", strip=" << strip;
   return true;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofPoint::ParseSuid(Int_t suid, Int_t &sector, Int_t &detector, Int_t &gap, Int_t &strip)
{
   sector   = GetSector(suid);
   detector = GetDetector(suid);
   gap      = GetGap(suid);
   strip    = GetStrip(suid);
}
//------------------------------------------------------------------------------------------------------------------------
Int_t MpdTofPoint::GetSuid72(Int_t sector, Int_t detector, Int_t strip)
{
   const Int_t nStrips = 24;

   Int_t gap = (strip - 1) / nStrips + 1; // gap range [1,3]
   strip     = (strip - 1) % nStrips + 1; // strip range [1,72]  -> [1,24]

   return GetSuid24(sector, detector, gap, strip);
}
//------------------------------------------------------------------------------------------------------------------------
Int_t MpdTofPoint::GetSuid24(Int_t sector, Int_t detector, Int_t gap, Int_t strip)
{
#ifdef DEBUG
   Int_t suid       = (sector << 24) + (detector << 16) + (gap << 8) + strip;
   Int_t sectorID   = GetSector(suid);
   Int_t detectorID = GetDetector(suid);
   Int_t gapID      = GetGap(suid);
   Int_t stripID    = GetStrip(suid);
   assert(sector == sectorID);
   assert(detector == detID);
   assert(gap == gapID);
   assert(strip == stripID);
   return suid;
#else
   return (sector << 24) + (detector << 16) + (gap << 8) + strip;
#endif
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofPoint::IsSameDetector(Int_t suid1, Int_t suid2)
{
#ifdef DEBUG
   assert(((suid1 & 0xFFFF0000) == (suid2 & 0xFFFF0000)) ==
          (GetSector(suid1) == GetSector(suid2) && GetDetector(suid1) == GetDetector(suid2)));
#endif
   return ((suid1 & 0xFFFF0000) == (suid2 & 0xFFFF0000));
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofPoint::IsSameStrip(Int_t suid1, Int_t suid2)
{
#ifdef DEBUG
   assert(((suid1 & 0xFFFF00FF) == (suid2 & 0xFFFF00FF)) ==
          (GetSector(suid1) == GetSector(suid2) && GetDetector(suid1) == GetDetector(suid2) &&
           GetStrip(suid1) == GetStrip(suid2)));
#endif
   return ((suid1 & 0xFFFF00FF) == (suid2 & 0xFFFF00FF));
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofPoint::IsSameGap(Int_t suid1, Int_t suid2)
{
#ifdef DEBUG
   assert(((suid1 & 0x0000FF00) == (suid2 & 0x0000FF00)) == (GetGap(suid1) == GetGap(suid2)));
#endif
   return ((suid1 & 0x0000FF00) == (suid2 & 0x0000FF00));
}
//------------------------------------------------------------------------------------------------------------------------
