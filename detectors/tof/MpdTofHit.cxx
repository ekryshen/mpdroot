#include <iostream>
#include <cmath>

#include "MpdTofUtils.h"
#include "MpdTofHit.h"

using namespace std;
//------------------------------------------------------------------------------------------------------------------------
MpdTofHit::MpdTofHit()
{
   fTime = fX = fY = fZ = fDx = fDy = fDz = NAN;
   fFlag                                  = 0;
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofHit::MpdTofHit(Int_t suid, TVector3 pos, TVector3 dpos, Int_t index, Double_t time, Int_t flag, Double_t ampl)
   : FairHit(suid, pos, dpos, index), fTime(time), fFlag(flag), fAmplitude(ampl)
{
   
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofHit::CheckSuid(Int_t suid) const
{
   for (Int_t i = 0, nLinks = GetNLinks(); i < nLinks; i++) {
      auto link = GetLink(i);
      if (MpdTofUtils::volumeUID == link.GetType() && suid == link.GetIndex()) return true;
   }

   return false;
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofHit::CheckTid(Int_t tid) const
{
   for (Int_t i = 0, nLinks = GetNLinks(); i < nLinks; i++) {
      auto link = GetLink(i);
      if (MpdTofUtils::mcTrackIndex == link.GetType() && tid == link.GetIndex()) return true;
   }

   return false;
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofHit::IsSameTid(const MpdTofHit &hit) const
{
   for (Int_t i = 0, nLinks = hit.GetNLinks(); i < nLinks; i++) {
      auto link = hit.GetLink(i);
      if (MpdTofUtils::mcTrackIndex == link.GetType() && CheckTid(link.GetIndex())) return true;
   }

   return false;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHit::getLinks(const MpdTofUtils::k_LinkType type, vector<Int_t> &v) const
{
   v.clear();

   for (Int_t i = 0, nLinks = GetNLinks(); i < nLinks; i++) {
      auto link = GetLink(i);
      if (type == link.GetType()) v.push_back(link.GetIndex());
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHit::Print(const char *comment, ostream &os) const
{
   if (comment != nullptr) os << comment;

   os << " [MpdTofHit] uid: " << fDetectorID << " pos.=("<<fX<<", "<<fY<<", "<<fZ
	<<") cm,  pos. err.=(" << fDx << ", " << fDy << ", " << fDz 
	<<") cm, time=" << fTime << " ns, flag="<<fFlag<<", ampl.="<<fAmplitude;

   if(GetNLinks() == 0) return; // for hits from MpdTofDigit( MC data is absent)

   os << " MCtid: ";
   for (Int_t i = 0, nLinks = GetNLinks(); i < nLinks; i++) {
      auto link = GetLink(i);
      if (link.GetType() == MpdTofUtils::mcTrackIndex) os << link.GetIndex() << ", ";
   }

   os << " MCpointId: ";
   for (Int_t i = 0, nLinks = GetNLinks(); i < nLinks; i++) {
      auto link = GetLink(i);
      if (link.GetType() == MpdTofUtils::mcPointIndex) os << link.GetIndex() << ", ";
   }

   os << " MCsuid: ";
   for (Int_t i = 0, nLinks = GetNLinks(); i < nLinks; i++) {
      auto link = GetLink(i);
      if (link.GetType() == MpdTofUtils::volumeUID) os << link.GetIndex() << ", ";
   }
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdTofHit);
