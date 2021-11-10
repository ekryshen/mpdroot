/// \class MpdTpcHit
/// 
/// Hit in MPD TPC
/// \author Alexander Zinchenko (LHEP, JINR, Dubna) - extension of TpcHit

//#include "TpcCommon.h"
#include "MpdTpcHit.h"
#include "MpdKalmanHit.h"
#include "MpdTpc2dCluster.h"
#include "FairLinkManager.h"

//---------------------------------------------------------------------------
MpdTpcHit::MpdTpcHit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index)
  : FairHit(detID, pos, dpos, index),
    fiPad(-1), fiBin(-1), fLayer(-1), fNdigits(0), fFlag(0), fQ(0), fStep(0), fLength(0), fLocalX(0), fLocalY(0), fLocalZ(0)
{

}

//---------------------------------------------------------------------------
Int_t MpdTpcHit::Compare(const TObject* hit) const
{
  /// "Compare" function to sort in descending order in radius

  MpdTpcHit *tpcHit = (MpdTpcHit*) hit;
  if (fLayer < tpcHit->GetLayer()) return 1;
  else if (fLayer > tpcHit->GetLayer()) return -1;
  else {
    if (GetR() < tpcHit->GetR()) return 1;
    else if (GetR() > tpcHit->GetR()) return -1;
    else return 0;
  }
}
//---------------------------------------------------------------------------

Int_t MpdTpcHit::GetTrackID() const
{
  // Returns track ID (from the link with the smallest weight)

  if (FairLinkManager::Instance() == NULL) return fIDs[0];

  FairMultiLinkedData links = GetLinksWithType(MpdTpcHit::MCTrackIndex);
  Int_t nLinks = links.GetNLinks();
  Int_t id = links.GetLink(0).GetIndex();
  Float_t w = links.GetLink(0).GetWeight();

  for (Int_t i = 1; i < nLinks; ++i) {
    FairLink link = links.GetLink(i);
    if (link.GetWeight() < w) {
      id = link.GetIndex();
      w = link.GetWeight();
    }
  }
  return id;
}

//---------------------------------------------------------------------------

void MpdTpcHit::SetFlags(const MpdTpc2dCluster *clus)
{
  // Transfer flags from the cluster 
  
  if (clus->Virtual()) fFlag |= MpdKalmanHit::kVirtual;
  if (clus->Overflows()) fFlag |= MpdKalmanHit::kOverflow;
  if (clus->Edge()) fFlag |= MpdKalmanHit::kEdge;
  if (clus->Flag() & 1) fFlag |= MpdKalmanHit::kMlem;
  if (clus->MultMax()) fFlag |= MpdKalmanHit::kMultMax;
  if (clus->SinglePad()) fFlag |= MpdKalmanHit::kSinglePad;
}

//---------------------------------------------------------------------------
ClassImp(MpdTpcHit)
