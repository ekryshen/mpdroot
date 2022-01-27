//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofMatching
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <bitset>
#include <cmath>

#include <TMath.h>

#include "MpdTofHit.h"
#include "MpdTofMatching.h"

#include "MpdTofMatchingData.h"

using namespace std;

ClassImp(MpdTofMatchingData);

//------------------------------------------------------------------------------------------------------------------------
MpdTofMatchingData::MpdTofMatchingData()
{
   fTime = fBeta = fMass2 = fLength = NAN;
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofMatchingData::MpdTofMatchingData(Int_t kfTrackId, Int_t tofHitId, Double_t weight, const MpdTofHit *pHit,
                                       Double_t length, Int_t nTrHits, const TVector3 &P, const TVector3 &point)
   : fKFTrackIndex(kfTrackId), fTofHitIndex(tofHitId), fWeight(weight), fLength(length), fNmbTrHits(nTrHits)
{
   fFlag     = pHit->GetFlag();
   fTime     = pHit->GetTime();
   fMomentum = P;
   pHit->Position(fHitPosition);
   fEstPoint = point;

   fdPhi = TMath::Sqrt(fEstPoint.X() * fEstPoint.X() + fEstPoint.Y() * fEstPoint.Y()) *
           (TMath::ATan2(fEstPoint.Y(), fEstPoint.X()) - TMath::ATan2(fHitPosition.Y(), fHitPosition.X()));
   // fdPhi = TMath::Sqrt( (fEstPoint.X()-fHitPosition.X())*(fEstPoint.X()-fHitPosition.X()) +
   // (fEstPoint.Y()-fHitPosition.Y())*(fEstPoint.Y()-fHitPosition.Y()) );
   fdZed = fEstPoint.Z() - fHitPosition.Z();

   const static Double_t c_vel = TMath::C();
   fBeta                       = fLength / (fTime * 1.e-7) / c_vel; // [cm/nc] -> m/c
   Double_t beta2              = fBeta * fBeta;
   Double_t gamma2             = 1. / (1. - beta2);
   fMass2                      = fMomentum.Mag2() / (gamma2 * beta2);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingData::Print(const char *comment, ostream &os) const
{
   if (comment != nullptr) os << comment;

   os << "[MpdTofMatchingData] pair<" << fKFTrackIndex << ", " << fTofHitIndex << ">, weight=" << fWeight
      << ", norm. weight=" << fNormWeight << ", dev.=" << GetDelta() << ", flag=" << std::bitset<sizeof(Int_t)>(fFlag)
      << ", Momentum(" << fMomentum.X() << "," << fMomentum.Y() << "," << fMomentum.Z() << "), hit position("
      << fHitPosition.X() << "," << fHitPosition.Y() << "," << fHitPosition.Z() << "), length=" << fLength;
}
//------------------------------------------------------------------------------------------------------------------------
