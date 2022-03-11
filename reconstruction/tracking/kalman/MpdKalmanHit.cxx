/// \class MpdKalmanHit
///
/// Hit object for Kalman tracking in the MPD detector
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdKalmanHit.h"
#include <TMath.h>

//__________________________________________________________________________
MpdKalmanHit::MpdKalmanHit()
   : TObject(), fDetectorID(0), fFlag(0), fLength(0.), fIndex(0), fHitType(kFixedP), fNofDim(2), fSignal(0.), fDist(0.),
     fEdge(99.)
{
   /// Default constructor

   for (Int_t i = 0; i < 2; ++i) fMeas[i] = fErr[i] = fCosSin[i] = 0;
   fCosSin[0] = 1.;
}

//__________________________________________________________________________
MpdKalmanHit::MpdKalmanHit(Int_t detID, Int_t nDim, HitType hitType, Double_t *meas, Double_t *err, Double_t *cosSin,
                           Double_t signal, Double_t dist, Int_t index, Double_t edge)
   : TObject(), fDetectorID(detID), fFlag(0), fLength(0.),
     // fIndex(index),
     fHitType(hitType), fNofDim(nDim), fSignal(signal), fDist(dist), fEdge(edge)
{
   /// Constructor

   for (Int_t i = 0; i < 2; ++i) {
      fMeas[i]   = meas[i];
      fErr[i]    = err[i];
      fCosSin[i] = cosSin[i];
   }
   SetIndex(index);
}

//__________________________________________________________________________
MpdKalmanHit::MpdKalmanHit(const MpdKalmanHit &hit)
   : TObject(hit), fDetectorID(hit.fDetectorID), fFlag(hit.fFlag), fLength(hit.fLength), fIndex(hit.fIndex),
     fHitType(hit.fHitType), fNofDim(hit.fNofDim), fSignal(hit.fSignal), fDist(hit.fDist), fEdge(hit.fEdge)
{
   /// copy constructor

   for (Int_t i = 0; i < 2; ++i) {
      fMeas[i]   = hit.fMeas[i];
      fErr[i]    = hit.fErr[i];
      fCosSin[i] = hit.fCosSin[i];
   }
}

//__________________________________________________________________________
MpdKalmanHit::~MpdKalmanHit()
{
   /// Destructor
}

//__________________________________________________________________________
Int_t MpdKalmanHit::Compare(const TObject *hit) const
{
   /// "Compare" function to sort in descending order in fDist

   MpdKalmanHit *kHit = (MpdKalmanHit *)hit;
   if (kHit->GetType() == GetType()) {
      // Check layers
      if (GetLayer() < kHit->GetLayer())
         return 1;
      else if (GetLayer() > kHit->GetLayer())
         return -1;
      if (GetType() == kFixedP) {
         // Sort according to sector number
         if (GetDetectorID() % 1000000 < kHit->GetDetectorID() % 1000000)
            return -1;
         else if (GetDetectorID() % 1000000 > kHit->GetDetectorID() % 1000000)
            return 1;
      }
   }
   if (TMath::Abs(GetDist()) < TMath::Abs(kHit->GetDist()))
      return 1;
   else if (TMath::Abs(GetDist()) > TMath::Abs(kHit->GetDist()))
      return -1;
   return 0;
}

//__________________________________________________________________________
void MpdKalmanHit::Print(Option_t *opt)
{
   /// Print hit info
}

//__________________________________________________________________________
Double_t MpdKalmanHit::GetPos() const
{
   /// Distance to (0,0,0)

   /*
   if (fHitType != kFixedP) return fDist;
   printf(" !!! Not implemented for kFixedP hits. Exit. \n");
   exit(0);
   */
   return fDist;
}

//__________________________________________________________________________
void MpdKalmanHit::SetIndex(Int_t indx)
{
   /// Add point index

   Int_t size = fIndex.GetSize();
   fIndex.Set(size + 1);
   fIndex[size] = indx;
}

ClassImp(MpdKalmanHit)
