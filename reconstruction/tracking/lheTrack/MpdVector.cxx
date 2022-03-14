/// \class MpdVector
///
/// Cellular automaton track object for the MPD inner tracking system
/// \author Alexander Zinchenko, Maxim Strelchenko (LHEP, JINR, Dubna)

#include "MpdVector.h"

#include <map>

//__________________________________________________________________________
MpdVector::MpdVector()
   : TObject(), fDetectorID(0), fFlag(1), fIndex(0), fHitType(kFixedP), fNofDim(2), fLength(0.), fSignal(0.), fDist(0.),
     fCode(""), fTrackPointer(nullptr), fKalmanHit(nullptr)

{
   /// Default constructor

   for (Int_t i = 0; i < 2; ++i) fMeas[i] = fErr[i] = fCosSin[i] = 0;
   // fCosSin[0] = 1.;
}

//__________________________________________________________________________
MpdVector::MpdVector(Int_t detID, Int_t nDim, HitType hitType, TVector3 &meas, Double_t *err, Double_t *cossin,
                     Double_t signal, Double_t dist, Int_t index, Int_t index1, Int_t trackNo, MpdVector *pointer,
                     MpdKalmanHit *hit)
   : TObject(), fDetectorID(detID), fFlag(1), fHitType(hitType), fNofDim(nDim), fLength(0.), fSignal(signal),
     fDist(dist), fMeas(meas), fTrackNo(trackNo), fCode(""), fTrackPointer(pointer), fKalmanHit(hit)
{
   /// Constructor

   for (Int_t i = 0; i < 2; ++i) {
      // fMeas[i] = meas[i];
      fErr[i]    = err[i];
      fCosSin[i] = cossin[i];
   }
   SetIndex(index);
   SetIndex(index1);
}

//__________________________________________________________________________
MpdVector::MpdVector(const MpdVector &track)
   : TObject(track), fDetectorID(track.fDetectorID), fFlag(track.fFlag), fLength(track.fLength), fIndex(track.fIndex),
     fHitType(track.fHitType), fNofDim(track.fNofDim), fSignal(track.fSignal), fDist(track.fDist), fMeas(track.fMeas),
     fTrackNo(track.fTrackNo), fCode(track.fCode), fTrackPointer(track.fTrackPointer), fKalmanHit(track.fKalmanHit)
{
   /// copy constructor

   for (Int_t i = 0; i < 2; ++i) {
      // fMeas[i] = hit.fMeas[i];
      fErr[i]    = track.fErr[i];
      fCosSin[i] = track.fCosSin[i];
   }
}

//__________________________________________________________________________
MpdVector::~MpdVector()
{
   /// Destructor
}

//__________________________________________________________________________
Int_t MpdVector::Compare(const TObject *hit) const
{
   /// "Compare" function to sort in descending order in fDist

   MpdVector *kHit = (MpdVector *)hit;
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
void MpdVector::Print(Option_t *opt)
{
   /// Print hit info
}

//__________________________________________________________________________
Double_t MpdVector::GetPos() const
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
void MpdVector::SetIndex(Int_t indx)
{
   /// Add point index

   Int_t size = fIndex.GetSize();
   fIndex.Set(size + 1);
   fIndex[size] = indx;
}

ClassImp(MpdVector);
