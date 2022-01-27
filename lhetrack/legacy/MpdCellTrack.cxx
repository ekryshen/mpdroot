/// \class MpdCellTrack
/// 
/// Cellular automaton track object for the MPD inner tracking system
/// \author Alexander Zinchenko, Maxim Strelchenko (LHEP, JINR, Dubna)

#include "MpdCellTrack.h" 

#include <map>

//__________________________________________________________________________
MpdCellTrack::MpdCellTrack() 
  : TObject(),
    fDetectorID(0),
    fFlag(1),
    fLength(0.),
    fIndex(0),
    fHitType(kFixedP),
    fNofDim(2),
    fSignal(0.),
    fDist(0.),
    fCode("")
 
{
  /// Default constructor
  
  //for (Int_t i = 0; i < 2; ++i) fMeas[i] = fErr[i] = fCosSin[i] = 0;
  // fCosSin[0] = 1.;
}
 
//__________________________________________________________________________
MpdCellTrack::MpdCellTrack(Int_t detID, Int_t nDim, HitType hitType, TVector3 &meas, Double_t *err, 
			   Double_t *cossin, Double_t signal, Double_t dist, 
			   Int_t index, Int_t index1, Int_t trackNo)
  : TObject(),
    fDetectorID(detID),
    fFlag(1),
    fLength(0.),
    fHitType(hitType),
    fNofDim(nDim),
    fSignal(signal),
    fDist(dist),
    fMeas(meas),
    fTrackNo(trackNo),
    fCode("")
{
  /// Constructor

  for (Int_t i = 0; i < 2; ++i) {
    //fMeas[i] = meas[i];
    fErr[i] = err[i];
    fCosSin[i] = cossin[i];
  } 
  SetIndex(index);
  SetIndex(index1);
}

//__________________________________________________________________________
MpdCellTrack::MpdCellTrack (const MpdCellTrack& track)
  : TObject(track),
    fDetectorID(track.fDetectorID),
    fFlag(track.fFlag),
    fLength(track.fLength),
    fIndex(track.fIndex),
    fHitType(track.fHitType),
    fNofDim(track.fNofDim),
    fSignal(track.fSignal),
    fDist(track.fDist),
    fMeas(track.fMeas),
    fTrackNo(track.fTrackNo),
    fCode(track.fCode)
{
  ///copy constructor

  for (Int_t i = 0; i < 2; ++i) {
    //fMeas[i] = hit.fMeas[i];
    fErr[i] = track.fErr[i];
    fCosSin[i] = track.fCosSin[i];
   }
}

//__________________________________________________________________________
MpdCellTrack::~MpdCellTrack() 
{
  /// Destructor
}

//__________________________________________________________________________
Int_t MpdCellTrack::Compare(const TObject* hit) const
{
  /// "Compare" function to sort in descending order in fDist

  MpdCellTrack *kHit = (MpdCellTrack*) hit;
  if (kHit->GetType() == GetType()) {
    // Check layers
    if (GetLayer() < kHit->GetLayer()) return 1;
    else if (GetLayer() > kHit->GetLayer()) return -1;
    if (GetType() == kFixedP) {
      // Sort according to sector number
      if (GetDetectorID() % 1000000 < kHit->GetDetectorID() % 1000000) return -1;
      else if (GetDetectorID() % 1000000 > kHit->GetDetectorID() % 1000000) return 1;
    }
  }
  if (TMath::Abs(GetDist()) < TMath::Abs(kHit->GetDist())) return 1;
  else if (TMath::Abs(GetDist()) > TMath::Abs(kHit->GetDist())) return -1;
  return 0;
}

//__________________________________________________________________________
void MpdCellTrack::Print(Option_t *opt)
{
  /// Print hit info

}

//__________________________________________________________________________
Double_t MpdCellTrack::GetPos() const
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
void MpdCellTrack::SetIndex(Int_t indx)
{
  /// Add point index

  Int_t size = fIndex.GetSize();
  fIndex.Set (size + 1);
  fIndex[size] = indx;
}

ClassImp(MpdCellTrack);

