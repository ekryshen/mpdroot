/// \class MpdEctKalmanTrack
///
/// Kalman filter track object for the MPD straw end-cap detector
/// Track parameters: 0: X - X-coordinate @ fixed Z
///                   1: Y - Y-coordinate
///                   2: Phi - local azimuthal angle
///                   3: Theta - local dip angle (angle w.r.t. the transverse plane)
///                   4: q/Pt - signed inverse Pt
/// \author Alexander Zinchenko (LHE, JINR, Dubna)

#include "MpdEctKalmanTrack.h"
#include "FairHit.h"

#include "MpdKalmanFilter.h"
#include "MpdKalmanTrack.h"
//#include "MpdKalmanHitZ.h"
#include "MpdKalmanHit.h"
//#include "TpcLheKalmanTrack.h"
#include "MpdTpcKalmanTrack.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <Riostream.h>

using namespace std;

//__________________________________________________________________________
MpdEctKalmanTrack::MpdEctKalmanTrack()
   : MpdKalmanTrack(), fLastLay(0), fTpcIndex(-1), fTofIndex(-1), fFlag(0), fMirrors(0), fMisses(0), fTrHits(0x0)
{
   /// Default constructor

   // SetType(MpdKalmanTrack::kFixedZ);
   SetType(MpdKalmanTrack::kEndcap);
}

//__________________________________________________________________________
MpdEctKalmanTrack::MpdEctKalmanTrack(Int_t tpcIndx, const MpdTpcKalmanTrack &track)
   : MpdKalmanTrack(track), fLastLay(0), fTpcIndex(tpcIndx), fTofIndex(-1), fFlag(0), fMirrors(0), fMisses(0),
     fTrHits(new TClonesArray("MpdKalmanHit"))
{
   /// Constructor from TPC track

   // MpdKalmanFilter::Instance()->Convert(track, this);
   SetDirection(MpdKalmanTrack::kOutward);
   /*
   SetWeight(*fCovar);
   Int_t iok;
   MpdKalmanFilter::Instance()->MnvertLocal(fWeight->GetMatrixArray(), 5, 5, 5, iok);
   fPos = fPosNew;
   *fParamNew = *fParam;
   */
   SetWeight(*fWeightAtHit);
   fPosNew = fPos = fPosAtHit;
   *fParamNew     = *fParamAtHit;
   *fParam        = *fParamNew;
   fLength        = fLengAtHit;
   // fHits->Clear();
   fNhits = 0;

   Int_t nHits = track.GetNofTrHits();
   for (Int_t i = 0; i < nHits; ++i) fHits->Add(track.GetTrHits()->UncheckedAt(i));
   SetFromTpc();
}

//__________________________________________________________________________
MpdEctKalmanTrack::MpdEctKalmanTrack(Int_t tofIndx, Int_t tpcIndx, FairHit *tof, TpcLheHit *tpc, TVector3 &vert)
   : MpdKalmanTrack(tof->GetZ(), vert), fLastLay(0), fTpcIndex(tpcIndx), fTofIndex(tofIndx), fFlag(0), fMirrors(0),
     fMisses(0),
     // fTrHits(new TClonesArray("MpdKalmanHitZ"))
     fTrHits(new TClonesArray("MpdKalmanHit"))
{
   /// Constructor from ETOF and TPC hits

   // SetTrackID(tof->GetTrackID());
   SetTrackID(-1);
   SetType(MpdKalmanTrack::kEndcap);
   // SetDirection(MpdKalmanTrack::kInward);
}

//__________________________________________________________________________
MpdEctKalmanTrack::MpdEctKalmanTrack(Int_t tofIndx, Int_t ectIndx,
                                     // MpdEtofPoint *tof, MpdKalmanHitZ *ect, TVector3 &vert)
                                     FairHit *tof, MpdKalmanHit *ect, TVector3 &vert)
   : MpdKalmanTrack(tof->GetZ(), vert), fLastLay(-1), fTpcIndex(ectIndx), fTofIndex(tofIndx), fFlag(0), fMirrors(0),
     fMisses(0), fTrHits(new TClonesArray("MpdKalmanHit"))
{
   /// Constructor from ETOF and ECT (CPC) hits

   // Add first hit
   if (ect) {
      fLastLay = ect->GetLayer();
      fHits->Add(ect);
   }

   SetType(MpdKalmanTrack::kEndcap);
   SetDirection(MpdKalmanTrack::kInward);
   SetFromEtof();
}

//__________________________________________________________________________
MpdEctKalmanTrack::MpdEctKalmanTrack(const MpdEctKalmanTrack &track)
   : MpdKalmanTrack(track), fLastLay(track.fLastLay), fTpcIndex(track.fTpcIndex), fTofIndex(track.fTofIndex),
     fFlag(track.fFlag), fMirrors(track.fMirrors), fMisses(track.fMisses),
     // fTrHits(new TClonesArray("MpdKalmanHitZ"))
     fTrHits(new TClonesArray("MpdKalmanHit"))
{
   /// copy constructor

   Int_t nHits = track.fTrHits->GetEntriesFast();
   for (Int_t i = 0; i < nHits; ++i) {
      // MpdKalmanHitZ *hit = (MpdKalmanHitZ*)(track.fTrHits->UncheckedAt(i));
      // new ((*fTrHits)[i]) MpdKalmanHitZ(*hit);
      MpdKalmanHit *hit = (MpdKalmanHit *)(track.fTrHits->UncheckedAt(i));
      new ((*fTrHits)[i]) MpdKalmanHit(*hit);
   }
}

//__________________________________________________________________________
MpdEctKalmanTrack &MpdEctKalmanTrack::operator=(const MpdEctKalmanTrack &track)
{
   /// Asignment operator

   // check assignement to self
   if (this == &track) return *this;

   // base class assignement
   MpdKalmanTrack::operator=(track);

   fLastLay  = track.fLastLay;
   fTpcIndex = track.fTpcIndex;
   fTofIndex = track.fTofIndex;
   fFlag     = track.fFlag;
   fMirrors  = track.fMirrors;
   fMisses   = track.fMisses;
   // fTrHits = new TClonesArray("MpdKalmanHitZ");
   fTrHits = new TClonesArray("MpdKalmanHit");
   // fLengAtHit = track.fLengAtHit;

   Int_t nHits = track.fTrHits->GetEntriesFast();
   for (Int_t i = 0; i < nHits; ++i) {
      // MpdKalmanHitZ *hit = (MpdKalmanHitZ*)(track.fTrHits->UncheckedAt(i));
      // new ((*fTrHits)[i]) MpdKalmanHitZ(*hit);
      MpdKalmanHit *hit = (MpdKalmanHit *)(track.fTrHits->UncheckedAt(i));
      new ((*fTrHits)[i]) MpdKalmanHit(*hit);
   }
   return *this;
}

//__________________________________________________________________________
MpdEctKalmanTrack::~MpdEctKalmanTrack()
{
   /// Destructor

   // if (fTrHits) fTrHits->Clear("C");
   if (fTrHits) fTrHits->Delete();
   delete fTrHits;
}

//__________________________________________________________________________
void MpdEctKalmanTrack::Reset()
{
   /// Reset track

   // if (fTrHits) fTrHits->Clear("C");
   if (fTrHits) fTrHits->Delete();
   delete fTrHits;
   Clear();
}

//__________________________________________________________________________
Double_t MpdEctKalmanTrack::GetRadNew()
{
   /// Get track new radial position

   // if (fTrackType == kFixedZ) {
   if (fTrackType == kEndcap) {
      Double_t x = GetParamNew(0), y = GetParamNew(1);
      Double_t r = x * x + y * y;
      return TMath::Sqrt(r);
   }
   return fPosNew;
}

//__________________________________________________________________________
/*
void TpcLheKalmanTrack::EvalCovar(const MpdKalmanHitR *hitOut, const MpdKalmanHitR *hitIn)
{
  /// Evaluate covariance matrix for track seed

  (*fWeight)(0,0) = hitOut->GetRphiErr() * hitOut->GetRphiErr(); // <RphiRphi>
  (*fWeight)(0,0) *= 4.; // extra factor of 2

  (*fWeight)(1,1) = hitOut->GetZerr() * hitOut->GetZerr(); // <zz>

  Double_t phiOut = hitOut->GetRphi() / hitOut->GetR();
  Double_t phiIn = hitIn->GetRphi() / hitIn->GetR();
  Double_t xIn = hitIn->GetR() * TMath::Cos(phiIn);
  Double_t yIn = hitIn->GetR() * TMath::Sin(phiIn);
  Double_t xOut = hitOut->GetR() * TMath::Cos(phiOut);
  Double_t yOut = hitOut->GetR() * TMath::Sin(phiOut);
  Double_t dist2 = (xOut-xIn)*(xOut-xIn) + (yOut-yIn)*(yOut-yIn);
  Double_t sinPhi = TMath::Sin ((*fParam)(2,0));
  Double_t cosPhi = TMath::Cos ((*fParam)(2,0));
  Double_t pOut = TMath::Cos(phiOut) * cosPhi + TMath::Sin(phiOut) * sinPhi;
  Double_t pIn = TMath::Cos(phiIn) * cosPhi + TMath::Sin(phiIn) * sinPhi;
  (*fWeight)(2,2) = (pOut * pOut + pIn * pIn) / dist2 * (*fWeight)(0,0); // <PhiPhi>
  (*fWeight)(2,2) *= 4.; // extra factor of 2

  Double_t tanThe = TMath::Tan((*fParam)(3,0));
  Double_t dRad = hitOut->GetR() - hitIn->GetR();
  Double_t denom = dRad * (1.+tanThe*tanThe);
  (*fWeight)(3,3) = (*fWeight)(1,1) * 2. / denom / denom; // <TheThe>

  //(*fWeight)(4,4) = ((*fParam)(4,0)*0.5) * ((*fParam)(4,0)*0.5); // error 50%
  //(*fWeight)(4,4) = ((*fParam)(4,0)*0.75) * ((*fParam)(4,0)*0.75); // error 75%
  (*fWeight)(4,4) = ((*fParam)(4,0)*2.) * ((*fParam)(4,0)*2.); // error 200%

  //fWeight->Print();
  //fWeight->Invert(); // weight matrix
  Int_t iok = 0;
  MpdKalmanFilter::Instance()->MnvertLocal(fWeight->GetMatrixArray(), 5, 5, 5, iok);
  //fWeight->Print();
}
*/

//__________________________________________________________________________
void MpdEctKalmanTrack::StartBack()
{
   /// Prepare for back tracing

   fHits->Sort(); // in descending order in radius
   fChi2     = 0.;
   fTrackDir = kInward;
   // if (fLastLay == 50) { cout << " back" << endl; Weight2Cov()->Print(); exit(0); }

   if (GetNofHits() == 1) return;
   /*
   for (Int_t i = 0; i < 5; ++i) {
     for (Int_t j = i; j < 5; ++j) {
       if (j == i) (*fWeight)(i,j) /= 100.;
       else (*fWeight)(i,j) = (*fWeight)(j,i) = 0.;
     }
   }
   */
   fWeight->Zero();
   // EvalCovar((MpdKalmanHitR*)fHits->UncheckedAt(0),(MpdKalmanHitR*)fHits->UncheckedAt(1));
}

//__________________________________________________________________________
void MpdEctKalmanTrack::Print(Option_t *opt)
{
   /// Print track info

   cout << " -------------------------------------------------------------------- " << endl;
   cout << " Track ID: " << GetTrackID() << "; type: " << GetType() << "; direction: " << GetDirection() << endl;
   cout << " Track position: ";
   if (GetType() == MpdKalmanTrack::kBarrel)
      cout << "(R-Phi, Z, R): " << GetParam(0) << " " << GetParam(1) << " " << GetPos() << endl;
   else
      cout << "(X, Y, Z): " << GetParam(0) << " " << GetParam(1) << " " << GetPos() << endl;
   cout << " Track Pt: " << 1. / GetParam(4) << endl;
   cout << " Track new position ";

   if (fParamNew && GetType() == MpdKalmanTrack::kBarrel)
      cout << "(R-Phi, Z, R): " << GetParamNew(0) << " " << GetParamNew(1) << " " << GetPosNew() << endl;
   else
      cout << "(X, Y, Z): " << GetParamNew(0) << " " << GetParamNew(1) << " " << GetPosNew() << endl;
   cout << " Track Pt: " << 1. / GetParamNew(4) << endl;
   cout << " Number of hits: " << GetHits()->GetEntriesFast() << " " << GetNofTrHits() << ", Chi2: " << GetChi2()
        << endl;
   cout << " -------------------------------------------------------------------- " << endl;
}

ClassImp(MpdEctKalmanTrack);
