/// \class MpdItsKalmanTrack
///
/// Kalman filter track object for the MPD central detector
/// Track parameters: 0: RPhi - coordinate in R-Phi direction
///                   1: Z - longitudinal coordinate
///                   2: Phi - local azimuthal angle
///                   3: Theta - local dip angle (angle w.r.t. the transverse plane)
///                   4: q/Pt - signed inverse Pt
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdItsKalmanTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
#include "MpdVector.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TVector3.h>
#include <TMath.h>
#include <Riostream.h>

//__________________________________________________________________________
MpdItsKalmanTrack::MpdItsKalmanTrack()
   : MpdKalmanTrack(), fTrHits(0x0),
     // fKHits(0x0),//attention!
     fNofIts(0), fChi2Its(0.0)
{
   /// Default constructor
}

//__________________________________________________________________________
MpdItsKalmanTrack::MpdItsKalmanTrack(MpdKalmanHit *hitOut, MpdKalmanHit *hitIn, TVector3 &vertex, Double_t pt)
   : MpdKalmanTrack(0., vertex), fTrHits(new TClonesArray("MpdKalmanHit")),
     // fKHits(new TClonesArray("MpdKalmanHit")), //attention!!!
     fNofIts(0), fChi2Its(0.0)
{
   /// Constructor from 2 hits

   Double_t rOut = 0, phiOut = 0, rIn = 0, phiIn = 0;
   if (hitIn->GetType() == MpdKalmanHit::kFixedP) {
      TVector3 posIn  = MpdKalmanFilter::Instance()->GetGeo()->GlobalPos(hitIn);
      TVector3 posOut = MpdKalmanFilter::Instance()->GetGeo()->GlobalPos(hitOut);
      rOut            = posOut.Pt();
      phiOut          = posOut.Phi();
      rIn             = posIn.Pt();
      phiIn           = posIn.Phi();
   } else {
      rOut   = hitOut->GetPos();
      phiOut = hitOut->GetMeas(0) / rOut;
      rIn    = hitIn->GetPos();
      phiIn  = hitIn->GetMeas(0) / rIn;
   }

   Double_t xIn       = rIn * TMath::Cos(phiIn);
   Double_t yIn       = rIn * TMath::Sin(phiIn);
   Double_t xOut      = rOut * TMath::Cos(phiOut);
   Double_t yOut      = rOut * TMath::Sin(phiOut);
   Double_t parOut[4] = {rOut, phiOut, xOut, yOut};
   Double_t parIn[4]  = {rIn, phiIn, xIn, yIn};

   fPos = fPosNew  = rOut;
   (*fParam)(0, 0) = rOut * phiOut; // R-Phi coordinate
   //(*fParam)(1,0) = posOut.Z(); // Z-coordinate
   (*fParam)(1, 0) = hitOut->GetMeas(1);                   // Z-coordinate
   (*fParam)(2, 0) = TMath::ATan2(yOut - yIn, xOut - xIn); // Track Phi
   //(*fParam)(3,0) = TMath::ATan2 (posOut.Z()-posIn.Z(), rOut-rIn); // Track Theta
   (*fParam)(3, 0) = TMath::ATan2(hitOut->GetMeas(1) - hitIn->GetMeas(1), rOut - rIn); // Track Theta
   (*fParam)(4, 0) = 1. / pt;                                                          // q/Pt

   *fParamNew = *fParam;
   fWeight->Zero();
   EvalCovar(hitOut, hitIn, parOut, parIn);

   fHits->Add(hitOut);
   fLastLay = hitOut->GetLayer(); // layer number
}

//__________________________________________________________________________
MpdItsKalmanTrack::MpdItsKalmanTrack(const MpdItsKalmanTrack &track)
   : MpdKalmanTrack(track), fTrHits(new TClonesArray("MpdKalmanHit")), fNofIts(track.fNofIts), fChi2Its(track.fChi2Its)
{
   /// copy constructor

   Int_t nHits = track.fTrHits->GetEntriesFast();
   for (Int_t i = 0; i < nHits; ++i) {
      MpdKalmanHit *hit = (MpdKalmanHit *)(track.fTrHits->UncheckedAt(i));
      new ((*fTrHits)[i]) MpdKalmanHit(*hit);
   }
}

//__________________________________________________________________________
MpdItsKalmanTrack &MpdItsKalmanTrack::operator=(const MpdItsKalmanTrack &track)
{
   /// Asignment operator

   // check assignement to self
   if (this == &track) return *this;

   // base class assignement
   MpdKalmanTrack::operator=(track);

   if (fTrHits) {
      fTrHits->Delete();
      delete fTrHits;
   }
   fTrHits  = new TClonesArray("MpdKalmanHit");
   fNofIts  = track.fNofIts;
   fChi2Its = track.fChi2Its;

   Int_t nHits = track.fTrHits->GetEntriesFast();
   for (Int_t i = 0; i < nHits; ++i) {
      MpdKalmanHit *hit = (MpdKalmanHit *)(track.fTrHits->UncheckedAt(i));
      new ((*fTrHits)[i]) MpdKalmanHit(*hit);
   }
   if (track.fHits == 0x0) {
      for (Int_t i = 0; i < nHits; ++i) {
         fHits->Add(fTrHits->UncheckedAt(i));
      }
   }
   return *this;
}

//__________________________________________________________________________
MpdItsKalmanTrack::MpdItsKalmanTrack(const MpdTpcKalmanTrack &track)
   : MpdKalmanTrack(track), fTrHits(new TClonesArray("MpdKalmanHit")), fNofIts(0), fChi2Its(0.0)
{
   /// constructor from TPC track

   TClonesArray *hits  = track.GetTrHits();
   Int_t         nHits = hits->GetEntriesFast();
   for (Int_t i = 0; i < nHits; ++i) {
      MpdKalmanHit *hit = (MpdKalmanHit *)(hits->UncheckedAt(i));
      new ((*fTrHits)[i]) MpdKalmanHit(*hit);
   }
}

//__________________________________________________________________________
MpdItsKalmanTrack::MpdItsKalmanTrack(const MpdVector &track, TVector3 &vec)
   : MpdKalmanTrack(0.0, vec), fTrHits(new TClonesArray("MpdKalmanHit")), fNofIts(0), fChi2Its(0.0)
{
   /// construct from last layer celltrack
   const MpdVector *temp = &track;

   for (Int_t i = 0; i < 5; ++i) {
      MpdKalmanHit *hit = temp->GetKalmanHit();
      temp              = temp->GetPrevTrackPointer();
      new ((*fTrHits)[i]) MpdKalmanHit(*hit);
      if (temp == NULL) break;
   }
   fTrHits->Sort();
}

//__________________________________________________________________________
MpdItsKalmanTrack::~MpdItsKalmanTrack()
{
   /// Destructor

   /// if (fTrHits) fTrHits->Clear("C");
   if (fTrHits) fTrHits->Delete(); // AZ
   delete fTrHits;
   fTrHits = NULL;
}

//__________________________________________________________________________
void MpdItsKalmanTrack::Reset()
{
   /// Reset track

   // if (fTrHits) fTrHits->Clear("C");
   if (fTrHits) fTrHits->Delete();
   delete fTrHits;
   fTrHits = NULL;
   Clear();
}

//__________________________________________________________________________
void MpdItsKalmanTrack::EvalCovar(const MpdKalmanHit *hitOut, const MpdKalmanHit *hitIn, Double_t *parOut,
                                  Double_t *parIn)
{
   /// Evaluate covariance matrix for track seed

   (*fWeight)(0, 0) = hitOut->GetErr(0) * hitOut->GetErr(0); // <RphiRphi>
   (*fWeight)(0, 0) *= 4.;                                   // extra factor of 2

   (*fWeight)(1, 1) = hitOut->GetErr(1) * hitOut->GetErr(1); // <zz>

   Double_t phiOut = parOut[1], phiIn = parIn[1];
   Double_t dx = parOut[2] - parIn[2], dy = parOut[3] - parIn[3];

   Double_t dist2   = dx * dx + dy * dy;
   Double_t sinPhi  = TMath::Sin((*fParam)(2, 0));
   Double_t cosPhi  = TMath::Cos((*fParam)(2, 0));
   Double_t pOut    = TMath::Cos(phiOut) * cosPhi + TMath::Sin(phiOut) * sinPhi;
   Double_t pIn     = TMath::Cos(phiIn) * cosPhi + TMath::Sin(phiIn) * sinPhi;
   (*fWeight)(2, 2) = (pOut * pOut + pIn * pIn) / dist2 * (*fWeight)(0, 0); // <PhiPhi>
   (*fWeight)(2, 2) *= 4.;                                                  // extra factor of 2

   Double_t tanThe  = TMath::Tan((*fParam)(3, 0));
   Double_t dRad    = parOut[0] - parIn[0];
   Double_t denom   = dRad * (1. + tanThe * tanThe);
   (*fWeight)(3, 3) = (*fWeight)(1, 1) * 2. / denom / denom; // <TheThe>
   //(*fWeight)(3,3) = (*fWeight)(1,1) * 20. / denom / denom; // <TheThe>

   //(*fWeight)(4,4) = ((*fParam)(4,0)*0.5) * ((*fParam)(4,0)*0.5); // error 50%
   //(*fWeight)(4,4) = ((*fParam)(4,0)*0.75) * ((*fParam)(4,0)*0.75); // error 75%
   (*fWeight)(4, 4) = ((*fParam)(4, 0) * 2.) * ((*fParam)(4, 0) * 2.); // error 200%

   // fWeight->Print();
   // fWeight->Invert(); // weight matrix
   Int_t iok = 0;
   MpdKalmanFilter::Instance()->MnvertLocal(fWeight->GetMatrixArray(), 5, 5, 5, iok);
   // fWeight->Print();
}

//__________________________________________________________________________
void MpdItsKalmanTrack::StartBack()
{
   /// Prepare for back tracing

   fHits->Sort(); // in descending order in radius
   fChi2     = 0.;
   fTrackDir = kInward;
   // if (fLastLay == 50) { cout << " back" << endl; Weight2Cov()->Print(); exit(0); }

   Int_t nHits = GetNofHits();
   if (nHits == 1) return;
   nHits *= nHits;
   //*
   for (Int_t i = 0; i < 5; ++i) {
      for (Int_t j = i; j < 5; ++j) {
         // if (j == i) (*fWeight)(i,j) /= 100.;
         if (j == i)
            (*fWeight)(i, j) /= nHits;
         else
            (*fWeight)(i, j) = (*fWeight)(j, i) = 0.;
      }
   }
   //*/
   // fWeight->Zero();
   // EvalCovar((MpdKalmanHitR*)fHits->UncheckedAt(0),(MpdKalmanHitR*)fHits->UncheckedAt(1));
}

//__________________________________________________________________________
Int_t MpdItsKalmanTrack::Compare(const TObject *track) const
{
   /// "Compare" function to sort in descending order in pt

   MpdKalmanTrack *trackKF = (MpdKalmanTrack *)track;
   Double_t        pt      = 1. / TMath::Abs(trackKF->GetParam(4));
   Double_t        ptthis  = 1. / TMath::Abs((*fParam)(4, 0));
   if (ptthis < pt)
      return 1;
   else if (ptthis > pt)
      return -1;
   return 0;
}

ClassImp(MpdItsKalmanTrack);
