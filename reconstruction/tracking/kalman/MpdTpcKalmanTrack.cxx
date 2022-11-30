/// \class MpdTpcKalmanTrack
///
/// Kalman filter track object for the MPD central detector
/// Track parameters: 0: RPhi - coordinate in R-Phi direction
///                   1: Z - longitudinal coordinate
///                   2: Phi - local azimuthal angle
///                   3: Theta - local dip angle (angle w.r.t. the transverse plane)
///                   4: q/Pt - signed inverse Pt
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
//#include "MpdTpcHit.h" // debug
//#include "TpcPoint.h" // debug

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TVector3.h>
#include <TMath.h>
#include <Riostream.h>

using namespace std;

inline TClonesArray *CreateTrHitsArray(Int_t size) {
  return size == 0 ? nullptr : new TClonesArray("MpdKalmanHit", size);
}

//__________________________________________________________________________
MpdTpcKalmanTrack::MpdTpcKalmanTrack(Int_t size)
    : MpdKalmanTrack(size), fTrHits(CreateTrHitsArray(size))
{
   /// Default constructor
}

//__________________________________________________________________________
MpdTpcKalmanTrack::MpdTpcKalmanTrack(MpdKalmanHit *hitOut, MpdKalmanHit *hitIn, TVector3 &vertex, Double_t pt,
                                     Double_t *params, Int_t size)
   : MpdKalmanTrack(0., vertex), fTrHits(CreateTrHitsArray(size))
{
   /// Constructor from 2 hits

   Double_t rOut      = hitOut->GetPos();
   Double_t phiOut    = hitOut->GetMeas(0) / rOut;
   Double_t rIn       = hitIn->GetPos();
   Double_t phiIn     = hitIn->GetMeas(0) / rIn;
   Double_t parOut[4] = {rOut, phiOut, 0., 0.};
   Double_t parIn[4]  = {rIn, phiIn, 0., 0.};
   EvalParams(hitOut, hitIn, parOut, parIn, pt, params);
}

//__________________________________________________________________________
MpdTpcKalmanTrack::MpdTpcKalmanTrack(MpdKalmanHit *hitOut, MpdKalmanHit *hitIn, TVector3 &vertex, TVector3 &posOut,
                                     TVector3 &posIn, Double_t pt, Double_t *params, Int_t size)
   : MpdKalmanTrack(0., vertex), fTrHits(CreateTrHitsArray(size))
{
   /// Constructor from 2 hits and 2 TVector3's (modular geometry)

   Double_t rOut = posOut.Pt(), phiOut = posOut.Phi();
   Double_t rIn = posIn.Pt(), phiIn = posIn.Phi();
   Double_t parOut[4] = {rOut, phiOut, 0., 0.};
   Double_t parIn[4]  = {rIn, phiIn, 0., 0.};
   EvalParams(hitOut, hitIn, parOut, parIn, pt, params);
}

//__________________________________________________________________________
MpdTpcKalmanTrack::MpdTpcKalmanTrack(const MpdTpcKalmanTrack &track)
   : MpdKalmanTrack(track), fTrHits(CreateTrHitsArray(track.fTrHits->GetEntriesFast()))
{
   /// copy constructor

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
}

//__________________________________________________________________________
MpdTpcKalmanTrack &MpdTpcKalmanTrack::operator=(const MpdTpcKalmanTrack &track)
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

   Int_t nHits = track.fTrHits->GetEntriesFast();
   fTrHits = CreateTrHitsArray(nHits);

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
MpdTpcKalmanTrack::~MpdTpcKalmanTrack()
{
   /// Destructor

   if (fTrHits) fTrHits->Delete();
   delete fTrHits;
   fTrHits = NULL;
}

//__________________________________________________________________________
void MpdTpcKalmanTrack::Reset()
{
   /// Reset track

   if (fTrHits) fTrHits->Delete();
   delete fTrHits;
   fTrHits = NULL;
   Clear();
}

//__________________________________________________________________________
void MpdTpcKalmanTrack::EvalParams(MpdKalmanHit *hitOut, const MpdKalmanHit *hitIn, Double_t *parOut, Double_t *parIn,
                                   Double_t pt, Double_t *params)
{
   /// Evaluate track params

   Double_t rIn = parIn[0], phiIn = parIn[1];
   Double_t rOut = parOut[0], phiOut = parOut[1];
   Double_t xIn  = rIn * TMath::Cos(phiIn);
   Double_t yIn  = rIn * TMath::Sin(phiIn);
   Double_t xOut = rOut * TMath::Cos(phiOut);
   Double_t yOut = rOut * TMath::Sin(phiOut);
   parOut[2]     = xOut;
   parOut[3]     = yOut;
   parIn[2]      = xIn;
   parIn[3]      = yIn;

   fPos = fPosNew  = rOut;
   (*fParam)(0, 0) = rOut * phiOut; // R-Phi coordinate
   //(*fParam)(1,0) = posOut.Z(); // Z-coordinate
   (*fParam)(1, 0) = hitOut->GetMeas(1);                   // Z-coordinate
   Double_t phi0   = TMath::ATan2(yOut - yIn, xOut - xIn); // Track Phi - approximation
   (*fParam)(2, 0) = phi0 + TMath::Sign(params[2] / 2, pt);
   //(*fParam)(3,0) = TMath::ATan2 (hitOut->GetMeas(1)-hitIn->GetMeas(1), rOut-rIn); // Track Theta
   (*fParam)(3, 0) = TMath::ATan2(hitOut->GetMeas(1) - hitIn->GetMeas(1), TMath::Abs(params[0])); // Track Theta
   (*fParam)(4, 0) = 1. / pt;                                                                     // q/Pt
   // Debugging
   /*
   TClonesArray *hits = (TClonesArray*) FairRootManager::Instance()->GetObject("TpcHit");
   TClonesArray *points = (TClonesArray*) FairRootManager::Instance()->GetObject("TpcPoint");
   MpdTpcHit *h = (MpdTpcHit*) hits->UncheckedAt(hitOut->GetIndex());
   TpcPoint *p = (TpcPoint*) points->UncheckedAt(h->GetRefIndex());
   cout << " Phi angle: " << TMath::ATan2(p->GetPy(),p->GetPx()) << " " << phi0 << " " << (*fParam)(2,0) << endl;
   */

   *fParamNew = *fParam;
   fWeight->Zero();
   EvalCovar(hitOut, hitIn, parOut, parIn);

   fHits->Add(hitOut);
   fLastLay = hitOut->GetLayer(); // layer number
}

//__________________________________________________________________________
void MpdTpcKalmanTrack::EvalCovar(const MpdKalmanHit *hitOut, const MpdKalmanHit *hitIn, Double_t *parOut,
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

   Double_t tanThe = TMath::Tan((*fParam)(3, 0));
   // 31.12.2018 Double_t dRad = parOut[0] - parIn[0];
   Double_t dRad    = TMath::Sqrt(dist2);
   Double_t denom   = dRad * (1. + tanThe * tanThe);
   (*fWeight)(3, 3) = (*fWeight)(1, 1) * 2. / denom / denom; // <TheThe>
   //(*fWeight)(3,3) = (*fWeight)(1,1) * 20. / denom / denom; // <TheThe>

   //(*fWeight)(4,4) = ((*fParam)(4,0)*0.5) * ((*fParam)(4,0)*0.5); // error 50%
   //(*fWeight)(4,4) = ((*fParam)(4,0)*0.75) * ((*fParam)(4,0)*0.75); // error 75%
   (*fWeight)(4, 4) = ((*fParam)(4, 0) * 2.) * ((*fParam)(4, 0) * 2.); // error 200%

   // for (Int_t jj = 0; jj < 5; ++jj) cout << TMath::Sqrt((*fWeight)(jj,jj)) << " ";
   // cout << endl;

   // fWeight->Print();
   // fWeight->Invert(); // weight matrix
   Int_t iok = 0;
   MpdKalmanFilter::Instance()->MnvertLocal(fWeight->GetMatrixArray(), 5, 5, 5, iok);
   // fWeight->Print();
}

//__________________________________________________________________________
void MpdTpcKalmanTrack::StartBack()
{
   /// Prepare for back tracing

   //#pragma omp critical
   {
      fHits->Sort(); // in descending order in radius
   }
   fChi2 = 0.;
   // fTrackDir = kInward;
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
Int_t MpdTpcKalmanTrack::Compare(const TObject *track) const
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
//__________________________________________________________________________
Bool_t MpdTpcKalmanTrack::GetRecoQuality(Double_t dist, Double_t percentage)
{
   /// returns kTRUE if number of hits closer to boundaries than dist divided by nHits is larger than percentage

   MpdKalmanHit *hit        = NULL;
   Int_t         nCloseHits = 0;
   Int_t         nTrHits    = fTrHits->GetEntries();
   if (nTrHits == 0) {
      cout << "NO KALMAN HITS ARE FOUND" << endl;
      return kTRUE;
   }
   for (Int_t itr = 0; itr < nTrHits; itr++) {
      hit = (MpdKalmanHit *)fTrHits->UncheckedAt(itr);
      if (TMath::Abs(hit->GetEdge()) < dist) nCloseHits++;
   }
   if ((Double_t)(nCloseHits / (Double_t)nTrHits) > percentage)
      return kTRUE;
   else
      return kFALSE;
}

ClassImp(MpdTpcKalmanTrack);
