/// \class MpdKalmanTrack
///
/// Kalman filter track object for the MPD detector
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdKalmanTrack.h"
#include "MpdKalmanFilter.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TObjArray.h>
#include <TVector3.h>
#include <TMath.h>
#include <Riostream.h>

// TpcLheKalmanFilter* TpcLheKalmanTrack::fgKF = 0x0;

inline TObjArray *CreateHitsArray(Int_t size) {
  return size == 0 ? nullptr : new TObjArray(size);
}

inline TMatrixD *CloneMatrix(TMatrixD *matrix, Int_t m, Int_t n) {
  return matrix ? new TMatrixD(*matrix) : new TMatrixD(m, n);
}

inline TMatrixDSym *CloneMatrix(TMatrixDSym *matrix, Int_t n) {
  return matrix ? new TMatrixDSym(*matrix) : new TMatrixDSym(n);
}

//__________________________________________________________________________
MpdKalmanTrack::MpdKalmanTrack(Int_t size)
   : TObject(), fID(0), fNhits(0), fTrackDir(kInward), fTrackType(kBarrel), fLastLay(0), fNofWrong(0), fNode(""),
     fNodeNew(""), fPartID(0), fPos(0.), fPosNew(0.), fPosAtHit(0.), fChi2(0.), fChi2Vertex(0.), fLength(0.),
     fLengAtHit(0.), fParam(0x0), fParamNew(0x0), fParamAtHit(0x0), fCovar(0x0), fWeight(0x0), fWeightAtHit(0x0),
     fHits(CreateHitsArray(size)), fFlag(kOk)
{
   /// Default constructor
}

//__________________________________________________________________________
MpdKalmanTrack::MpdKalmanTrack(Double_t pos, TVector3 &vertex, Int_t size)
   : TObject(), fID(0), fNhits(0), fTrackDir(kInward), fTrackType(kBarrel), fLastLay(0), fNofWrong(0), fNode(""),
     fNodeNew(""), fPartID(0), fPos(pos), fPosNew(pos), fPosAtHit(0.), fChi2(0.), fChi2Vertex(0.), fLength(0.),
     fLengAtHit(0.), fParam(new TMatrixD(5, 1)), fParamNew(new TMatrixD(5, 1)), fParamAtHit(new TMatrixD(5, 1)),
     fCovar(new TMatrixDSym(5)), fWeight(new TMatrixDSym(5)), fWeightAtHit(new TMatrixDSym(5)), fVertex(vertex),
     fHits(CreateHitsArray(size)), fFlag(kOk)
{
   /// Constructor from a track position and vertex

   // fHits->SetOwner(kFALSE);
}

//__________________________________________________________________________
MpdKalmanTrack::MpdKalmanTrack(const MpdKalmanTrack &track)
   : TObject(track),
     fID(track.fID),
     fNhits(track.fNhits),
     fTrackDir(track.fTrackDir),
     fTrackType(track.fTrackType),
     fLastLay(track.fLastLay),
     fNofWrong(track.fNofWrong),
     fNode(track.fNode),
     fNodeNew(track.fNodeNew),
     fPartID(track.fPartID),
     fPos(track.fPos),
     fPosNew(track.fPosNew),
     fPosAtHit(track.fPosAtHit),
     fChi2(track.fChi2),
     fChi2Vertex(track.fChi2Vertex),
     fLength(track.fLength),
     fLengAtHit(track.fLengAtHit),
     fParam(CloneMatrix(track.fParam, 5, 1)),
     fParamNew(CloneMatrix(track.fParamNew, 5, 1)),
     fParamAtHit(CloneMatrix(track.fParamAtHit, 5, 1)),
     fCovar(CloneMatrix(track.fCovar, 5)),
     fWeight(CloneMatrix(track.fWeight, 5)),
     fWeightAtHit(CloneMatrix(track.fWeightAtHit, 5)),
     fVertex(track.fVertex),
     fHits(CreateHitsArray(track.fHits->GetEntriesFast())),
     fStepMap(track.fStepMap),
     fFlag(track.fFlag)
{
   /// copy constructor

   if (track.fHits == 0x0) return;
   Int_t nHits = track.fHits->GetEntriesFast();
   for (Int_t i = 0; i < nHits; ++i) {
      fHits->Add(track.fHits->UncheckedAt(i));
   }
}

//__________________________________________________________________________
MpdKalmanTrack &MpdKalmanTrack::operator=(const MpdKalmanTrack &track)
{
   /// Asignment operator

   // check assignement to self
   if (this == &track) return *this;

   // base class assignement
   TObject::operator=(track);

   Clear();
   fID          = track.fID;
   fNhits       = track.fNhits;
   fTrackDir    = track.fTrackDir;
   fTrackType   = track.fTrackType;
   fLastLay     = track.fLastLay;
   fNofWrong    = track.fNofWrong;
   fNode        = track.fNode;
   fNodeNew     = track.fNodeNew;
   fPartID      = track.fPartID;
   fPos         = track.fPos;
   fPosNew      = track.fPosNew;
   fPosAtHit    = track.fPosAtHit;
   fChi2        = track.fChi2;
   fChi2Vertex  = track.fChi2Vertex;
   fLength      = track.fLength;
   fLengAtHit   = track.fLengAtHit;
   fParam       = CloneMatrix(track.fParam, 5, 1);
   fParamNew    = CloneMatrix(track.fParamNew, 5, 1);
   fParamAtHit  = CloneMatrix(track.fParamAtHit, 5, 1);
   fCovar       = CloneMatrix(track.fCovar, 5);
   fWeight      = CloneMatrix(track.fWeight, 5);
   fWeightAtHit = CloneMatrix(track.fWeightAtHit, 5);
   fVertex      = track.fVertex;
// fHits        = CreateHitsArray(track.fHits->GetEntriesFast());
   fHits        = CreateHitsArray(70);
   fStepMap     = track.fStepMap;
   fFlag        = track.fFlag;

   if (track.fHits == 0x0) return *this;
   Int_t nHits = track.fHits->GetEntriesFast();
   for (Int_t i = 0; i < nHits; ++i) {
      fHits->Add(track.fHits->UncheckedAt(i));
   }
   return *this;
}

//__________________________________________________________________________
MpdKalmanTrack::~MpdKalmanTrack()
{
   /// Destructor
   delete fParam;
   delete fParamNew;
   delete fParamAtHit;
   fParam = fParamNew = fParamAtHit = NULL;
   delete fCovar;
   fCovar = NULL;
   delete fWeight;
   delete fWeightAtHit;
   fWeight = fWeightAtHit = NULL;
   delete fHits;
   fHits = NULL;
}

//__________________________________________________________________________
void MpdKalmanTrack::Clear()
{
   /// Reset track
   delete fParam;
   delete fParamNew;
   delete fParamAtHit;
   fParam = fParamNew = fParamAtHit = NULL;
   delete fCovar;
   fCovar = NULL;
   delete fWeight;
   delete fWeightAtHit;
   fWeight = fWeightAtHit = NULL;
   delete fHits;
   fHits = NULL;
   fStepMap.clear();
}

//__________________________________________________________________________
TMatrixDSym *MpdKalmanTrack::Weight2Cov()
{
   /// Get covariance matrix from weight matrix
   *fCovar = *fWeight;
   // fCovar->Invert();
   Int_t iok = 0;
   MpdKalmanFilter::Instance()->MnvertLocal(fCovar->GetMatrixArray(), 5, 5, 5, iok);
   return fCovar;
}

//__________________________________________________________________________
void MpdKalmanTrack::ReSetWeight()
{
   /// Set weight matrix from covariance
   *fWeight  = *fCovar;
   Int_t iok = 0;
   MpdKalmanFilter::Instance()->MnvertLocal(fWeight->GetMatrixArray(), 5, 5, 5, iok);
}

//__________________________________________________________________________
Double_t MpdKalmanTrack::DCA() const
{
   /// Return signed track DCA

   if (fPosNew < 1.e-6) return fPosNew;
   Float_t  phiTr = GetParam(2);
   TVector3 vecTr(TMath::Cos(phiTr), TMath::Sin(phiTr), 0.0);
   Float_t  phi = GetParam(0) / fPosNew;
   TVector3 vec(TMath::Cos(phi), TMath::Sin(phi), 0.0);
   return fPosNew * TMath::Sign(1.0, vec.Cross(vecTr).Z());
}

//__________________________________________________________________________
void MpdKalmanTrack::PrintNode(Bool_t nNew) const
{
   /// Print node name (for debugging)

   if (nNew)
      printf("%s\n", fNodeNew.Data());
   else
      printf("%s\n", fNode.Data());
}

///__________________________________________________________________________

ClassImp(MpdKalmanTrack);
