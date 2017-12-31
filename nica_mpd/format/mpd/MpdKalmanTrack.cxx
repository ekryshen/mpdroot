/// \class MpdKalmanTrack
/// 
/// Kalman filter track object for the MPD detector
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdKalmanTrack.h"
//#include "MpdKalmanFilter.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TObjArray.h>
#include <TVector3.h>
#include <TMath.h>
#include <Riostream.h>

//TpcLheKalmanFilter* TpcLheKalmanTrack::fgKF = 0x0;

//__________________________________________________________________________
MpdKalmanTrack::MpdKalmanTrack() 
  : TObject(),
    fID(0),
    fNhits(0),
    fTrackDir(kInward),
    fTrackType(kBarrel),
    fLastLay(0),
    fNofWrong(0),
    fNode(""),
    fNodeNew(""),
    fPartID(0),
    fPos(0.),
    fPosNew(0.),
    fPosAtHit(0.),
    fChi2(0.),
    fChi2Vertex(0.),
    fLength(0.),
    fLengAtHit(0.),
    fParam(0x0),
    fParamNew(0x0),
    fParamAtHit(0x0),
    fCovar(0x0),
    fWeight(0x0),
    fWeightAtHit(0x0),
    fHits(0x0)
{
  /// Default constructor
}

//__________________________________________________________________________
MpdKalmanTrack::MpdKalmanTrack(Double_t pos, TVector3 &vertex) 
  : TObject(),
    fID(0),
    fNhits(0),
    fTrackDir(kInward),
    fTrackType(kBarrel),
    fLastLay(0),
    fNofWrong(0),
    fNode(""),
    fNodeNew(""),
    fPartID(0),
    fPos(pos),
    fPosNew(pos),
    fPosAtHit(0.),
    fChi2(0.),
    fChi2Vertex(0.),
    fLength(0.),
    fLengAtHit(0.),
    fParam(new TMatrixD(5,1)),
    fParamNew(new TMatrixD(5,1)),
    fParamAtHit(new TMatrixD(5,1)),
    fCovar(new TMatrixDSym(5)),
    fWeight(new TMatrixDSym(5)),
    fWeightAtHit(new TMatrixDSym(5)),
    fVertex(vertex),
    fHits(new TObjArray(70))
{
  /// Constructor from a track position and vertex

  //fHits->SetOwner(kFALSE);
}

//__________________________________________________________________________
MpdKalmanTrack::MpdKalmanTrack (const MpdKalmanTrack& track)
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
    fParam(new TMatrixD(*(track.fParam))),
    //fParamNew(new TMatrixD(*(track.fParamNew))),
    fParamAtHit(new TMatrixD(*(track.fParamAtHit))),
    fCovar(new TMatrixDSym(*(track.fCovar))),
    //fWeight(new TMatrixDSym(*(track.fWeight))),
    fWeightAtHit(new TMatrixDSym(*(track.fWeightAtHit))),
    fVertex(track.fVertex),
    fHits(new TObjArray(70)),
    fStepMap(track.fStepMap)
{
  ///copy constructor

  if (track.fParamNew) fParamNew = new TMatrixD(*(track.fParamNew));
  else fParamNew = new TMatrixD(5,1);
  if (track.fWeight) fWeight = new TMatrixDSym(*(track.fWeight));
  else fWeight = new TMatrixDSym(5);

  if (track.fHits == 0x0) return;
  Int_t nHits = track.fHits->GetEntriesFast();
  for (Int_t i = 0; i < nHits; ++i) {
    fHits->Add(track.fHits->UncheckedAt(i));
  }
}

//__________________________________________________________________________
MpdKalmanTrack & MpdKalmanTrack::operator=(const MpdKalmanTrack& track)
{
  /// Asignment operator

  // check assignement to self
  if (this == &track) return *this;

  // base class assignement
  TObject::operator=(track);

  fID = track.fID;
  fNhits = track.fNhits;
  fTrackDir = track.fTrackDir;
  fTrackType = track.fTrackType;
  fLastLay = track.fLastLay;
  fNofWrong = track.fNofWrong;
  fNode = track.fNode;
  fNodeNew = track.fNodeNew;
  fPartID = track.fPartID;
  fPos = track.fPos;
  fPosNew = track.fPosNew;
  fPosAtHit = track.fPosAtHit;
  fChi2 = track.fChi2;
  fChi2Vertex = track.fChi2Vertex;
  fLength = track.fLength;
  fLengAtHit = track.fLengAtHit;
  if(fParam) delete fParam;
  fParam = new TMatrixD(*(track.fParam));
  if(fParamAtHit) delete fParamAtHit;
  fParamAtHit = new TMatrixD(*(track.fParamAtHit));
  if(fCovar) delete fCovar;
  fCovar = new TMatrixDSym(*(track.fCovar));
  if(fWeightAtHit) delete fWeightAtHit;
  fWeightAtHit = new TMatrixDSym(*(track.fWeightAtHit));
  fVertex = track.fVertex;
  fHits = new TObjArray(70);
  fStepMap = track.fStepMap;
  if(fParamNew) delete fParamNew;
  if (track.fParamNew) fParamNew = new TMatrixD(*(track.fParamNew));
  else fParamNew = new TMatrixD(5,1);
  if(fWeight) delete fWeight;
  if (track.fWeight) fWeight = new TMatrixDSym(*(track.fWeight));
  else fWeight = new TMatrixDSym(5);

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
  delete fCovar; 
  delete fWeight; 
  delete fWeightAtHit; 
  delete fHits; 
}

//__________________________________________________________________________
void MpdKalmanTrack::Clear()
{
  /// Reset track
  delete fParam;
  delete fParamNew;
  delete fParamAtHit;
  delete fCovar;
  delete fWeight;
  delete fWeightAtHit;
  delete fHits;
  fStepMap.clear();
}

//__________________________________________________________________________
TMatrixDSym* MpdKalmanTrack::Weight2Cov()
{
  /// Get covariance matrix from weight matrix
  *fCovar = *fWeight;
  //fCovar->Invert();
  Int_t iok = 0;
  //MpdKalmanFilter::Instance()->MnvertLocal(fCovar->GetMatrixArray(), 5, 5, 5, iok);
  return fCovar;
}

//__________________________________________________________________________
void MpdKalmanTrack::ReSetWeight()
{
  /// Set weight matrix from covariance
  *fWeight = *fCovar;
  Int_t iok = 0;
  //MpdKalmanFilter::Instance()->MnvertLocal(fWeight->GetMatrixArray(), 5, 5, 5, iok);
}

//__________________________________________________________________________
Double_t MpdKalmanTrack::DCA() const
{
  /// Return signed track DCA

  if (fPosNew < 1.e-6) return fPosNew;
  Float_t phiTr = GetParam(2);
  TVector3 vecTr(TMath::Cos(phiTr), TMath::Sin(phiTr), 0.0);
  Float_t phi = GetParam(0) / fPosNew;
  TVector3 vec(TMath::Cos(phi), TMath::Sin(phi), 0.0);
  return fPosNew * TMath::Sign(1.0,vec.Cross(vecTr).Z());
}

//__________________________________________________________________________
void MpdKalmanTrack::PrintNode(Bool_t nNew) const
{
  /// Print node name (for debugging)

  if (nNew) printf("%s\n",fNodeNew.Data()); 
  else printf("%s\n",fNode.Data());
}

///__________________________________________________________________________

ClassImp(MpdKalmanTrack)