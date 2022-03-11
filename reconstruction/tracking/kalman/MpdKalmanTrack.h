#ifndef MPDKALMANTRACK_H
#define MPDKALMANTRACK_H

/// \ingroup rec
/// \class MpdKalmanTrack
/// \brief Kalman track in MPD
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TObjArray.h>
#include <TString.h>
#include <TVector3.h>
#include <map>

class MpdKalmanTrack : public TObject {

public:
   enum TrackDir { kInward, kOutward };
   enum TrackType { kBarrel, kEndcap };
   enum TrackFlag { kOk, kNoOut };

public:
   MpdKalmanTrack();                                       ///< Default ctor
   MpdKalmanTrack(Double_t pos, TVector3 &vertex);         ///< Ctor from position and vertex
   virtual ~MpdKalmanTrack();                              ///< Destructor
   MpdKalmanTrack(const MpdKalmanTrack &track);            ///< copy constructor
   MpdKalmanTrack &operator=(const MpdKalmanTrack &track); ///< assignment operator

   Int_t     GetTrackID() const { return fID; }         ///< get track ID
   TrackDir  GetDirection() const { return fTrackDir; } ///< Get track direction
   TrackType GetType() const { return fTrackType; }     ///< Get track type
   TrackFlag GetFlag() const { return fFlag; }          ///< Get track flag
   Int_t     GetLastLay() const { return fLastLay; }    ///< get the outermost layer
   Int_t     GetNofWrong() const { return fNofWrong; }  ///< number of wrong hits
   TString   GetNode() const { return fNode; }          ///< node path
   TString   GetNodeNew() const { return fNodeNew; }    ///< node path
   Double_t  GetPartID() const { return fPartID; }      ///< get particle ID
   Double_t  GetDedx(Double_t coef = 6.036e-3) const
   {
      return fPartID * coef;
   }                                                      ///< get particle ID (De/dx converted from ADC counts to keV)
   Double_t     GetPos() const { return fPos; }           ///< Get current position
   Double_t     GetPosNew() const { return fPosNew; }     ///< Get propagated position
   Double_t     GetPosAtHit() const { return fPosAtHit; } ///< Get position of outermost hit
   Double_t     GetZ() const { return fTrackType == kEndcap ? fPos : GetParam(1); } ///< Z-position
   Double_t     GetChi2() const { return fChi2; }                                   ///< Get track Chi2
   Double_t     GetChi2Vertex() const { return fChi2Vertex; }                       ///< Get Chi2 to vertex
   Double_t     GetLength() const { return fLength; }                               ///< Get track length
   Double_t     GetLengAtHit() const { return fLengAtHit; }                         ///< Get track length at last hit
   TMatrixD    *GetParam() const { return fParam; }                                 ///< get param. matrix
   TMatrixD    *GetParamNew() const { return fParamNew; }                           ///< get new param. matrix
   TMatrixD    *GetParamAtHit() const { return fParamAtHit; }                       ///< get param. matrix at last hit
   Double_t     GetParam(Int_t i) const { return (*fParam)(i, 0); }                 ///< get param. value
   Double_t     GetParamNew(Int_t i) const { return (*fParamNew)(i, 0); }           ///< get new param. value
   TMatrixDSym *GetCovariance() const { return fCovar; }                            ///< Get covariance matrix
   TMatrixDSym *Weight2Cov();                                   ///< Get covariance matrix from weight matrix
   TMatrixDSym *GetWeight() const { return fWeight; }           ///< Get weight matrix
   TMatrixDSym *GetWeightAtHit() const { return fWeightAtHit; } ///< Get weight matrix at last hit
   Int_t        GetNofHits() const { return fHits ? fHits->GetEntriesFast() : fNhits; } ///< Get number of track hits
   TObjArray   *GetHits() const { return fHits; }                                       ///< Get track hits container
   TVector3    &GetVertex() { return fVertex; }                                         ///< AZ

   void         SetTrackID(Int_t id) { fID = id; }                         ///< set track ID
   void         SetNofHits(Int_t nHits) { fNhits = nHits; }                ///< Set number of track hits
   void         SetDirection(TrackDir dir) { fTrackDir = dir; }            ///< set track direction
   void         SetType(TrackType type) { fTrackType = type; }             ///< set track type
   void         SetFlag(TrackFlag flag) { fFlag = flag; }                  ///< set track flag
   void         SetLastLay(Int_t lay) { fLastLay = lay; }                  ///< set last layer
   void         SetNofWrong(Int_t nWrong) { fNofWrong = nWrong; }          ///< set number of wrong hits
   void         SetNode(TString node) { fNode = node; }                    ///< set node path
   void         SetNodeNew(TString node) { fNodeNew = node; }              ///< set node path
   void         SetPartID(Double_t id) { fPartID = id; }                   ///< set particle ID
   void         SetPos(Double_t pos) { fPos = pos; }                       ///< set position
   void         SetPosNew(Double_t posNew) { fPosNew = posNew; }           ///< set propagated position
   void         SetPosAtHit(Double_t posHit) { fPosAtHit = posHit; }       ///< set position of outermost hit
   void         SetChi2(Double_t chi2) { fChi2 = chi2; }                   ///< set chi2
   void         SetChi2Vertex(Double_t chi2) { fChi2Vertex = chi2; }       ///< set chi2 to vertex
   void         SetLength(Double_t length) { fLength = length; }           ///< set length
   void         SetLengAtHit(Double_t length) { fLengAtHit = length; }     ///< set length at last hit
   void         UpdateLength(Double_t length) { fLength += length; }       ///< update length
   void         SetParam(TMatrixD &param) { *fParam = param; }             ///< set param. matrix
   void         SetParamNew(TMatrixD &param) { *fParamNew = param; }       ///< set new param. matrix
   void         SetParamAtHit(TMatrixD &param) { *fParamAtHit = param; }   ///< set param. matrix at last hit
   void         SetParam(Int_t i, Double_t val) { (*fParam)(i, 0) = val; } ///< set param. value
   void         SetCovariance(TMatrixDSym &cov) { *fCovar = cov; }         ///< set covar. matrix
   void         SetWeight(TMatrixD &weight) { fWeight->SetMatrixArray(weight.GetMatrixArray()); } ///< set weight matrix
   void         SetWeight(TMatrixDSym &weight) { *fWeight = weight; }                             ///< set weight matrix
   void         SetWeightAtHit(TMatrixDSym &weight) { *fWeightAtHit = weight; } ///< set weight matrix at last hit
   void         ReSetWeight();                                                  ///< set weight from covariance
   virtual void StartBack(){};                                                  ///< prepare for back tracing
   virtual void Clear();                                                        ///< reset track (similar to destructor)
   void         PrintNode(Bool_t nNew = kFALSE) const;                          ///< printout node name
   std::map<Int_t, Double_t> &GetSteps() { return fStepMap; }                   ///< steps between points

   // Getters of "physics" track parameters
   Double_t Pt() const { return TMath::Abs(GetParam(4)) > 0.01 ? 1. / TMath::Abs(GetParam(4)) : 100.; }
   Double_t Phi() const { return GetParam(2); }
   Double_t Theta() const { return TMath::PiOver2() - GetParam(3); }
   Double_t Momentum() const { return TMath::Abs(Pt() / TMath::Cos(GetParam(3))); }
   Int_t    Charge() const { return TMath::Nint(TMath::Sign(1., -GetParam(4))); }
   TVector3 Momentum3() const
   {
      return TVector3(Pt() * TMath::Cos(Phi()), Pt() * TMath::Sin(Phi()), Momentum() * TMath::Cos(Theta()));
   }
   Double_t DCA() const; ///< signed DCA
   // Bool_t IsSortable() const { return kTRUE; }
   // Int_t Compare(const TObject* track) const; ///< sort in descending order in Pt

   // private:
protected:
   Int_t                     fID;          ///< trackID
   Int_t                     fNhits;       ///< number of hits (just for convenience of tree browsing)
   TrackDir                  fTrackDir;    ///< track direction (inward/outward)
   TrackType                 fTrackType;   ///< track type (kBarrel/kEndcap)
   Int_t                     fLastLay;     ///< the outermost layer number achieved
   Int_t                     fNofWrong;    ///< number of wrong hits (different from the track ID)
   TString                   fNode;        ///< node path (for local track parameters)
   TString                   fNodeNew;     ///< node path (for local track parameters)
   Double32_t                fPartID;      ///< particle ID
   Double32_t                fPos;         ///< current track radial position
   Double32_t                fPosNew;      ///< propagated track radial position
   Double32_t                fPosAtHit;    ///< radial position of outermost hit
   Double32_t                fChi2;        ///< track Chi2
   Double32_t                fChi2Vertex;  ///< track Chi2-distance to primary vertex
   Double32_t                fLength;      ///< track length
   Double32_t                fLengAtHit;   ///< track length at last hit
   TMatrixD                 *fParam;       ///< pointer to track parameter matrix
   TMatrixD                 *fParamNew;    ///<! pointer to track parameter matrix after propagation
   TMatrixD                 *fParamAtHit;  ///< pointer to track parameter matrix at last hit
   TMatrixDSym              *fCovar;       ///< pointer to covariance matrix
   TMatrixDSym              *fWeight;      ///<! pointer to weight matrix
   TMatrixDSym              *fWeightAtHit; ///< pointer to weight matrix at last hit
   TVector3                  fVertex;      ///< track vertex
   TObjArray                *fHits;        ///<! array of hit pointers assigned to the track
   std::map<Int_t, Double_t> fStepMap;     ///<! steps between points on track
   TrackFlag                 fFlag;        ///< track flag

   ClassDef(MpdKalmanTrack, 5);
};
#endif
