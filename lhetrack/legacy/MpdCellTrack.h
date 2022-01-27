#ifndef MPD_CELLTRACK_H 
#define MPD_CELLTRACK_H 

/// \class MpdCellTrack
/// 
/// Cellular automaton track object for the MPD inner tracking system
/// \author Alexander Zinchenko, Maxim Strelchenko (LHEP, JINR, Dubna)

#include "MpdItsKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
#include "MpdStsHit.h"
#include "FairTask.h"

#include <TArrayI.h> 
#include "TClonesArray.h"
#include "TMath.h" 
#include "TString.h"
#include "TVector3.h" 

class MpdKalmanHit;
class MpdKalmanTrack;
class MpdItsKalmanTrack;

class MpdCellTrack :public TObject 
{
 public:
  enum HitType {kFixedP, kFixedR, kFixedZ};
public:
 /** Constructor **/
 MpdCellTrack(); 
 MpdCellTrack (Int_t detID, Int_t nDim, HitType hitType, TVector3 &meas, Double_t *err, Double_t *cossin, Double_t signal, Double_t dist, Int_t index, Int_t index1, Int_t trackNo); 
 MpdCellTrack (const MpdCellTrack& hit); ///< copy constructor

 /** Destructor **/
 virtual ~MpdCellTrack();

 void SetDetectorID(Int_t detID) { fDetectorID = detID; } ///< set detector ID
 void SetNofDim(Int_t nDim) { fNofDim = nDim; } ///< set number of measur. / point
 void SetType(HitType type) { fHitType = type; } ///< Set hit type
 void SetFlag(Int_t flag) { fFlag = flag; } ///< Set flag
 void SetMirror() { fFlag += 1000*TMath::Sign(1,fFlag); } ///< set mirror flag
 void SetLength(Double_t leng) { fLength = leng; } ///< set track length
 void SetMeas(Int_t indx, Double_t meas) { fMeas[indx] = meas; } ///< set measurement 
 void SetErr(Int_t indx, Double_t err) { fErr[indx] = err; } ///< set measurement error
 void SetCosSin(Int_t indx, Double_t cos) { fCosSin[indx] = cos; } ///< set cos(angle) or sin(angle)
 void SetSignal(Double_t signal) { fSignal = signal; } ///< set signal value
 void SetDist(Double_t dist) { fDist = dist; } ///< set distance
 void SetPos(Double_t pos) { fDist = pos; } ///< set distance
 void SetIndex(Int_t indx); ///< Add index of detector hit
 void SetPrevTrack(Int_t indx) { fTrackNo = indx; } ///< Set track number from previous Layer
 void SetCode(TString code) { fCode = code; } ///< set track code
 void SetCode(Int_t icode) { fCode += icode; fCode += "-"; } ///< set track code

 Int_t GetDetectorID() const { return fDetectorID; } ///< get detector ID
 Int_t GetLayer() const { return Int_t(fDetectorID/1000000); } ///< get layer No.
 HitType GetType() const { return fHitType; } ///< Get hit type
 Int_t GetNofDim() const { return fNofDim; } ///< get number of measur. per point 
 Int_t GetIndex(Int_t indx = 0) const { return fIndex[indx]; } ///< Get index of detector hit
 TArrayI* Index() { return &fIndex; } ///< Get index array of detector hit
 Int_t GetFlag() const { return fFlag; } ///< Get flag
 Bool_t IsMirror() const { return fFlag / 1000; } ///< Get mirror flag
 Double_t GetLength() const { return fLength; } ///< Get track length
 Double_t GetErr(Int_t indx) const { return fErr[indx]; } ///< Get measurement error
 Double_t GetCosSin(Int_t indx) const { return fCosSin[indx]; } ///< Cos (Sin)
 Double_t GetPhi() const { return TMath::ATan2(fCosSin[1],fCosSin[0]); } ///< angle
 Double_t GetSignal() const { return fSignal; } ///< Signal value
 Double_t GetDist() const { return fDist; } ///< Distance to interaction point
 Double_t GetPos() const; ///< Distance to interaction point
 TString GetCode() const { return fCode; } ///< track code

  //Double_t GetXY(Int_t indx) const { return fXY[indx]; } ///< get wire X or Y
 Bool_t IsSortable() const { return kTRUE; }
 Int_t Compare(const TObject* hit) const; ///< sort in descending order in detector ID, "Compare" function for sorting
 void Print(Option_t *opt); ///< print hit info
 TVector3 GetMeas() const { return fMeas; } ///< Get measurement
 Int_t GetPrevTrack() const { return fTrackNo; } ///< Get track number from previous Layer

 private:
 
 Int_t fDetectorID; ///< detector ID
 Int_t fFlag; ///< flag
 TArrayI fIndex; ///< MC point indices
 HitType fHitType; ///< hit type
 Int_t fNofDim; ///< number of measurements per point
 Double32_t fLength; ///< track length (temporary entry)
 Double32_t fErr[2]; ///< measurement errors
 Double32_t fCosSin[4]; ///< rotation factors (for stereo measurements)
 Double32_t fSignal; ///< signal
 Double32_t fDist; ///< distance to interaction point
 TVector3 fMeas;///< measurements
 Int_t fTrackNo; ///< track number from previous layer 
 TString fCode; ///< track code (sequence of 2D-hit indices)

 ClassDef(MpdCellTrack, 1); 
};

#endif
