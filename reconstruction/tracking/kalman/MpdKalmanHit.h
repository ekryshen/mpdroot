#ifndef MPDKALMANHIT_H
#define MPDKALMANHIT_H

/// \ingroup rec
/// \class MpdKalmanHit
/// \brief MPD hit for Kalman tracking
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

#include <TArrayI.h>
#include <TMath.h>
#include <TObject.h>

class MpdKalmanHit : public TObject {

public:
   enum HitType { kFixedP, kFixedR, kFixedZ };
   enum HitFlag {
      kUsed      = 1,
      kMirror    = 2,
      kMlem      = 4,
      kOverflow  = 8,
      kEdge      = 16,
      kVirtual   = 32,
      kSinglePad = 64,
      kSinglePix = 128,
      kMultMax   = 256
   };

public:
   MpdKalmanHit(); ///< Default ctor
   MpdKalmanHit(Int_t detID, Int_t nDim, HitType hitType, Double_t *meas, Double_t *err, Double_t *cosSin,
                Double_t signal, Double_t dist, Int_t index, Double_t edge = 99.); ///< Ctor
   MpdKalmanHit(const MpdKalmanHit &hit);                                          ///< copy constructor
   virtual ~MpdKalmanHit();                                                        ///< Destructor
   Int_t    GetDetectorID() const { return fDetectorID; }                          ///< get detector ID
   Int_t    GetLayer() const { return Int_t(fDetectorID / 1000000); }              ///< get layer No.
   HitType  GetType() const { return fHitType; }                                   ///< Get hit type
   Int_t    GetNofDim() const { return fNofDim; }                                  ///< get number of measur. per point
   Int_t    GetIndex(Int_t indx = 0) const { return fIndex[indx]; }                ///< Get index of detector hit
   TArrayI *Index() { return &fIndex; }                                            ///< Get index array of detector hit
   Int_t    GetFlag() const { return fFlag; }                                      ///< Get flag
   // Bool_t IsMirror() const { return fFlag / 1000; } ///< Get mirror flag
   Bool_t   IsMirror() const { return fFlag & kMirror; }                    ///< Get mirror flag
   Double_t GetLength() const { return fLength; }                           ///< Get track length
   Double_t GetMeas(Int_t indx) const { return fMeas[indx]; }               ///< Get measurement
   Double_t GetErr(Int_t indx) const { return fErr[indx]; }                 ///< Get measurement error
   Double_t GetCosSin(Int_t indx) const { return fCosSin[indx]; }           ///< Cos (Sin)
   Double_t GetPhi() const { return TMath::ATan2(fCosSin[1], fCosSin[0]); } ///< angle
   Double_t GetSignal() const { return fSignal; }                           ///< Signal value
   Double_t GetDist() const { return fDist; }                               ///< Distance to interaction point
   Double_t GetPos() const;                                                 ///< Distance to interaction point
   Double_t GetEdge() const { return fEdge; }                               ///< Distance to sector boundary
   // Double_t GetXY(Int_t indx) const { return fXY[indx]; } ///< get wire X or Y

   void SetDetectorID(Int_t detID) { fDetectorID = detID; } ///< set detector ID
   void SetNofDim(Int_t nDim) { fNofDim = nDim; }           ///< set number of measur. / point
   void SetType(HitType type) { fHitType = type; }          ///< Set hit type
   void SetFlag(Int_t flag) { fFlag = flag; }               ///< Set flag
   // void SetMirror() { fFlag += 1000*TMath::Sign(1,fFlag); } ///< set mirror flag
   void SetMirror() { fFlag |= kMirror; }                            ///< set mirror flag
   void SetLength(Double_t leng) { fLength = leng; }                 ///< set track length
   void SetMeas(Int_t indx, Double_t meas) { fMeas[indx] = meas; }   ///< set measurement
   void SetErr(Int_t indx, Double_t err) { fErr[indx] = err; }       ///< set measurement error
   void SetCosSin(Int_t indx, Double_t cos) { fCosSin[indx] = cos; } ///< set cos(angle) or sin(angle)
   void SetSignal(Double_t signal) { fSignal = signal; }             ///< set signal value
   void SetDist(Double_t dist) { fDist = dist; }                     ///< set distance
   void SetPos(Double_t pos) { fDist = pos; }                        ///< set distance
   void SetIndex(Int_t indx);                                        ///< Add index of detector hit
   void SetEdge(Double_t edge) { fEdge = edge; }                     ///< set distance to sector boundary

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *hit) const; ///< sort in descending order in detector ID
   void   Print(Option_t *opt);              ///< print hit info

private:
   Int_t fDetectorID; ///< detector ID
   Int_t fFlag;       ///< flag
   // Int_t fIndex; ///< ? MC point index ?
   TArrayI fIndex;   ///< hit indices
   HitType fHitType; ///< hit type
   Int_t   fNofDim;  ///< number of measurements per point
   // Bool_t fMirror; ///< flag for mirror hit (left-right ambiguity)
   Double32_t fLength;    ///< track length (temporary entry)
   Double32_t fMeas[2];   ///< measurements
   Double32_t fErr[2];    ///< measurement errors
   Double32_t fCosSin[2]; ///< rotation factors (for stereo measurements)
   Double32_t fSignal;    ///< signal
   Double32_t fDist;      ///< distance to interaction point
   Double32_t fEdge;      ///< distance to sector boundary
   // Double_t fXY[2]; ///< X and Y of some point on the wire (for stereo measurements)

   ClassDef(MpdKalmanHit, 4);
};
#endif
