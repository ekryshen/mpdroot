#ifndef MPDSTSHIT_H
#define MPDSTSHIT_H

#include "FairHit.h"
#include <TArrayI.h>
#include <TMath.h>
#include <TObject.h>

class MpdStsHit : public FairHit {

public:
   enum Shifts { kSensorS = 1, kLadderS = 6, kLayerS = 11, kSecTypeS = 14 };
   enum Masks { kSensorM = 31, kLadderM = 31, kLayerM = 7, kSecTypeM = 7 };

public:
   /** Default constructor **/
   MpdStsHit();

   /** Constructor with hit parameters (1)**/
   MpdStsHit(Int_t detectorID, TVector3 pos, TVector3 dpos, Double_t signal, Int_t refIndex, Int_t flag); // 14.08

   /** Constructor with hit parameters (2) [not the flag]**/
   MpdStsHit(Int_t detectorID, TVector3 pos, TVector3 dpos, Double_t signal, Int_t refIndex); // 14.08

   /**copy constructor **/
   // MpdStsHit (const MpdStsHit& hit); ///< copy constructor

   /** Destructor **/
   virtual ~MpdStsHit();

   void Print(const Option_t *opt = 0) const;

   /** Accessors **/
   Int_t GetFlag() const { return fFlag; };
   Int_t GetStrip() const { return fStrip; };
   Int_t GetTrackID() const { return fTrackID; };

   //                                               --- Sector layout --- --- Sensor layout ---
   // Int_t    Layer()      const { return GetUniqueID() ? ((fDetectorID >> 11) & 7) : ((fDetectorID >> 11) & 3) + 1; }
   // Int_t    Ladder()     const { return GetUniqueID() ? ((fDetectorID >> 6) & 31) : (fDetectorID >> 6) & 31; }
   // Int_t    Sensor()     const { return GetUniqueID() ? ((fDetectorID >> 1) & 31) : (fDetectorID >> 1) & 31; }
   Int_t SectorType() const { return (fDetectorID >> kSecTypeS) & kSecTypeM; }
   Int_t Layer() const { return GetUniqueID() ? ((fDetectorID >> kLayerS) & kLayerM) : ((fDetectorID >> 11) & 3) + 1; }
   Int_t Ladder() const { return GetUniqueID() ? ((fDetectorID >> kLadderS) & kLadderM) : (fDetectorID >> 6) & 31; }
   Int_t Sensor() const { return GetUniqueID() ? ((fDetectorID >> kSensorS) & kSensorM) : (fDetectorID >> 1) & 31; }
   Int_t Detector() const { return Sensor(); }
   Int_t Side() const { return fDetectorID % 2; }
   Int_t Sector() const { return Sensor(); }
   Int_t Module() const { return Sector(); }

   // Sector layout
   Int_t GetSectorType(Int_t detID) const { return (detID >> kSecTypeS) & kSecTypeM; }
   Int_t GetLayer(Int_t detID) const { return (detID >> kLayerS) & kLayerM; }
   Int_t GetLadder(Int_t detID) const { return (detID >> kLadderS) & kLadderM; }
   Int_t GetSensor(Int_t detID) const { return (detID >> kSensorS) & kSensorM; }
   Int_t GetSide(Int_t detID) const { return detID % 2; }

   Double_t GetSignal() const { return fSignal; }
   Double_t GetLocalX() const { return fLocalX; }

   /** Modifiers **/
   void  SetFlag(Int_t flag) { fFlag = flag; };
   void  SetStrip(Int_t strip) { fStrip = strip; }
   void  SetSignal(Double_t signal) { fSignal = signal; }
   void  SetLocalX(Double_t x) { fLocalX = x; }
   void  SetDetId(Int_t detID) { fDetectorID = detID; }                                /// 111
   Int_t SetDetId(Int_t sectorType, Int_t layer, Int_t ladder, Int_t det, Int_t side); ///< helper for the tracking
   Int_t SetLadder(Int_t detID, Int_t layer) const { return detID & ~(kLadderM << kLadderS) | (layer << kLadderS); }
   Int_t SetSensor(Int_t detID, Int_t sensor) const { return detID & ~(kSensorM << kSensorS) | (sensor << kSensorS); }
   void  SetTrackID(Int_t trackID) { fTrackID = trackID; }

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *hit) const; ///< sort in ascending order of detector layer

protected:
   Int_t    fTrackID;
   Int_t    fFlag;   // Flag for general purposes [TDC, event tagging...]
   Int_t    fStrip;  // Strip number
   Double_t fSignal; // Signal (energy deposit)
   Double_t fLocalX; // Local coordinate (in rotated coordinate system)

   ClassDef(MpdStsHit, 1);
};

#endif
