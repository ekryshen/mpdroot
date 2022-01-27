#ifndef MPDItsHit5spd_H
#define MPDItsHit5spd_H

#include "FairHit.h"
#include <TArrayI.h>
#include <TMath.h>
#include <TObject.h>

class MpdItsHit5spd : public FairHit {

public:
   //  enum Shifts {kSensorS = 1, kLadderS = 6, kLayerS = 11, kSecTypeS = 14};
   //  enum Masks {kSensorM = 31, kLadderM = 31, kLayerM = 7, kSecTypeM = 7};
   enum Shifts { kSensorS = 0, kLadderS = 7, kLayerS = 13, kSecTypeS = 16 };
   enum Masks { kSensorM = 127, kLadderM = 63, kLayerM = 7, kSecTypeM = 7 };

public:
   /** Default constructor **/
   MpdItsHit5spd();

   /** Constructor with hit parameters (1)**/
   MpdItsHit5spd(Int_t detectorID, TVector3 pos, TVector3 dpos, Double_t signal, Int_t refIndex, Int_t flag); // 14.08

   /** Constructor with hit parameters (2) [not the flag]**/
   MpdItsHit5spd(Int_t detectorID, TVector3 pos, TVector3 dpos, Double_t signal, Int_t refIndex); // 14.08

   /**copy constructor **/
   // MpdItsHit_5spd (const MpdItsHit_5spd& hit); ///< copy constructor

   /** Destructor **/
   virtual ~MpdItsHit5spd();

   void Print(const Option_t *opt = 0) const;

   /** Accessors **/
   Int_t GetFlag() const { return fFlag; };
   Int_t GetTrackID() const { return fTrackID; };

   //  Int_t    SectorType() const { return (fDetectorID >> kSecTypeS) & kSecTypeM; }
   Int_t Layer() const { return GetUniqueID() ? ((fDetectorID >> kLayerS) & kLayerM) : ((fDetectorID >> 13) & 3) + 1; }
   Int_t Ladder() const { return GetUniqueID() ? ((fDetectorID >> kLadderS) & kLadderM) : (fDetectorID >> 7) & 31; }
   Int_t Sensor() const { return GetUniqueID() ? ((fDetectorID >> kSensorS) & kSensorM) : (fDetectorID >> 0) & 127; }
   Int_t Detector() const { return Sensor(); }
   Int_t Sector() const { return Sensor(); }
   Int_t Module() const { return Sector(); }

   // Sector layout
   //  Int_t    GetSectorType(Int_t detID) const { return (detID >> kSecTypeS) & kSecTypeM; }
   Int_t GetLayer(Int_t detID) const { return (detID >> kLayerS) & kLayerM; }
   Int_t GetLadder(Int_t detID) const { return (detID >> kLadderS) & kLadderM; }
   Int_t GetSensor(Int_t detID) const { return (detID >> kSensorS) & kSensorM; }

   Double_t GetSignal() const { return fSignal; }
   Double_t GetLocalX() const { return fLocalX; }

   /** Modifiers **/
   void  SetFlag(Int_t flag) { fFlag = flag; };
   void  SetCol(Int_t col) { fCol = col; }
   void  SetRow(Int_t row) { fRow = row; }
   void  SetSignal(Double_t signal) { fSignal = signal; }
   void  SetLocalX(Double_t x) { fLocalX = x; }
   void  SetDetId(Int_t detID) { fDetectorID = detID; }
   Int_t SetDetId(Int_t layer, Int_t ladder, Int_t det); ///< helper for the tracking
   Int_t SetLadder(Int_t detID, Int_t layer) const { return (detID & ~(kLadderM << kLadderS)) | (layer << kLadderS); }
   Int_t SetSensor(Int_t detID, Int_t sensor) const { return (detID & ~(kSensorM << kSensorS)) | (sensor << kSensorS); }
   void  SetTrackID(Int_t trackID) { fTrackID = trackID; }

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *hit) const; ///< sort in ascending order of detector layer

protected:
   Int_t    fTrackID;
   Int_t    fFlag;   // Flag for general purposes [TDC, event tagging...]
   Int_t    fCol;    // Column number
   Int_t    fRow;    // Row number
   Double_t fSignal; // Signal (energy deposit)
   Double_t fLocalX; // Local coordinate (in rotated coordinate system)

   ClassDef(MpdItsHit5spd, 1);
   ;
};

#endif
