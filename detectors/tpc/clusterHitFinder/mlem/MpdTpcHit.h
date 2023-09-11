#ifndef _MPDTPCHIT_H_
#define _MPDTPCHIT_H_

/// \ingroup tpc
/// \class MpdTpcHit
/// \brief Hit in MPD TPC
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - extension of TpcHit

#include "MpdKalmanHit.h"

#include "AbstractTpcHit.h"
#include <TObject.h>
#include <vector>

class MpdTpc2dCluster;

class MpdTpcHit : public AbstractTpcHit {

public:
   enum k_LinkType { PointIndex = 1, MCTrackIndex = 2 };

public:
   MpdTpcHit()
      : fiPad(-1), fiBin(-1), fLayer(-1), fNdigits(0), fFlag(0), fQ(0), fStep(0), fLength(0), fLocalX(0), fLocalY(0),
        fLocalZ(0)
   {
   }

   MpdTpcHit(Int_t iPad, Int_t iBin)
      : fiPad(iPad), fiBin(iBin), fLayer(-1), fNdigits(0), fFlag(0), fQ(0), fStep(0), fLength(0), fLocalX(0),
        fLocalY(0), fLocalZ(0)
   {
   }

   MpdTpcHit(Int_t detID, TVector3 posHit, TVector3 posHitErr, Int_t pointIndx);

   //     TpcHit(const TpcHit &hit) : FairHit((const FairHit &)hit) {
   //         fiBin = hit.fiBin;
   //         fiPad = hit.fiPad;
   //         fQ = hit.fQ;
   //         fLocalX = hit.fLocalX;
   //         fLocalY = hit.fLocalY;
   //         fLocalZ = hit.fLocalZ;
   //     }

   virtual ~MpdTpcHit() {}

   /* Interface implementation */
   int    GetClusterID() const { return GetRefIndex(); }
   int    GetSector() const { return fSector; }
   double GetPadCoordinate() const { return fLocalX; }
   double GetTimeBinCoordinate() const { return fLocalZ; }
   double GetDriftTime() const { return fDriftTime; }
   float  GetTotalSignal() const { return (float)fQ; }

   /** Accessors **/
   Int_t    GetModular() const { return GetUniqueID(); }
   Int_t    GetPad() const { return fiPad; }
   Int_t    GetBin() const { return fiBin; }
   Int_t    GetLayer() const { return fLayer; }
   Int_t    GetFlag() const { return fFlag; }
   Double_t GetQ() const { return fQ; }
   Double_t GetStep() const { return fStep; }     // step during MC transport
   Double_t GetLength() const { return fLength; } // track length
   // For backward compatibility
   // Int_t GetTrackID() const { return GetLinksWithType(Int_t type).GetLink(Int_t pos); }
   // Int_t GetTrackID() const { return GetLinksWithType(MpdTpcHit::MCTrackIndex).GetLink(0).GetIndex(); }
   Int_t    GetTrackID() const;
   Double_t GetR() const { return fLocalY; }
   Double_t GetRphi() const { return fLocalX; }
   Double_t GetEnergyLoss() const { return fQ; }

   Double_t                            GetLocalX() const { return fLocalX; }
   Double_t                            GetLocalY() const { return fLocalY; }
   Double_t                            GetLocalZ() const { return fLocalZ; }
   void                                LocalPosition(TVector3 &pos) const { pos.SetXYZ(fLocalX, fLocalY, fLocalZ); }
   Double_t                            GetRMS(Int_t ixz = 0) const { return fXZrms[ixz]; }
   Int_t                               GetNdigits() const { return fNdigits; }
   Int_t                               GetNtracks() const { return fIDs.size(); }
   std::vector<Int_t>                 &GetIDs() { return fIDs; }
   std::vector<std::pair<int, float>> &GetTrackIDs() { return vpTrackIDs; }

   /**Get different flags **/
   Int_t IsOverflow() const { return (fFlag & MpdKalmanHit::kOverflow); }
   Int_t IsVirtual() const { return (fFlag & MpdKalmanHit::kVirtual); }
   Int_t IsEdge() const { return (fFlag & MpdKalmanHit::kEdge); }
   Int_t IsMlem() const { return (fFlag & MpdKalmanHit::kMlem); }
   Int_t IsMultMax() const { return (fFlag & MpdKalmanHit::kMultMax); }
   Int_t IsSinglePad() const { return (fFlag & MpdKalmanHit::kSinglePad); }
   Int_t IsSinglePix() const { return (fFlag & MpdKalmanHit::kSinglePix); }

   /** Modifiers **/
   void SetModular(Int_t imod) { SetUniqueID(imod); }
   void SetPad(Int_t ipad) { fiPad = ipad; }
   void SetBin(Int_t ibin) { fiBin = ibin; }
   void SetSector(int sector) { fSector = sector; }
   void SetDriftTime(double driftTime) { fDriftTime = driftTime; }
   void SetLayer(Int_t lay) { fLayer = lay; }
   void SetQ(Double_t q) { fQ = q; }
   void SetStep(Double_t step) { fStep = step; }
   void SetLength(Double_t length) { fLength = length; }
   void SetFlag(Int_t flag) { fFlag = flag; }
   void SetSinglePix() { fFlag |= MpdKalmanHit::kSinglePix; }
   // For backward compatibility
   void SetR(Double_t r) { fLocalY = r; }
   void SetRphi(Double_t rphi) { fLocalX = rphi; }
   void SetEnergyLoss(Double_t edep) { fQ = edep; }

   void SetLocalX(Double_t x) { fLocalX = x; }
   void SetLocalY(Double_t y) { fLocalY = y; }
   void SetLocalZ(Double_t z) { fLocalZ = z; }
   void SetLocalXYZ(Double_t x, Double_t y, Double_t z)
   {
      fLocalX = x;
      fLocalY = y;
      fLocalZ = z;
   }
   void SetLocalPosition(const TVector3 &pos)
   {
      fLocalX = pos.X();
      fLocalY = pos.Y();
      fLocalZ = pos.Z();
   }
   void SetRMS(Double_t rms, Int_t ixz = 0) { fXZrms[ixz] = rms; }
   void SetNdigits(Int_t ndigs) { fNdigits = ndigs; }
   void AddID(Int_t id) { fIDs.push_back(id); }
   void AddTrackID(int id, float charge) { vpTrackIDs.push_back(std::make_pair(id, charge)); }
   void SetFlags(const MpdTpc2dCluster *clus);

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *hit) const; // "Compare" function for sorting

private:
   Int_t                              fiPad;
   Int_t                              fiBin;
   Int_t                              fLayer;
   Int_t                              fNdigits; // number of digits in the hit
   Int_t                              fFlag;
   double                             fDriftTime; // physical drift time of electron to pad area
   int                                fSector;
   std::vector<Int_t>                 fIDs;       // track IDs with the highest charge contribution
   std::vector<std::pair<int, float>> vpTrackIDs; // track IDs with its' charge contribution

   Double32_t fQ;
   Double32_t fStep;
   Double32_t fLength;

   // Sector coordinate system
   Double32_t fLocalX;
   Double32_t fLocalY;
   Double32_t fLocalZ;
   Double32_t fXZrms[2]; // RMS of the hit (group of digits) along pad and drift directions

   ClassDef(MpdTpcHit, 5);
};

#endif // _MPDTPCHIT_H_
