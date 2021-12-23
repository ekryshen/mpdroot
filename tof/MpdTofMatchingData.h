//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_MATCHING_DATA_H
#define __MPD_TOF_MATCHING_DATA_H 1

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofMatchingData
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <TObject.h>
#include <TVector3.h>
//------------------------------------------------------------------------------------------------------------------------
class MpdTofHit;
class MpdTofMatchingData : public TObject {
   Int_t fKFTrackIndex = -1; //  KF Track index( branchname: "TpcKalmanTrack", classname: "MpdTpcKalmanTrack")
   Int_t fTofHitIndex  = -1; // if > 0 TofHit index( branchname: "TOFHit", classname: "MpdTofHit")

   Double_t fBeta;            // calc.
   Double_t fMass2;           // calc., [GeV^2]
   Double_t fWeight     = 0.; // matching weight
   Double_t fNormWeight = 0.; // normalized matching weight [0., 1.]
   Double_t fLength;          // Track length [cm], 	{KFTrack refit}
   Int_t    fNmbTrHits = 0;   // number of track hits	{KFTrack refit}
   TVector3 fMomentum;        // Momentum [GeV/c]	{KFTrack copy}
   Double_t fdPhi;            // VR
   Double_t fdZed;            // VR
   Double_t fTime;            // Time since event start [ns]	{TofHit copy}

public:
   // CAUTION: transparent data, MUST be used only at run-time, for debug purpose
   Int_t    fFlag = 0;               //! copy of MpdTofHit::fFlag instance.
   TVector3 fEstPoint;               //! extrapolated point on pad plate(x,y,z), [cm] {KFTrack refit}
   TVector3 fHitPosition;            //! Position of hit [cm] {TofHit copy}
   bool     fIsTrueMatching = false; //!
   bool     fHadTofSignal   = false; //!
   bool     fBestParameter  = false; //!
   double   fBayesPid[4] =
      {}; //! species probability: 0 - proton, 1 - pion, 2 - kaon, 3 - electron(filled by [MpdTofBayesPid::Exec])

   MpdTofMatchingData();
   MpdTofMatchingData(Int_t kfTrackId, Int_t tofHitId, Double_t weight, const MpdTofHit *, Double_t length,
                      Int_t nTrHits, const TVector3 &Momentum, const TVector3 &estPoint);

   bool operator==(const MpdTofMatchingData &rhs)
   {
      return (this->fKFTrackIndex == rhs.fKFTrackIndex && this->fTofHitIndex == rhs.fTofHitIndex);
   }

   void     Print(const char *comment = nullptr, std::ostream &os = std::cout) const;
   Int_t    GetKFTrackIndex(void) const { return fKFTrackIndex; };
   Int_t    GetTofHitIndex(void) const { return fTofHitIndex; };
   Double_t GetTrackLength(void) const { return fLength; };
   Int_t    GetNmbTrHits(void) const { return fNmbTrHits; };
   Double_t GetBeta(void) const { return fBeta; };
   Double_t GetMass2(void) const { return fMass2; };
   TVector3 GetMomentum(void) const { return fMomentum; };
   Double_t GetWeight(void) const { return fWeight; };
   Double_t GetNormWeight(void) const { return fNormWeight; };
   void     SetNormWeight(Double_t v) { fNormWeight = v; };
   Double_t GetTime(void) const { return fTime; };
   Double_t GetdPhi(void) const { return fdPhi; };
   Double_t GetdZed(void) const { return fdZed; };

   // CAUTION: getters for transparent(NOT serialized) data
   Double_t GetDelta(void) const { return (fEstPoint - fHitPosition).Mag(); };

   ClassDef(MpdTofMatchingData, 10)
};
//------------------------------------------------------------------------------------------------------------------------

#endif
