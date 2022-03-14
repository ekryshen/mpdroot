#ifndef MPDECTKALMANTRACK_H
#define MPDECTKALMANTRACK_H

/// \ingroup rec
/// \class MpdEctKalmanTrack
/// \brief Kalman track in MPD straw end-cap detector
///
/// \author Alexander Zinchenko, LHE JINR Dubna

#include "MpdKalmanTrack.h"
#include "MpdKalmanHit.h"
#include <TClonesArray.h>
//#include <TVector3.h>
// class MpdEtofPoint;
class FairHit;
class TpcLheHit;
class MpdTpcKalmanTrack;

class MpdEctKalmanTrack : public MpdKalmanTrack {

public:
   MpdEctKalmanTrack();          ///< Default ctor
   virtual ~MpdEctKalmanTrack(); ///< Destructor
   // MpdEctKalmanTrack(Int_t tpcIndx, const TpcLheKalmanTrack &tpcTrack); ///< Ctor from the TPC track
   MpdEctKalmanTrack(Int_t tpcIndx, const MpdTpcKalmanTrack &tpcTrack); ///< Ctor from TPC track
   MpdEctKalmanTrack(Int_t tofIndx, Int_t tpcIndx, FairHit *tof, TpcLheHit *tpc,
                     TVector3 &vert); ///< Ctor from ETOF and TPC hits
   // MpdEctKalmanTrack(Int_t tofIndx, Int_t ectIndx, MpdEtofPoint *tof, MpdKalmanHitZ *ect, TVector3 &vert); ///< Ctor
   // from ETOF and ECT hits
   MpdEctKalmanTrack(Int_t tofIndx, Int_t ectIndx, FairHit *tof, MpdKalmanHit *ect,
                     TVector3 &vert);                            ///< Ctor from ETOF and ECT hits
   MpdEctKalmanTrack(const MpdEctKalmanTrack &track);            ///< copy constructor
   MpdEctKalmanTrack &operator=(const MpdEctKalmanTrack &track); // assignment operator

   Int_t         GetLastLay() const { return fLastLay; }                    ///< get the outermost layer
   Int_t         GetNofTrHits() const { return fTrHits->GetEntriesFast(); } ///< get number of track hits
   Int_t         GetMirrors() const { return fMirrors; }                    /// get number of mirror hits
   Int_t         GetMisses() const { return fMisses; }                      /// get number of missing tubes
   Int_t         GetTpcIndex() const { return fTpcIndex; }                  /// get index of TPC track or hit
   Int_t         GetTofIndex() const { return fTofIndex; };                 /// get index of TOF hit
   TClonesArray *GetTrHits() const { return fTrHits; }                      ///< get track hits
   Double_t      GetRadNew();                                               ///< get track new radial position
   Double_t     *GetParam1() { return fParam1; }                            ///< get params at first hit
   Int_t         GetFlag() const { return fFlag; }                          /// get flag
   Int_t         IsFromTpc() const { return fFlag & 1; }                    /// is track propagated from TPC?
   Int_t         IsFromEtof() const { return fFlag & 2; }                   /// is track propagated from ETOF?
   Int_t         IsMatchTpc() const { return fFlag & 4; }                   /// is track matched with TPC track?
   Int_t         IsMatchEtof() const { return fFlag & 8; }                  /// is track matched with ETOF hit?

   void SetLastLay() { fLastLay = ((MpdKalmanHit *)GetHits()->Last())->GetLayer(); } ///< set last layer
   void SetLastLay(Int_t lay) { fLastLay = lay; }                                    ///< set last layer
   void SetMirrors(Int_t nMirr) { fMirrors = nMirr; }                                /// set number of mirror hits
   void SetMisses(Int_t nMiss) { fMisses = nMiss; }                                  /// set number of missing tubes
   void SetParam1()
   {
      for (Int_t i = 0; i < 5; ++i) fParam1[i] = GetParamNew(i);
   }                                                ///< set params at first hit
   void SetTofIndex(Int_t tof) { fTofIndex = tof; } /// set TOF hit index
   void SetTpcIndex(Int_t tpc) { fTpcIndex = tpc; } /// set TPC track index
   void SetFlag(Int_t flag) { fFlag = flag; }       /// set flag
   void SetFromTpc() { fFlag |= 1; }                /// track propagated from TPC
   void SetFromEtof() { fFlag |= 2; }               /// track propagated from ETOF
   void SetMatchTpc() { fFlag |= 4; }               /// track matched with TPC track
   void SetMatchEtof() { fFlag |= 8; }              /// track matched with ETOF hit
   void StartBack();                                ///< prepare for back tracing
   void Print(Option_t *opt);                       ///< print track info

   // Bool_t IsSortable() const { return kTRUE; }
   // Int_t Compare(const TObject* track) const; ///< sort in descending order in Pt
   void Reset(); ///< reset track (similar to destructor)

private:
   // void EvalCovar(const MpdKalmanHitR *hitOut, const MpdKalmanHitR *hitIn); ///< evaluate covar. matrix

   Int_t         fLastLay;   ///< the outermost layer number achieved
   Int_t         fTpcIndex;  ///< TPC track or hit index
   Int_t         fTofIndex;  ///< TOF hit index
   Int_t         fFlag;      ///< track flag
   Int_t         fMirrors;   ///< number of mirror hits (due to left-right ambiguity)
   Int_t         fMisses;    ///< number of missing tubes
   TClonesArray *fTrHits;    ///< track hits
   Double_t      fParam1[5]; //!< track param at first hit

   ClassDef(MpdEctKalmanTrack, 2);
};
#endif
