//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_MATCHING_H
#define __MPD_TOF_MATCHING_H 1

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofMatching
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <assert.h>

#include <TObject.h>
#include <TVector3.h>
#include <TClonesArray.h>
#include <TGeoMatrix.h>
#include <TStopwatch.h>

#include "MpdTpcKalmanTrack.h"
#include "MpdTofHit.h"
#include "MpdTofMatchingData.h"

#include "FairTask.h"
//------------------------------------------------------------------------------------------------------------------------
struct ltPt {
   bool operator()(const MpdKalmanTrack *p1, const MpdKalmanTrack *p2) const //  sorted by Pt descending
   {
      return (p1->Pt() > p2->Pt());
   }
};
class MpdTofPoint;

typedef std::multimap<const MpdKalmanTrack *, Int_t, ltPt> TmPt; // <MpdKalmanTrack*, KfIndex>  sorted by Pt descending
typedef std::multimap<Int_t, MpdTofPoint *>                TmmT2P; // pair< MCtrackID, MpdTofPoint*>

typedef std::multimap<Int_t, std::pair<MpdTofHit *, Int_t>> TmmD2H; // pair< detUID, pair<MpdTofHit*, hitIndex> >
typedef std::map<Int_t, std::pair<MpdTofHit *, Int_t>>      TmS2H;  // pair< stripUID, pair<MpdTofHit*, hitIndex> >
typedef std::multimap<Int_t, std::pair<MpdTofHit *, Int_t>> TmmT2H; // pair< trackID, pair<MpdTofHit*, hitIndex> >
//------------------------------------------------------------------------------------------------------------------------
struct mcInfo //   one-to-many relationship (one kftrack -> many mcpoints)
{
   // track info
   Int_t    trackIndex, pid;
   TVector3 vertMomentum, vertPosition;

   // mc points info
   std::vector<MpdTofPoint *> vPoints;

   bool TofTouch;

   mcInfo() : trackIndex(-1), pid(-1), TofTouch(false) {}

   // return MpdTofPoint number for same tid (kftrack--->mctid<---mcPoint)
   size_t Npoints() const { return vPoints.size(); }
   // return true, if mcPoint exist for this suid
   bool GetPosition(TVector3 &pos, Int_t suid) const;
   // return pair<stripUID, MCtrackId> of closest mcPoint for detector unique ID
   std::pair<Int_t, Int_t> GetClosestPointOnDetector(TVector3 &pos, TVector3 &mom, const TVector3 &from,
                                                     Int_t duid) const;
   // return true, if mcPoint exist for this strip
   bool IsSameStrip(Int_t suid) const;
   // find mcpoint by tid, if suid != -1, ALSO check suid
   const MpdTofPoint *FindPoint(Int_t tid, Int_t suid = -1) const;
};
//------------------------------------------------------------------------------------------------------------------------
class MpdTofMatchingQA;
class LWeightMatrix {
   double *pData    = nullptr; // mData[kf trackId, hitId] = weight
   size_t  fNTracks = 0;       // number of kf tracks
   size_t  fNHits   = 0;       // number of hits
public:
   typedef std::map<long, MpdTofMatchingData> TmCandidate; // pair< hash(kf trackId, hitId), MpdTofMatchingData >
   TmCandidate                                mCandidates;

   // debug info
   size_t                  fNCandTofTracks = 0, fNCandTrueMatch = 0, fNFoundTrueMatch = 0, fNFoundMisMatch = 0;
   MpdTofMatchingQA       *pMatchingQA = nullptr;
   std::map<Int_t, mcInfo> fmT2MC; // pair< kf trackIndex, mcInfo>

   double        SumWeightsForHit(size_t hitId) const;                              // [0,...N-1]
   inline size_t GetTrackID(long uid) const { return ((uid & 0xFFFF0000) >> 16); }; // 4 bytes used from 8
   inline size_t GetHitID(long uid) const { return (uid & 0x0000FFFF); };

   inline long hash(size_t trackId, size_t hitId) const
   {
      long   uid   = (trackId << 16) + hitId;
      size_t track = GetTrackID(uid);
      size_t hit   = GetHitID(uid);
      assert(track == trackId);
      assert(hit == hitId);
      return uid;
   }

public:
   LWeightMatrix(MpdTofMatchingQA *);
   ~LWeightMatrix();

   void AddWeight(size_t trackId, size_t hitId, double value);
   void SetWeight(size_t trackId, size_t hitId, double value);

   inline double GetWeight(size_t trackId, size_t hitId) const
   {
      assert(trackId < fNTracks);
      assert(hitId < fNHits);
      return *(pData + trackId * fNHits + hitId);
   }

   inline size_t GetRank(size_t trackId, size_t &hitIndex) const
   {
      size_t nRank = 0;
      hitIndex     = -1;
      for (size_t hitId = 0; hitId < fNHits; hitId++) {
         if (GetWeight(trackId, hitId) > 0.) {
            nRank++;
            hitIndex = hitId;
         }
      }
      return nRank;
   }

   inline size_t GetRank(size_t trackId, size_t *hitIndex) const
   {
      size_t nRank = 0;
      for (size_t hitId = 0; hitId < fNHits; hitId++) {
         if (GetWeight(trackId, hitId) > 0.) {
            if (nRank < 49) hitIndex[nRank] = hitId;
            nRank++;
         }
      }
      return nRank;
   }

   inline void RemoveHit(size_t hitId)
   {
      for (size_t trackId = 0; trackId < fNTracks; trackId++) SetWeight(trackId, hitId, 0.);
   }

   void                      Normalize(bool copyNormWeight = true);
   void                      Process(TmPt tids, const TClonesArray *aKfTracks);
   void                      CountCandidates(const TmPt &tids);
   void                      CountFoundCandidates();
   void                      Boost1thRankCand(double mult);
   void                      Reset(size_t trackSize, size_t hitSize);
   MpdTofMatchingData       *AddCandidate(const MpdTofMatchingData &data);
   mcInfo                   *saveMCinfo(Int_t trackId, const mcInfo &mc);
   const mcInfo             *getMCinfo(Int_t trackId) const;
   void                      MoveEntries(TClonesArray *aTofMatchings) const;
   const MpdTofMatchingData *FindCandidate(size_t trackId, size_t hitId) const;
   std::pair<size_t, size_t> GetMatchFound() const { return {fNFoundTrueMatch, fNFoundMisMatch}; }
   std::pair<size_t, size_t> GetMatchCand() const { return {fNCandTofTracks, fNCandTrueMatch}; }
   void                      Print(const char *comment = nullptr, const TClonesArray *aKfTracks = nullptr,
                                   const TClonesArray *aTofHits = nullptr, const TmmT2H *map = nullptr, std::ostream &os = std::cout) const;
};
//------------------------------------------------------------------------------------------------------------------------
class TRandom2;
class MpdKalmanFilter;
class MpdTofMatchingQA;
class LStrip;
class LRectangle;

class MpdTofMatching : public FairTask {
   TClonesArray *aMcPoints     = nullptr; //! <--- MC input
   TClonesArray *aMcTracks     = nullptr; //! <--- MC input
   TClonesArray *aTofHits      = nullptr; //! <--- input TOF hits
   TClonesArray *aTPCkfTracks  = nullptr; //! <--- input KF TPC tracks
   TClonesArray *aTofMatchings = nullptr; //! ---> output

   MpdKalmanFilter  *pKF      = nullptr; //!
   LWeightMatrix    *fWeights = nullptr; //!
   TRandom2         *pRandom  = nullptr; //!
   MpdTofMatchingQA *pQA      = nullptr; //!
   TStopwatch       *pTimer   = nullptr; //!

   Bool_t         fDoTest, fIsMcRun = false, fDoMCtest;
   const Double_t fTofBarrelRadius = 152.;                  // [cm]
   Double_t       fThreshZ = 5., fThreshPhi = 5.;           // [cm]
   size_t         fNTofTracksEvent = 0, fNTofTracksRun = 0; // number of Tof touched kf tracks per event & run

   // refit kalman track before propagation
   MpdTpcKalmanTrack RefitTrack(const MpdTpcKalmanTrack *) const;
   // propagate kalman track to cylinder surface.
   TVector3 EstTrackOnR(const MpdTpcKalmanTrack *tr) const;
   // propagate kalman track to plane surface, return error code, 0 - otherwise.
   int EstTrackOnPlate(const MpdTpcKalmanTrack &tr, const TVector3 &point, const TVector3 &perp, TVector3 &pos,
                       Double_t &length, TVector3 &P) const;
   // propagate kalman track to plane surface with smearing track parameters, return true if ok.
   bool EstMarkerTrackOnPlate(size_t mode, const MpdTpcKalmanTrack *tr, const TVector3 &point, const TVector3 &perp,
                              TVector3 &pos, TVector3 &P) const;
   void FillWeights(const TmPt &tids, const TmmD2H &dets, const TmS2H &mHits, const TmmT2P &);

public:
   MpdTofMatching(const char *name = "TOF Matching", Int_t verbose = 1, Bool_t DoTest = false,
                  const char *flnm = "QA.MpdTofMatching.root");
   virtual ~MpdTofMatching();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *option);
   virtual void       Finish();

   static mcInfo GetMcInfo(const TmmT2P &p2ts, const MpdKalmanTrack *pKfTrack, TClonesArray *aMcTracks);

   // mapping pair<MCtrackID, MpdTofPoint*> and sorting by time
   void MapMcT2P(TClonesArray *aMcPoints, TmmT2P &mmMCpoints) const;
   void Print(const MpdTpcKalmanTrack *, const char *comment = nullptr, std::ostream &os = std::cout) const;
   void Report(Int_t verbose, std::ostream &os = std::cout) const;

   ClassDef(MpdTofMatching, 4); // MpdTofMatching
};
//------------------------------------------------------------------------------------------------------------------------
#endif
