//----------------------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_LMATCHING_FILTER_H
#define __MPD_LMATCHING_FILTER_H

#include <assert.h>
#include <vector>
#include <map>
//#include<debug/map>
#include <set>
#include <iostream>

#include <TClonesArray.h>
#include <TH2D.h>

#include "FairLogger.h"
#include "FairMCTrack.h"

#include "MpdTofMatchingQA.h"
#include "MpdTofMatchingData.h"
#include "MpdTofHit.h"
#include "MpdTofPoint.h"

//#define DoTest
struct CandHistory {
   enum kFlags { kNothing = 0, k1thRank = (1 << 1), kNoCompetitor = (1 << 2), kHaveCompetitor = (1 << 3) };
   Int_t flag;

   std::vector<MpdTofMatchingData *> vCompetitors;
   MpdTofMatchingData *              pCandidate;

   void AddFlag(Int_t f) { flag = flag | f; }
   void ResetFlag(Int_t f)
   {
      Clear();
      flag = f;
   }

   void Clear()
   {
      flag = kNothing;
      vCompetitors.clear();
      pCandidate = nullptr;
   }

   void Print()
   {
      switch (flag) {
      case k1thRank | kNoCompetitor:
         std::cout << "\n  k1thRank & kNoCompetitor pair<" << pCandidate->GetKFTrackIndex() << ", "
                   << pCandidate->GetTofHitIndex() << "> deltaZ=" << pCandidate->GetDelta1()
                   << " deltaPhi=" << pCandidate->GetDelta2();
         break;
      case k1thRank | kHaveCompetitor:
         std::cout << "\n  k1thRank & kHaveCompetitor pair<" << pCandidate->GetKFTrackIndex() << ", "
                   << pCandidate->GetTofHitIndex() << "> deltaZ=" << pCandidate->GetDelta1()
                   << " deltaPhi=" << pCandidate->GetDelta2();
         for (int i = 0; i < vCompetitors.size(); i++) {
            MpdTofMatchingData *ptr = vCompetitors[i];
            std::cout << "\n\t  competitiors pair<" << ptr->GetKFTrackIndex() << ", " << ptr->GetTofHitIndex()
                      << "> deltaZ=" << ptr->GetDelta1() << " deltaPhi=" << ptr->GetDelta2();
         }
         break;
      };
   }
};
//----------------------------------------------------------------------------------------------------------------------------------------
class LMatchingFilter {
public:
   template <typename T>
   using TmmIntT  = std::multimap<Int_t, T>;
   using TmmCand  = TmmIntT<MpdTofMatchingData *>;
   using candIter = TmmCand::iterator;
   using linkIter = TmmIntT<candIter>::iterator;

private:
   TClonesArray *              aMdata; //  Matching data container;
   TmmIntT<MpdTofMatchingData> mcTrueMdata,
      mcMaybeMdata; // pair<kfTrackIndex, MpdTofMatchingData>, filled by mcFindTrueMachings method, if used MC run data.
   TmmIntT<MpdTofMatchingData *> mmCandidate, mmCandidateSnapshot, mmAccepted,
      mmAcceptedSnapshot;                      // pair<KfIndex, MpdTofMatchingData*>
   TmmIntT<candIter> mmLinks, mmLinksSnapshot; // pair<hitIndex, candIter>
   TmmIntT<Int_t>    mmRanks;                  // pair<number of candidate hits for current KFtrack, KfIndex>

   MpdTofMatchingQA *pQA;
   int               fVerbose;
   bool              isOwner;
   bool              isMCdata; // true, if used MC run data. (TofHits and KalmanTracks have links to MC origins.)
   double            fChi2;
   CandHistory       candHist;

   // for debug
   void Commit(void);
   void Status(const char *comment = nullptr, std::ostream &os = std::cout) const;
   void PrintChanges(const TmmIntT<candIter> &src, const TmmIntT<candIter> &checked,
                     std::ostream &os = std::cout) const;
   void PrintChanges(const TmmCand &src, const TmmCand &checked, std::ostream &os = std::cout)
      const; // cycle by src entries, print missed into the checked container entries.
   bool CheckAdequacy(
      const char *comment = nullptr); // check mmCandidate & mmLinks entries accordance; return true if pass test

   Int_t    ProcessSingleTracks(void);
   void     ProcessNRankTrack(Int_t rank);
   void     UpdateRanks();
   void     RemoveLink(Int_t hitIndex, Int_t KfIndex); // remove link to pair<KfIndex, hitIndex>
   void     RemoveCandidate(Int_t hitIndex, Int_t KfIndex);
   candIter FindCandidate(TmmCand *map, Int_t hitIndex,
                          Int_t KfIndex);     // find candidate iter to map pair<KfIndex, hitIndex>
   void     AcceptCandidate(candIter itCand); // move candidate to accepted, clean links and cross candidates
   candIter FindClosest1hRankCandidate(const std::pair<linkIter, linkIter> &range, Double_t &minDelta);
   double   FindMinDeltaForTrack(
        Int_t KfIndex, Int_t disabledHitIndex, candIter &minCand,
        bool &IsWithoutCompetition); // Search through all the  KfIndex candidates except disabledHitIndex

   inline bool IsSameRank(Int_t rank, Int_t KfIndex) // return true, if KfIndex cantidate  have the "rank" rank
   {
      auto range = mmRanks.equal_range(rank);
      for (auto iter = range.first; iter != range.second; ++iter)
         if (iter->second == KfIndex) return true;
      return false;
   }

   inline Int_t GetMinRank() const
   {
      auto iter = mmRanks.begin();
      assert(iter != mmRanks.end());
      return iter->first;
   }
   inline Int_t GetMinRankAfter(Int_t rank) const
   {
      auto iter = mmRanks.upper_bound(rank);
      if (iter != mmRanks.end())
         return iter->first;
      else
         return -1;
   }

public:
   LMatchingFilter(MpdTofMatchingQA *, int verbose = 0);
   ~LMatchingFilter();

   void SetContainer(TClonesArray *array);
   void AddCandidate(const MpdTofMatchingData &data);

   void  Reset();
   Int_t Processing(Int_t nKFTracks, double &chi2);
   Int_t UpdateContainer();

   void Dump(const char *comment = nullptr, std::ostream &os = std::cout) const;
   void PrintRanks(const char *comment = nullptr, std::ostream &os = std::cout) const;

   // method used MC data
   template <typename kfTrackT>
   void mcFindTrueMachings(TClonesArray *aPoints, TClonesArray *aHits, TClonesArray *aMCTracks, TClonesArray *aKfTracks,
                           bool print = false); // work correctly only if used MC run data.

   bool mcIsTrueMatching(const MpdTofMatchingData *);
   bool mcIsMaybeMatching(const MpdTofMatchingData *);
};
//----------------------------------------------------------------------------------------------------------------------------------------

LMatchingFilter::LMatchingFilter(MpdTofMatchingQA *ptr, int verbose)
   : pQA(ptr), fVerbose(verbose), isOwner(true), isMCdata(false)
{
   aMdata = new TClonesArray("MpdTofMatchingData", 1000);
}
//----------------------------------------------------------------------------------------------------------------------------------------
LMatchingFilter::~LMatchingFilter()
{
   if (isOwner) delete aMdata;
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::AddCandidate(const MpdTofMatchingData &data)
{
   MpdTofMatchingData *pData =
      new ((*aMdata)[aMdata->GetEntriesFast()]) MpdTofMatchingData(data); // save data to the TClonesArray container

   candIter it = mmCandidate.insert(
      std::make_pair(pData->GetKFTrackIndex(), pData)); // insert matching candidate into the mmCandidate

   mmLinks.insert(std::make_pair(pData->GetTofHitIndex(), it)); // insert back link to matching candidate
}
//----------------------------------------------------------------------------------------------------------------------------------------
Int_t LMatchingFilter::UpdateContainer()
{
   // select accepted candidates
   std::set<MpdTofMatchingData *> sAccepted;
   for (const auto &cit : mmAccepted) sAccepted.insert(cit.second);

   for (int entry = 0, Nentries = aMdata->GetEntriesFast(); entry < Nentries;
        entry++) // cycle by MpdTofMatchingData at TClonesArray container
   {
      MpdTofMatchingData *pData = (MpdTofMatchingData *)aMdata->At(entry);

      if (sAccepted.find(pData) == sAccepted.end()) // no exist -> not accepted
         aMdata->RemoveAt(entry);
   }

   aMdata->Compress();

   return aMdata->GetEntriesFast();
}
//----------------------------------------------------------------------------------------------------------------------------------------
LMatchingFilter::candIter LMatchingFilter::FindCandidate(TmmCand *map, Int_t hitIndex, Int_t KfIndex)
{
   auto rangeCands = map->equal_range(KfIndex);
   for (candIter iter = rangeCands.first; iter != rangeCands.second; ++iter) // cycle by hits for this KfIndex
   {
      if (hitIndex == iter->second->GetTofHitIndex()) return iter;
   }

   assert(false);
   return map->end();
}
#include "MpdDetectorList.h"
//----------------------------------------------------------------------------------------------------------------------------------------
template <typename kfTrackT>
void LMatchingFilter::mcFindTrueMachings(TClonesArray *aPoints, TClonesArray *aHits, TClonesArray *aMCTracks,
                                         TClonesArray *aKfTracks, bool print)
{
   // Cleanup or initialize container
   mcTrueMdata.clear();
   mcMaybeMdata.clear();

   TmmIntT<kfTrackT *>            mTracks;
   TmmIntT<MpdTofHit *>           mHits;
   std::multimap<kfTrackT *, int> mTracksInv; // pair <Mpd?KalmanTrack*, kftrackIndex>
   std::map<MpdTofHit *, int>     mHitsInv;   // pair <MpdTofHit*, hitIndex>

   // Sorting MpdTpcKalmanTracks by tid
   Int_t nKFTracks = aKfTracks->GetEntriesFast();
   for (size_t KfIndex = 0; KfIndex < nKFTracks; KfIndex++) // cycle by MpdTpcKalmanTrack*
   {
      auto  pKfTrack = (kfTrackT *)aKfTracks->UncheckedAt(KfIndex);
      Int_t tid      = pKfTrack->GetTrackID();

      mTracks.insert(make_pair(tid, pKfTrack));
      mTracksInv.insert(make_pair(pKfTrack, KfIndex));
   }

   // Sorting MpdTofHits by tid
   Int_t nHits = aHits->GetEntriesFast();
   for (size_t hitIndex = 0; hitIndex < nHits; hitIndex++) // cycle by MpdTofHit*
   {
      auto  pHit   = (MpdTofHit *)aHits->UncheckedAt(hitIndex);
      auto  pPoint = (MpdTofPoint *)aPoints->UncheckedAt(pHit->GetRefIndex());
      Int_t tid    = pPoint->GetTrackID();

      mHits.insert(make_pair(tid, pHit));
      mHitsInv.insert(make_pair(pHit, hitIndex));
   }

   // Fill "True" &  "maybe" matching pairs
   Int_t tid = 0; // dummy value, init later
   for (auto itTrack = mTracks.begin(), itTrackEnd = mTracks.end(); itTrack != itTrackEnd;
        itTrack = mTracks.upper_bound(tid)) // cycle by MpdKalmanTrack* with unique tid
   {
      tid = itTrack->first;

      auto mcTrack = (const FairMCTrack *)aMCTracks->UncheckedAt(itTrack->second->GetTrackID());
      if (!mcTrack->GetNPoints(kTOF)) continue; // don't have Tof hit

      auto itHit = mHits.find(tid);
      if (itHit == mHits.end()) continue; // don't have hit for same tid;

      size_t counterT = mTracks.count(tid);
      size_t counterH = mHits.count(tid);

      //  container = mcTrueMdata, if relation one to one (both kftracks to mctrack and hits to mctrack)
      auto container = (counterT == 1 && counterH == 1) ? &mcTrueMdata : &mcMaybeMdata;

      for (size_t t = 0; t < counterT; t++, itTrack++) {
         kfTrackT *track        = itTrack->second;
         Int_t     kfTrackIndex = (mTracksInv.find(track))->second;

         for (size_t h = 0; h < counterH; h++, itHit++) {
            MpdTofHit *hit      = itHit->second;
            Int_t      hitIndex = (mHitsInv.find(hit))->second;

            container->insert(make_pair(kfTrackIndex, MpdTofMatchingData(kfTrackIndex, hitIndex)));
         }
      }
   }

   if (print) {
      std::cout << "\n\t True Matching data: ----------------------------------------------------------------------";
      for (const auto &it : mcTrueMdata) it.second.Print();

      std::cout << "\n\t Maybe Matching data: ----------------------------------------------------------------------";
      for (const auto &it : mcMaybeMdata) it.second.Print();
   }
}
//----------------------------------------------------------------------------------------------------------------------------------------
bool LMatchingFilter::mcIsTrueMatching(const MpdTofMatchingData *pSample)
{
   auto it = mcTrueMdata.find(pSample->GetKFTrackIndex()); // pair<kfindex, MdataT>

   if (it == mcTrueMdata.end()) return false; // no one pair for same kfindex exist;

   return (it->second == (*pSample)); // check by equality operator
   return false;
}
//----------------------------------------------------------------------------------------------------------------------------------------
bool LMatchingFilter::mcIsMaybeMatching(const MpdTofMatchingData *pSample)
{
   Int_t kfindex = pSample->GetKFTrackIndex();
   auto  it      = mcMaybeMdata.find(kfindex); // pair<kfindex, MdataT>

   if (it == mcMaybeMdata.end()) return false; // no one pair for same kfindex exist;

   int counter = mcMaybeMdata.count(kfindex);
   for (int i = 0; i < counter; i++, it++) // cycle by pairs for same kfindex
   {
      if (it->second == (*pSample)) return true;
   }

   return false;
}
//----------------------------------------------------------------------------------------------------------------------------------------
Int_t LMatchingFilter::Processing(Int_t nKFTracks, double &chi2)
{
   fChi2         = 0.;
   Int_t candNmb = mmCandidate.size();

#ifdef DoTest
   std::cout << "\n -I- LMatchingFilter::Processing --------------------------------------------------------- "
                "mmCandidate.size()="
             << mmCandidate.size() << std::flush;
#endif
   Int_t nAccepted, iterNmb = 0;

   if (mmCandidate.empty()) goto finish;

   UpdateRanks();

#ifdef DoTest
   PrintRanks(" FIRST ITER ");
   Commit();
   Dump();
#endif
   ProcessSingleTracks(); // process only single hit (1th rank) candidate tracks
   iterNmb++;

#ifdef DoTest
   Status("AFTER FIRST ProcessSingleTracks");
#endif
   if (mmCandidate.empty()) goto finish;

newIteration:

   iterNmb++;

   UpdateRanks();

#ifdef DoTest
   PrintRanks(" new ITER ");
   int nn = GetMinRankAfter(10);
   if (nn > 0) std::cout << "\n MInRank	= " << nn << " count()= " << mmRanks.count(nn);
   if (mmCandidate.size() <= 20) Dump();
#endif

   nAccepted = ProcessSingleTracks();

   if (nAccepted == 0) {
      ProcessNRankTrack(GetMinRank());
   }

   if (mmCandidate.empty()) goto finish;

   if (iterNmb > 1000) {
      FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN, " <TofMatchingFilter::Processing> Too many tries.");
      goto finish;
   }

   goto newIteration;

finish:
   chi2 = fChi2;

   if (pQA) pQA->FillCandidateNumber(candNmb, iterNmb);

#ifdef DoTest
   std::cout << "\n -I-  [TofMatchingFilter::Processing] FINAL SIZE: mmLinks=" << mmLinks.size()
             << " mmCandidate=" << mmCandidate.size() << " mmRanks=" << mmRanks.size() << "\n";
#endif

   if (fVerbose) std::cout << " -I- [TofMatchingFilter::Processing] Finished with " << iterNmb << " iterations.\n";

   return mmAccepted.size();
}
//----------------------------------------------------------------------------------------------------------------------------------------
int LMatchingFilter::ProcessSingleTracks() //  only single hit (1th rank) candidate tracks processing
{
   auto range = mmRanks.equal_range(1);
   if (std::distance(range.first, range.second) == 0) return 0; // no 1th rank candidates

#ifdef DoTest
   std::cout << "\n LMatchingFilter::ProcessSingleTracks_size=" << std::distance(range.first, range.second);
#endif

   Int_t MatchingOK = 0;
   for (auto iter = range.first; iter != range.second; ++iter) // cycle by ONLY 1th rank candidates
   {
      Int_t    KfIndex = iter->second;
      candIter itCand  = mmCandidate.find(KfIndex);
#ifdef DoTest
      std::cout << "\nProcessSingleTracks_KfIndex=" << KfIndex;
#endif

      if (itCand == mmCandidate.end()) continue; // already removed by competitor track candidate
#ifdef DoTest
      std::cout << " EXIST";
#endif
      Int_t hitIndex = itCand->second->GetTofHitIndex();
      Int_t counter  = mmLinks.count(hitIndex); // number of competitor tracks for this hit

      if (counter == 1) // don't have  competitor tracks
      {

         ///			candHist.ResetFlag(CandHistory< MdataT>::k1thRank | CandHistory< MdataT>::kNoCompetitor);
         ///			candHist.pCandidate = itCand->second;

         AcceptCandidate(itCand);
         MatchingOK++;
      } else // accept best(closest) 1th rank candidate
      {
         Double_t delta;
         candIter it = FindClosest1hRankCandidate(mmLinks.equal_range(hitIndex), delta);
         ///			candHist.pCandidate = it->second;

         AcceptCandidate(it);
         MatchingOK++;
      }
   }

#ifdef DoTest
   std::cout << "\n LMatchingFilter::ProcessSingleTracks 1th rank cands=" << std::distance(range.first, range.second)
             << " MatchingOK=" << MatchingOK;
#endif

   return MatchingOK;
}
//----------------------------------------------------------------------------------------------------------------------------------------
double LMatchingFilter::FindMinDeltaForTrack(Int_t KfIndex, Int_t disabledHitIndex, candIter &minCand,
                                             bool &IsWithoutCompetition)
{
   Double_t delta, minDelta = 1.e+10; // big value
   auto     rangeCands  = mmCandidate.equal_range(KfIndex);
   minCand              = mmCandidate.end();
   IsWithoutCompetition = true;

   assert(std::distance(rangeCands.first, rangeCands.second) != 0);

   for (candIter iter = rangeCands.first; iter != rangeCands.second; ++iter) // cycle by candidates for this KfIndex
   {
      int hitindex  = iter->second->GetTofHitIndex();
      int nmbTrHits = iter->second->GetNmbTrHits();

      if (disabledHitIndex != hitindex) // skip <KfIndex, disabledHitIndex> candidate
      {
         delta = iter->second->GetDelta();

         if (delta < minDelta) // new local minimum
         {
            IsWithoutCompetition = (1 == mmLinks.count(hitindex)) ? true : false;

            minCand  = iter;
            minDelta = delta;
         }
      }
   }

   assert(minCand != mmCandidate.end()); // don't have local minimum ????
   return minDelta;
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::ProcessNRankTrack(Int_t rank) // process  Nth rank  ONE track
{
   auto rangeRanks = mmRanks.equal_range(rank);
   if (std::distance(rangeRanks.first, rangeRanks.second) == 0) return; // no Nth rank candidates

   // select random(first) KfIndex from Nrank list
   Int_t KfIndex = rangeRanks.first->second; // track index for processing

#ifdef DoTest
   std::cout << "\n  ---------------ProcessNRankTrack---- rank=" << rank
             << "  mmRanks.count(rank)=" << std::distance(rangeRanks.first, rangeRanks.second) << " KfIndex=" << KfIndex
             << "  mmCandidate.count(KfIndex)=" << mmCandidate.count(KfIndex);
#endif
   bool     IsWithoutCompetition = true;
   candIter cand, bestCand;

   Double_t minDelta         = FindMinDeltaForTrack(KfIndex, -100, bestCand, IsWithoutCompetition);
   int      maxKfTrackWeight = bestCand->second->GetNmbTrHits();

   if (IsWithoutCompetition) {
      assert(cand != mmCandidate.end());
      AcceptCandidate(bestCand);
   } else {
      Int_t hitIndex   = bestCand->second->GetTofHitIndex();
      auto  rangeLinks = mmLinks.equal_range(hitIndex);

      for (linkIter it = rangeLinks.first; it != rangeLinks.second; ++it) // cycle by tracks for this hitIndex
      {
         if (it->second == bestCand) continue; // skip himself

         Int_t kfindex = it->second->first;
         if (IsSameRank(rank, kfindex)) // check only same rank candidates
         {
            Double_t delta         = FindMinDeltaForTrack(kfindex, -100, cand, IsWithoutCompetition);
            int      kftrackweight = cand->second->GetNmbTrHits();

            if (kftrackweight > maxKfTrackWeight) // new track have bigger weight -> set as best candidate
            {
               bestCand = cand;
               minDelta = delta;
            } else if (kftrackweight == maxKfTrackWeight) {
               if (delta < minDelta) // new track have same weight but smaller delta -> set as best candidate
               {
                  bestCand = cand;
                  minDelta = delta;
               }
            }
         }
      } // cycle by tracks for this hitIndex

      assert(cand != mmCandidate.end());
      AcceptCandidate(bestCand);
   }
}
//----------------------------------------------------------------------------------------------------------------------------------------
LMatchingFilter::candIter LMatchingFilter::FindClosest1hRankCandidate(const std::pair<linkIter, linkIter> &range,
                                                                      Double_t &                           minDelta)
{
   ///////////	candHist.ResetFlag(CandHistory< MdataT>::k1thRank | CandHistory< MdataT>::kHaveCompetitor);

   candIter retvalue = mmCandidate.end();
   minDelta          = 1.e+10; // big value

   for (linkIter it = range.first; it != range.second; it++) // cycle by competitor candidates
   {
      Double_t delta = it->second->second->GetDelta();

      /////////////		candHist.vCompetitors.push_back(it->second->second);

      if (delta < minDelta) // best (closest) candidate
      {
         retvalue = it->second;
         minDelta = delta;
      }
   }

   assert(retvalue != mmCandidate.end());
   return retvalue;
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::AcceptCandidate(candIter itCand)
{
   bool isTrueMatching      = mcIsTrueMatching(itCand->second);
   bool isMaybeTrueMatching = mcIsMaybeMatching(itCand->second);

#ifdef DoTest
   if (!(isTrueMatching || isMaybeTrueMatching)) // mismatching accepted !!!
   {
      ///		candHist.Print();
      // itCand->second->Print();
   }
   Commit();
#endif
   // Accept candidate
   Int_t acceptedKfIndex = itCand->first, acceptedHitIndex = itCand->second->GetTofHitIndex();

   mmAccepted.insert(make_pair(acceptedKfIndex, itCand->second)); // pair< KfIndex, MdataT*>

   double delta = itCand->second->GetDelta();
   fChi2 += delta * delta;

   auto rangeCand = mmCandidate.equal_range(acceptedKfIndex);
   for (candIter iter = rangeCand.first, iterEnd = rangeCand.second; iter != iterEnd;
        iter++) // all candidates for acceptedKfIndex
   {
      int hitindex = iter->second->GetTofHitIndex();

      // remove all links to hits for this KfIndex	(track erase)
      RemoveLink(hitindex, acceptedKfIndex);
   }

   auto rangeLinks = mmLinks.equal_range(acceptedHitIndex);
   for (linkIter iter = rangeLinks.first, iterEnd = rangeLinks.second; iter != iterEnd;
        iter++) // all candidates for acceptedHitIndex
   {
      int kfindex = iter->second->first;

      // remove all candidate  tracks for this hitIndex	(hit erase)
      RemoveCandidate(acceptedHitIndex, kfindex);
   }

   // remove all candidates for  this KfIndex
   rangeCand = mmCandidate.equal_range(acceptedKfIndex);
   mmCandidate.erase(rangeCand.first, rangeCand.second);

   // remove all links for this hitIndex hit	(hit erase)
   rangeLinks = mmLinks.equal_range(acceptedHitIndex);
   mmLinks.erase(rangeLinks.first, rangeLinks.second);

#ifdef DoTest
   bool pass = CheckAdequacy("AFTER ACCEPT");
   if (!pass) {
      Status(" change AFTER ACCEPT");
      std::cout << std::flush;
   }
   assert(pass);
#else
   assert(CheckAdequacy("AFTER ACCEPT"));
#endif
}
//----------------------------------------------------------------------------------------------------------------------------------------
bool LMatchingFilter::CheckAdequacy(const char *comment)
{
   TmmIntT<candIter> mmLinksTmp;
   //	for(const auto& entry : mmCandidate)   mmLinksTmp.insert(make_pair(entry.second->GetTofHitIndex(),
   //mmCandidate.find(entry.first)));
   for (candIter it = mmCandidate.begin(), itEnd = mmCandidate.end(); it != itEnd; it++)
      mmLinksTmp.insert(make_pair(it->second->GetTofHitIndex(), it));

   bool equalFront =
      true; // (mmLinksTmp == mmLinks); // iterators point to another mmaps, so operator == return always true
   if (mmLinksTmp.size() == mmLinks.size()) {
      for (const auto &iter : mmLinksTmp) {
         bool founded    = false;
         auto rangeLinks = mmLinks.equal_range(iter.first);
         for (auto it = rangeLinks.first; it != rangeLinks.second; ++it) // cycle by tracks for this hitIndex
         {
            if (it->second->second == iter.second->second) // MdataT*(mmLinksTmp) == MdataT* (mmLinks)
            {
               founded = true;
               break;
            }
         }

         if (false == founded) {
            equalFront = false;
            break;
         }
      }
   } else
      equalFront = false;

   TmmCand mmCandidateTmp;
   for (const auto &entry : mmLinks) mmCandidateTmp.insert(make_pair(entry.second->first, entry.second->second));

   bool equalBack = true; // (mmCandidateTmp == mmCandidate); // because of a different order filling the entries may be
                          // mismatch, so operator == may be return wrong false

   if (mmCandidateTmp.size() == mmCandidate.size()) {
      for (const auto &iter : mmCandidateTmp) {
         bool founded = false;
         auto range   = mmCandidate.equal_range(iter.first);

         for (auto it = range.first; it != range.second; ++it) // cycle by candidates for same KfIndex
         {
            if (it->second == iter.second) // MdataT*(mmCandidateTmp) == MdataT* (mmCandidate)
            {
               founded = true;
               break;
            }
         }

         if (false == founded) {
            equalFront = false;
            break;
         }
      }
   } else
      equalBack = false;

   if (!equalFront || !equalBack) // print report, if it's problems
   {
      std::cout << "\n CheckAdequacy: ";
      if (comment) std::cout << comment;
      std::cout << " mmCandidate.size=" << mmCandidate.size() << " mmLinks.size=" << mmLinks.size()
                << " (mmLinks ---> mmCandidate) = " << equalBack << ", (mmCandidate ---> mmLinks) = " << equalFront;
   }

   return (equalFront && equalBack);
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::RemoveLink(Int_t hitIndex, Int_t KfIndex)
{
   auto range = mmLinks.equal_range(hitIndex); // links for same  hitIndex

   for (linkIter it = range.first, itEnd = range.second; it != itEnd;) // iterating multimap && erasing
   {
      if (it->second->first == KfIndex)
         mmLinks.erase(it++);
      else
         ++it;
   }
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::RemoveCandidate(Int_t hitIndex, Int_t KfIndex)
{
   auto range = mmCandidate.equal_range(KfIndex);

   for (candIter it = range.first, itEnd = range.second; it != itEnd;) // iterating multimap && erasing
   {
      if (it->second->GetTofHitIndex() == hitIndex)
         mmCandidate.erase(it++);
      else
         ++it;
   }
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::UpdateRanks()
{
   mmRanks.clear();

   for (auto it = mmCandidate.begin(), end = mmCandidate.end(); it != end;
        it = mmCandidate.upper_bound(it->first)) // cycle by unique KfIndex
   {
      mmRanks.insert(make_pair(mmCandidate.count(it->first), it->first)); // pair <count(KfIndex), KfIndex>
   }
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::Reset()
{
   mmCandidate.clear();
   mmAccepted.clear();
   mmLinks.clear();

   candHist.Clear();
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::SetContainer(TClonesArray *array)
{
   aMdata  = array;
   isOwner = false;
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::Commit()
{
   mmCandidateSnapshot = mmCandidate;
   mmAcceptedSnapshot  = mmAccepted;

   // mmLinksSnapshot = mmLinks; ERROR: DON'T USE THIS:  mmlinks used iterators to mmCandidate, so it may be singular
   // iterator!!!!
   mmLinksSnapshot.clear();
   candIter iter;

   for (const auto &it : mmLinks) {
      Int_t hitIndex = it.first;
      iter           = FindCandidate(
         &mmCandidateSnapshot, hitIndex,
         it.second->second->GetKFTrackIndex()); //  mmLinksSnapshot used iterators to mmCandidateSnapshot NOW

      mmLinksSnapshot.insert(make_pair(hitIndex, iter));
   }
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::Status(const char *comment, std::ostream &os) const
{
   os << "\n -------------------------------------------------STATUS-----  ";
   if (comment != nullptr) os << comment;

   os << "\n mmCandidate added:";
   PrintChanges(mmCandidate, mmCandidateSnapshot, os);

   os << "\n mmCandidate removed:";
   PrintChanges(mmCandidateSnapshot, mmCandidate, os);

   os << "\n mmAccepted added:";
   PrintChanges(mmAccepted, mmAcceptedSnapshot, os);

   os << "\n mmAccepted removed:";
   PrintChanges(mmAcceptedSnapshot, mmAccepted, os);

   os << "\n mmLinks added:";
   PrintChanges(mmLinks, mmLinksSnapshot, os);

   os << "\n mmLinks removed:";

   PrintChanges(mmLinksSnapshot, mmLinks, os);
   os << "\n ------------------------------------------------------------" << std::flush;
}
//----------------------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::PrintChanges(const TmmCand &src, const TmmCand &checked, std::ostream &os) const
{
   if (src == checked) {
      os << "\t   nothing.";
      return;
   }

   Int_t KfIndex, hitIndex;
   bool  exist;
   for (const auto &it : src) {
      KfIndex    = it.first;
      auto range = checked.equal_range(KfIndex); // cands for same  KfIndex

      exist    = false;
      hitIndex = it.second->GetTofHitIndex();

      for (auto iter = range.first; iter != range.second;
           ++iter) // cycle by same KfIndex entries from checked container
      {
         if (hitIndex == iter->second->GetTofHitIndex()) // same hitIndex found
         {
            exist = true;
            break;
         }
      }

      if (!exist) os << "\n    KfIndex=" << KfIndex << ", hitIndex=" << hitIndex;
   }
}
//------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::PrintChanges(const TmmIntT<candIter> &src, const TmmIntT<candIter> &checked,
                                   std::ostream &os) const
{
   if (src == checked) {
      os << "\t   nothing.";
      return;
   }

   Int_t KfIndex, hitIndex;
   bool  exist;
   for (const auto &it : src) {
      hitIndex   = it.first;
      auto range = checked.equal_range(hitIndex); // links for same  hitIndex

      exist   = false;
      KfIndex = it.second->second->GetKFTrackIndex();

      for (auto iter = range.first; iter != range.second;
           ++iter) // cycle by same hitIndex entries from checked container
      {
         if (KfIndex == iter->second->second->GetKFTrackIndex()) // same KfIndex found
         {
            exist = true;
            break;
         }
      }

      if (!exist) os << "\n    KfIndex=" << KfIndex << ", hitIndex=" << hitIndex;
   }
}
//------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::Dump(const char *comment, std::ostream &os) const
{
   os << "\n ---------------------- DUMP --------------- ";
   if (comment) os << comment;
   os << " ------------   entries: " << mmCandidate.size();

   Int_t KfIndex, hitIndex;
   bool  exist;
   int   n = 0;
   for (auto it = mmCandidate.begin(); it != mmCandidate.end();) {
      KfIndex = it->first;
      os << "\n " << ++n << ") KfIndex=" << KfIndex << " --- ";
      int counter = mmCandidate.count(KfIndex);

      for (int hit = 0; hit < counter; hit++) {
         hitIndex = it->second->GetTofHitIndex();
         os << " hitIndex=" << hitIndex;

         auto iter = mmLinks.find(hitIndex);
         if (iter != mmLinks.end()) {
            int counter2 = mmLinks.count(hitIndex);
            if (counter2 > 1) os << " (";
            for (int track = 0; track < counter2; track++) {
               if (counter2 > 1) os << " " << iter->second->second->GetKFTrackIndex() << " ";
               ++iter;
            }
            if (counter2 > 1) os << ") ";
         }
         ++it;
      }
   }
   os << "\n ---------------------- DUMP --------------------------------------" << std::flush;
}
//------------------------------------------------------------------------------------------------------------------------
void LMatchingFilter::PrintRanks(const char *comment, std::ostream &os) const
{
   os << std::endl;
   if (comment) os << comment;
   os << "\n mmRanks.count(n) 1= " << mmRanks.count(1) << ", 2= " << mmRanks.count(2) << ", 3= " << mmRanks.count(3)
      << ", 4= " << mmRanks.count(4) << ", 5= " << mmRanks.count(5) << ", 6= " << mmRanks.count(6)
      << ", 7= " << mmRanks.count(7) << ", 8= " << mmRanks.count(8) << ", 9= " << mmRanks.count(9)
      << ", 10= " << mmRanks.count(10) << "  cand. total: " << mmCandidate.size();
}
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_LMATCHING_FILTER_H
