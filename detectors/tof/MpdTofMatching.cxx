//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofMatching
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <limits>

#include <TMath.h>
#include <TFile.h>
#include <TRandom2.h>
#include <TRotation.h>
#include <TStorage.h>
#include <TEfficiency.h>

#include "FairRootManager.h"
#include "MpdMCTrack.h"
#include "FairRunAna.h"
#include "FairLogger.h"

#include "MpdTofPoint.h"
#include "MpdTof.h"
#include "MpdKalmanFilter.h"
#include "MpdTofMatchingQA.h"
#include "MpdTofGeoUtils.h"

#include "MpdTofMatching.h"

#ifdef _OPENMP
#include "omp.h"
#include <sys/time.h>
#endif

using namespace std;
//------------------------------------------------------------------------------------------------------------------------
bool mcInfo::GetPosition(TVector3 &pos, Int_t suid) const
{
   for (const auto &cit : vPoints) // cycle by vector<MpdTofPoint*> (selected kftrack -> mctrack -> mcpoints)
   {
      if (suid == cit->GetDetectorID()) {
         cit->Position(pos);
         return true;
      }
   }

   return false;
}
//------------------------------------------------------------------------------------------------------------------------
pair<Int_t, Int_t> mcInfo::GetClosestPointOnDetector(TVector3 &pos, TVector3 &mom, const TVector3 &from,
                                                     Int_t duid) const
{
   pair<Int_t, Int_t> found = {-1, -1}; // {suid, mctid}
   if (vPoints.empty()) return found;

   double   dev = 1.e+10;
   TVector3 point, P;

   for (const auto &cit : vPoints) // cycle by vector<MpdTofPoint*> (selected kftrack -> mctrack -> mcpoints)
   {
      if (MpdTofPoint::IsSameDetector(duid, cit->GetSuid())) {
         cit->Position(point);
         double newdev = (point - from).Mag();
         if (newdev < dev) {
            dev = newdev;
            pos = point;
            cit->Momentum(mom);
            found.first  = cit->GetSuid();
            found.second = cit->GetTrackID();
         }
      }
   }

   return found;
}
//------------------------------------------------------------------------------------------------------------------------
bool mcInfo::IsSameStrip(Int_t suid) const
{
   if (vPoints.empty()) return false;

   for (const auto &cit : vPoints) // cycle by vector<MpdTofPoint*> (selected kftrack -> mctrack -> mcpoints)
   {
      if (MpdTofPoint::IsSameStrip(cit->GetDetectorID(), suid)) return true;
   }

   return false;
}
//------------------------------------------------------------------------------------------------------------------------
const MpdTofPoint *mcInfo::FindPoint(Int_t tid, Int_t suid) const
{
   if (vPoints.empty()) return nullptr;

   for (const auto &cit : vPoints) // cycle by vector<MpdTofPoint*> (selected kftrack -> mctrack -> mcpoints)
   {
      if (cit->GetTrackID() == tid) {
         if (-1 == suid)
            return cit;
         else if (cit->GetDetectorID() == suid)
            return cit;
      }
   }

   return nullptr;
}
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
LWeightMatrix::LWeightMatrix(MpdTofMatchingQA *ptr) : pData(nullptr), pMatchingQA(ptr) {}
//------------------------------------------------------------------------------------------------------------------------
LWeightMatrix::~LWeightMatrix()
{
   delete[] pData;
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::SetWeight(size_t trackId, size_t hitId, double value) // [0,...N-1]
{
   assert(trackId < fNTracks);
   assert(hitId < fNHits);
   *(pData + trackId * fNHits + hitId) = value;
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::AddWeight(size_t trackId, size_t hitId, double value)
{
   assert(trackId < fNTracks);
   assert(hitId < fNHits);
   *(pData + trackId * fNHits + hitId) += value;
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::Normalize(bool copyNormWeight)
{
   for (size_t hitId = 0; hitId < fNHits; hitId++) {
      double sum = SumWeightsForHit(hitId);

      if (sum > 0.) {
         for (size_t trackId = 0; trackId < fNTracks; trackId++) {
            double norm_weight = GetWeight(trackId, hitId) / sum;
            SetWeight(trackId, hitId, norm_weight);

            if (copyNormWeight && norm_weight > 0.) {
               auto it = mCandidates.find(hash(trackId, hitId));
               assert(it != mCandidates.end());
               it->second.SetNormWeight(norm_weight);
            }
         }
      }
   }
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::Boost1thRankCand(double mult)
{
   if (mult == 1.) return; // do nothing

   for (size_t trackId = 0; trackId < fNTracks; trackId++) {
      size_t nHits = 0, singleHitIndex;
      double singleHitWeight;

      for (size_t hitId = 0; hitId < fNHits; hitId++) {
         double weight = GetWeight(trackId, hitId);
         if (weight > 0.) {
            nHits++;
            singleHitIndex  = hitId;
            singleHitWeight = weight;
         }
         if (nHits > 1) break;
      }

      // boost 1th rank candidate.
      if (nHits == 1) SetWeight(trackId, singleHitIndex, singleHitWeight * mult);
   }
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::Reset(size_t tracks, size_t hits)
{
   if (tracks * hits > fNTracks * fNHits) // buffer size increased --> MUST recreate buffer
   {
      delete[] pData;
      pData = new double[tracks * hits];
   }

   fNTracks = tracks;
   fNHits   = hits;

   // cleanup
   memset(pData, 0, sizeof(double) * tracks * hits);
   mCandidates.clear();
   fmT2MC.clear();
}
//------------------------------------------------------------------------------------------------------------------------
double LWeightMatrix::SumWeightsForHit(size_t hitId) const
{
   assert(hitId < fNHits);
   double sum = 0.;

   for (size_t trackId = 0; trackId < fNTracks; trackId++) sum += GetWeight(trackId, hitId);

   return sum;
}
//------------------------------------------------------------------------------------------------------------------------
mcInfo *LWeightMatrix::saveMCinfo(Int_t trackId, const mcInfo &mc)
{
   auto ret = fmT2MC.insert(make_pair(trackId, mc));

   return &(ret.first->second);
}
//------------------------------------------------------------------------------------------------------------------------
const mcInfo *LWeightMatrix::getMCinfo(Int_t trackId) const
{
   auto it = fmT2MC.find(trackId);
   if (it == fmT2MC.end()) return nullptr;

   return &(it->second);
}
//----------------------------------------------------------------------------------------------------------------------------------------
MpdTofMatchingData *LWeightMatrix::AddCandidate(const MpdTofMatchingData &data)
{
   const long key = hash(data.GetKFTrackIndex(), data.GetTofHitIndex());
   auto       ret = mCandidates.insert(make_pair(key, data));
   return &(ret.first->second);
}
//----------------------------------------------------------------------------------------------------------------------------------------
const MpdTofMatchingData *LWeightMatrix::FindCandidate(size_t trackId, size_t hitId) const
{
   auto it = mCandidates.find(hash(trackId, hitId));

   if (it != mCandidates.cend()) return &(it->second);

   return nullptr;
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::Print(const char *comment, const TClonesArray *aKfTracks, const TClonesArray *aHits,
                          const TmmT2H *mmT2H, std::ostream &os) const
{
   auto prec = os.precision(5);
   os << "\n [LWeightMatrix::Print]-------------------------------------------------------------------->>> ";
   if (comment != nullptr) os << comment;

   for (size_t trackId = 0; trackId < fNTracks; trackId++) // cycle by kf track index
   {
      auto pKfTrack = (MpdTpcKalmanTrack *)aKfTracks->UncheckedAt(trackId);
      auto mctid    = pKfTrack->GetTrackID();
      os << "\nKfTid=" << trackId << " tid=" << mctid;

      // find max weight hit index
      size_t iMaxWeight = -1;
      double maxWeight  = 0.;
      for (size_t hitId = 0; hitId < fNHits; hitId++) // cycle by hit index
      {
         double weight = GetWeight(trackId, hitId);
         if (weight > maxWeight) {
            maxWeight  = weight;
            iMaxWeight = hitId;
         }
      }

      auto info    = fmT2MC.find(trackId); // pair<kf trackIndex, mcInfo>
      int  nPoints = 0;

      if (info != fmT2MC.end()) {
         nPoints = info->second.Npoints();
         os << "  (" << nPoints << ")\t";
      }

      if (mmT2H && info != fmT2MC.end()) {
         auto range = mmT2H->equal_range(info->second.trackIndex); // pair< trackID, pair<MpdTofHit*, hitIndex> >
         for (auto i = range.first; i != range.second; ++i) os << " hid=" << i->second.second;
      }

      os << "\t\t";

      if (aHits)
         for (size_t hitId = 0; hitId < fNHits; hitId++) // cycle by hit index
         {
            double weight = GetWeight(trackId, hitId);

            auto  pHit = (MpdTofHit *)aHits->At(hitId);
            Int_t suid = pHit->GetDetectorID();

            auto cand = FindCandidate(trackId, hitId); // MpdTofMatchingData

            Int_t sector, detector, gap, strip;
            MpdTofPoint::ParseSuid(suid, sector, detector, gap, strip);

            bool IsTrueMatching = pHit->CheckTid(mctid);

            if (weight > 0. || IsTrueMatching) {
               os << " hitId=" << hitId << ",w=[" << weight << "]";
               os << "{" << sector << "," << detector << "," << gap << "," << strip << "}";
               if (cand) os << "<" << cand->GetDelta() << ">";
            }

            if (IsTrueMatching) os << " *T*";

            if (iMaxWeight == hitId) // hit with max weight
            {
               if (IsTrueMatching)
                  os << "<<< ";
               else if (nPoints != 0)
                  os << " @@MISMATCH@@";
            }
         }
   }

   os << "\n [LWeightMatrix::Print]--------------------------------------------------------------------<<< ";
   os.precision(prec);
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::CountFoundCandidates()
{
   if (pMatchingQA) {
      size_t nFoundTrueMatch = 0, nFoundMisMatch = 0;

      for (auto cand = mCandidates.cbegin(), candEnd = mCandidates.cend(); cand != candEnd; cand++) {
         pMatchingQA->FillParameter(cand->second.fIsTrueMatching, cand->second.GetWeight(),
                                    cand->second.GetNormWeight());

         if (cand->second.fIsTrueMatching)
            nFoundTrueMatch++;
         else
            nFoundMisMatch++;
      }

      fNFoundTrueMatch += nFoundTrueMatch;
      fNFoundMisMatch += nFoundMisMatch;
   }
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::CountCandidates(const TmPt &tids)
{
   size_t nCandTofTracks = 0, nCandTrueMatch = 0;

   for (const auto &iter : tids) // cycle by TPC KF tracks(Pt descending)
   {
      Int_t KfIndex      = iter.second;
      bool  HadTofSignal = false, IsTrueMatching = false;

      for (size_t hitId = 0; hitId < fNHits; hitId++) // cycle by hits of KfIndex track
      {
         long uid = hash(KfIndex, hitId);

         auto it = mCandidates.find(uid); // pair< hash(trackId, hitId), MpdTofMatchingData >
         if (it != mCandidates.end()) {
            if (it->second.fHadTofSignal) HadTofSignal = true;
            if (it->second.fIsTrueMatching) IsTrueMatching = true;
         }
      }
      if (HadTofSignal) nCandTofTracks++;
      if (IsTrueMatching) nCandTrueMatch++;

   } // cycle by TPC KF tracks

   fNCandTrueMatch += nCandTrueMatch;
   fNCandTofTracks += nCandTofTracks;
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::Process(TmPt tids, const TClonesArray *aKfTracks)
{
   CountCandidates(tids);

   for (const auto &iter : tids) // cycle by TPC KF tracks(Pt descending)
   {
      Int_t KfIndex = iter.second;
      //	auto pKfTrack = (const MpdTpcKalmanTrack*) iter.first;
      //	auto mcInfo = getMCinfo(KfIndex);

      size_t hitIndex[50];
      size_t rank = GetRank(KfIndex, hitIndex);

      double                bestParameter = 1.E+5; // big value
      TmCandidate::iterator bestCand;
      size_t                bestHitId = -1;

      for (size_t i = 0; i < rank; i++) {
         auto it = mCandidates.find(hash(KfIndex, hitIndex[i])); // pair<hash(kf trackId, hitId), MpdTofMatchingData>
         if (it != mCandidates.end()) {
            double delta = it->second.GetDelta();

            if (delta < bestParameter) {
               bestCand      = it;
               bestParameter = delta;
               bestHitId     = hitIndex[i];
            }
         } else
            assert(false); // == "ERROR: MpdTofMatchingData candidate MUST exist."
      }

      if (bestHitId != -1) // accepted candidate
      {
         // set flag of accepted candidate.
         bestCand->second.fBestParameter = true;

         // erase this hit at all competitive candidates
         RemoveHit(bestHitId);
      }

   } // cycle by TPC KF tracks

   // clean rejected candidates
   for (auto cand = mCandidates.cbegin(), candEnd = mCandidates.cend();
        cand != candEnd;) // pair< hash(kftrackId, hitId), MpdTofMatchingData >
   {
      // fill debug histos
      if (pMatchingQA) pMatchingQA->FillHitDeviation(cand->second);

      if (cand->second.fBestParameter == false) {
         cand = mCandidates.erase(cand);
      } else
         cand++;
   }

   CountFoundCandidates();
}
//------------------------------------------------------------------------------------------------------------------------
void LWeightMatrix::MoveEntries(TClonesArray *aTofMatchings) const
{
   for (const auto &cit : mCandidates) {
      new ((*aTofMatchings)[aTofMatchings->GetEntriesFast()]) MpdTofMatchingData(cit.second);
   }
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdTofMatching);
//------------------------------------------------------------------------------------------------------------------------
MpdTofMatching::MpdTofMatching(const char *name, Int_t verbose, Bool_t test, const char *flnm)
   : FairTask(name, verbose), fDoTest(test)
{
   pRandom = new TRandom2;

   if (fDoTest) {
      pQA    = new MpdTofMatchingQA(flnm, false);
      pTimer = new TStopwatch;
   }

   fWeights = new LWeightMatrix(pQA);
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofMatching::~MpdTofMatching()
{
   delete pTimer;
   delete fWeights;
   delete pRandom;
   delete pQA;
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdTofMatching::Init()
{
   LOG(DEBUG2) << "[MpdTofMatching::Init] Begin initialization.";

   aMcPoints    = (TClonesArray *)FairRootManager::Instance()->GetObject("TOFPoint");
   aMcTracks    = (TClonesArray *)FairRootManager::Instance()->GetObject("MCTrack");
   aTofHits     = (TClonesArray *)FairRootManager::Instance()->GetObject("TOFHit");
   aTPCkfTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("TpcKalmanTrack");

   if (aMcPoints && aMcTracks) fIsMcRun = true;

   if (!aTofHits || !aTPCkfTracks) {
      LOG(ERROR) << "Branch not found! (TOFHit or TpcKalmanTrack)";
      return kERROR;
   }

   pKF = MpdKalmanFilter::Instance("KF", "KF");

   // Create and register output array
   aTofMatchings = new TClonesArray("MpdTofMatchingData");
   FairRootManager::Instance()->Register("TOFMatching", "Tof", aTofMatchings, kTRUE);

   MpdTofGeoUtils::Instance()->ParseTGeoManager(nullptr, false);
   MpdTofGeoUtils::Instance()->FindNeighborStrips(
      0.8, nullptr, false); // 0.8 [cm] <--- thresh. distance between neighbor strips,  (see h1TestDistance histo)

   LOG(DEBUG) << "[MpdTofMatching::Init] Initialization finished succesfully.";
   return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatching::Finish()
{
   if (pQA) pQA->Finish();
   if (pTimer) {
      Printf("\nMpdTofMatching TStopwatch report: ");
      pTimer->Print("m");
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatching::Exec(Option_t *option)
{
   if (pTimer) pTimer->Start(false);

   fDoMCtest = fDoTest && fIsMcRun;

   // Reset event
   aTofMatchings->Clear();

   Int_t nTofPoints = -1, nMCTracks = -1;
   if (fIsMcRun) {
      nTofPoints = aMcPoints->GetEntriesFast();
      nMCTracks  = aMcTracks->GetEntriesFast();
   }
   Int_t nTofHits  = aTofHits->GetEntriesFast();
   Int_t nKFTracks = aTPCkfTracks->GetEntries();

   LOG(DEBUG2) << "[MpdTofMatching::Exec] points= " << nTofPoints << ", hits= " << nTofHits
               << ", mc tracks= " << nMCTracks << ", kf tracks= " << nKFTracks;
   // ---------------------------------------------------------------------------------------->>> Select (! propagated
   // to Etof) & sort by Pt  the KfTpcTracks
   std::set<Int_t> sEndcapTracks;
   // ----------------------------------------------------------------------------------------
   TmPt mPt;
   for (Int_t index = 0; index < nKFTracks; index++) // cycle by TPC KF tracks
   {
      if (sEndcapTracks.find(index) != sEndcapTracks.end()) continue; // postpone for matching with ETof

      auto pKfTrack = (MpdTpcKalmanTrack *)aTPCkfTracks->UncheckedAt(index);
      mPt.insert(make_pair(pKfTrack, index));
   }

   // ---------------------------------------------------------------------------------------->>> Sorting & Mapping
   // points to MC tracks
   TmmT2P mmT2P; // pair< MCtrackID, MpdTofPoint*>
   if (fDoMCtest) MapMcT2P(aMcPoints, mmT2P);

   // ---------------------------------------------------------------------------------------->>> Mapping hits to
   // detectorUIDs
   TmmD2H mmD2H; // pair< duid, pair<MpdTofHit*, hitIndex> >
   TmS2H  mS2H;  // pair< suid, pair<MpdTofHit*, hitIndex> >
   TmmT2H mmT2H; // pair< tid, pair<MpdTofHit*, hitIndex> >

   vector<Int_t> v;
   for (Int_t hitIndex = 0; hitIndex < nTofHits; hitIndex++) // cycle by tof hits
   {
      auto  pTofHit = (MpdTofHit *)aTofHits->At(hitIndex);
      Int_t suid    = pTofHit->GetDetectorID();

      mmD2H.insert(make_pair((suid & 0xFFFF0000),
                             make_pair(pTofHit, hitIndex))); // convert suid to duid (reset stripID to 0 and gap to 0)
      mS2H.insert(make_pair(suid, make_pair(pTofHit, hitIndex)));

      pTofHit->getLinks(MpdTofUtils::mcTrackIndex, v); // vector<tid>
      for (auto it = v.begin(), itEnd = v.end(); it != itEnd; it++)
         mmT2H.insert(make_pair(*it, make_pair(pTofHit, hitIndex)));
   }
   // ----------------------------------------------------------------------------------------
   if (fDoMCtest) pQA->FillIdealMatching(mPt, mmT2H, aMcTracks);

   // resize matrix & clean by zero weight
   fWeights->Reset(mPt.size(), mS2H.size());

   // fill matching weight matrix.
   FillWeights(mPt, mmD2H, mS2H, mmT2P);

   if (fDoMCtest)
      pQA->FillNtracks(aMcTracks, fNTofTracksEvent); // MUST be call after FillWeights(calc. fNTofTracksEvent)

   // fWeights->Print("BEFORE", aTPCkfTracks, aTofHits, &mmT2H);

   //	fWeights->Boost1thRankCand(1000.);
   //	fWeights->Normalize();
   /// fWeights->Print("AfterNormalize", aTPCkfTracks, aTofHits, &mmT2H);

   fWeights->Process(mPt, aTPCkfTracks);
   // fWeights->Print("AFTER Process", aTPCkfTracks, aTofHits, &mmT2H);

   fWeights->MoveEntries(aTofMatchings);

   // MpdTof::Dump(aMcPoints, aTofHits, aMcTracks);

   if (fDoMCtest) pQA->FillMatchingEfficiency(mPt, aTofMatchings, aTofHits, aMcTracks);

   // pQA->tidsTofTouch.Print("\n tidsTofTouch ------------------------>>>");

   set<int> matchedTracks; // MCtid
   for (Int_t i = 0, N = aTofMatchings->GetEntriesFast(); i < N; i++) {
      auto match    = (MpdTofMatchingData *)aTofMatchings->At(i);
      int  KFtid    = match->GetKFTrackIndex();
      auto pKfTrack = (MpdTpcKalmanTrack *)aTPCkfTracks->UncheckedAt(KFtid);

      matchedTracks.insert(pKfTrack->GetTrackID());
   }

   set<int> kalmanTracks; // MCtid
   for (Int_t i = 0, N = aTPCkfTracks->GetEntriesFast(); i < N; i++) {
      auto pKfTrack = (MpdTpcKalmanTrack *)aTPCkfTracks->UncheckedAt(i);

      kalmanTracks.insert(pKfTrack->GetTrackID());
   }

   TVector3 momentum;
   for (Int_t tid = 0, Ntracks = aMcTracks->GetEntriesFast(); tid < Ntracks; tid++) // cycle by mc tracks
   {
      // if(()->GetNPoints(kTOF)) NmcTracks++;

      auto track = (MpdMCTrack *)aMcTracks->UncheckedAt(tid);

      if (track->GetMotherId() == -1) // primary mc track
      {
         bool IsKalman  = (kalmanTracks.find(tid) != kalmanTracks.end());
         bool IsMatched = (matchedTracks.find(tid) != matchedTracks.end());

         track->GetMomentum(momentum);

         if (pQA) pQA->FillTrackEfficiency(IsKalman, IsMatched, track->GetPdgCode(), momentum);
      }
   }

   Report(fVerbose);

   if (pTimer) pTimer->Stop();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatching::Report(Int_t verbose, ostream &os) const
{
   auto prec = os.precision(3);

   if (verbose > 0) {
      Int_t nMatchings = aTofMatchings->GetEntries();
      os << "\n -I- [MpdTofMatching::Exec] nTofTracks=" << fNTofTracksEvent << ", Matchings = " << nMatchings;
      if (fNTofTracksEvent != 0) os << "  (" << ((double)nMatchings) / fNTofTracksEvent * 100. << "%)";
   }

   if (verbose > 1) {
      auto cand  = fWeights->GetMatchCand();
      auto found = fWeights->GetMatchFound();
      if (fNTofTracksRun != 0) {
         os << "\n Cand.: tofTracks = (" << ((double)cand.first) / fNTofTracksRun * 100. << "%), TrueMatch = ("
            << ((double)cand.second) / fNTofTracksRun * 100. << "%)";
         os << "; Found: TrueMatch = (" << ((double)found.first) / fNTofTracksRun * 100. << "%), MisMatch = ("
            << ((double)found.second) / fNTofTracksRun * 100. << "%)";
      }
   }

   os.precision(prec);
}
//------------------------------------------------------------------------------------------------------------------------

#include "MpdTofGeoUtils.h"
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatching::FillWeights(const TmPt &mPt, const TmmD2H &mmD2H, const TmS2H &mS2H, const TmmT2P &mmT2P)
{
   const double Zerror   = 10.; // 10. [cm] // FIXME:  should be dependent on the parameters of the KFtrack
   const double PhiError = 0.1; // 0.1 // [rads] // FIXME:  should be dependent on the parameters of the KFtrack

   // if(pQA) pQA->tidsTofTouch.Reset();

   mcInfo *pMCdata = {};
   double  trLength;       //, trackLength;
                           // Int_t 		charge;
   vector<Tinterval> list; // list of strips overlapped probe point
   fNTofTracksEvent = 0;

   for (const auto &it : mPt) // cycle by TPC KF tracks
   {
      auto  pKfTrack = (const MpdTpcKalmanTrack *)it.first;
      Int_t KfIndex  = it.second;
      // Int_t tid = pKfTrack->GetTrackID();

      if (fDoMCtest) {
         pMCdata = fWeights->saveMCinfo(KfIndex, GetMcInfo(mmT2P, pKfTrack, aMcTracks));
         pQA->FillKfTracks(pMCdata->TofTouch, pKfTrack->Momentum3().Eta(), pKfTrack->Momentum());
      }

      if (pMCdata && pMCdata->TofTouch) fNTofTracksEvent++;

      TVector3 estPointOnCyl = EstTrackOnR(pKfTrack); // Estimate point on cylinder

      if (MpdTofGeoUtils::Instance()->IsPointInsideDetectors(estPointOnCyl, list, Zerror,
                                                             PhiError)) // !empty detector list
      {
         if (pQA) pQA->FillPointInsideDetector(estPointOnCyl);

         MpdTpcKalmanTrack ReFittedTrack(RefitTrack(pKfTrack));

         for (const auto &detector : list) // cycle by overlapped detector list
         {
            Int_t duid      = detector.value->volumeUID;
            auto  detCenter = detector.value->center;
            auto  detPerp   = detector.value->perp;

            // Propagate the KF track to the Tof detector plane; return true, if propagation successful
            TVector3 estPointOnPlate, hitPosition, Momentum;
            if (0 == EstTrackOnPlate(ReFittedTrack, detCenter, detPerp, estPointOnPlate, trLength,
                                     Momentum)) // Estimate point on detector plate
            {
               /*					// Propagate the KF  marker tracks to the detector plane; return true, if all propagations
                  successful bool retvalue = true; TVector3 est[4], P[4]; for(size_t i=0;i<4;i++) retvalue = retvalue &&
                  EstMarkerTrackOnPlate(i, pKfTrack, detCenter, detPerp, est[i], P[i]);

                              if(pQA) pQA->FillSmearedPoints(estPointOnPlate, est, P);
               */
               if (pMCdata && pQA) {
                  TVector3 pos, mom;
                  if (pMCdata->GetClosestPointOnDetector(pos, mom, estPointOnPlate, duid).first != -1)
                     pQA->FillMcPointDeviation(pos, estPointOnPlate, mom, Momentum, ReFittedTrack.GetNofHits());
               }

               auto range = mmD2H.equal_range(duid);
               for (auto iter = range.first; iter != range.second; iter++) // cycle by  estimated detector hits
               {
                  MpdTofHit *pHit     = iter->second.first;
                  Int_t      hitIndex = iter->second.second;
                  pHit->Position(hitPosition);

                  double devHit, devHitZ, devHitR, devHitPhi;
                  MpdTof::GetDelta(hitPosition, estPointOnPlate, devHit, devHitZ, devHitR, devHitPhi);

                  devHitZ   = abs(devHitZ);
                  devHitPhi = abs(devHitPhi);

                  if (devHitZ < fThreshZ && devHitPhi < fThreshPhi) // inside matching window
                  {
                     double weight = 1. / (estPointOnPlate - hitPosition).Mag();

                     fWeights->AddWeight(KfIndex, hitIndex, weight);
                     auto ptr =
                        fWeights->AddCandidate(MpdTofMatchingData(KfIndex, hitIndex, weight, pHit, trLength,
                                                                  pKfTrack->GetNofTrHits(), Momentum, estPointOnPlate));

                     if (pMCdata) {
                        ptr->fIsTrueMatching = pMCdata->IsSameStrip(pHit->GetDetectorID());
                        ptr->fHadTofSignal   = pMCdata->TofTouch;
                     }

                  } // inside matching window
               }
            }
         }
      }

   } // cycle by TPC KF tracks

   fNTofTracksRun += fNTofTracksEvent;
}
//------------------------------------------------------------------------------------------------------------------------
mcInfo MpdTofMatching::GetMcInfo(const TmmT2P &mmT2P, const MpdKalmanTrack *pKfTrack, TClonesArray *aMcTracks)
{
   mcInfo mcData;

   mcData.trackIndex = pKfTrack->GetTrackID();
   auto mcTrack      = (const MpdMCTrack *)aMcTracks->UncheckedAt(mcData.trackIndex);

   mcData.pid = mcTrack->GetPdgCode();
   auto track = const_cast<MpdMCTrack *>(mcTrack);
   track->GetMomentum(mcData.vertMomentum);
   track->GetStartVertex(mcData.vertPosition);

   auto range = mmT2P.equal_range(mcData.trackIndex); //	pair< MCtrackID, MpdTofPoint*>
   for (auto it = range.first; it != range.second; ++it)
      mcData.vPoints.push_back(it->second); // cycle by mcPoints for trackIndex mcTrack

   mcData.TofTouch = !mcData.vPoints.empty();

   return mcData;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatching::MapMcT2P(TClonesArray *aPoints, TmmT2P &mmT2P) const //  TmmT2P = pair< MCtrackID, MpdTofPoint*>
{
   for (Int_t index = 0, nTofPoints = aPoints->GetEntriesFast(); index < nTofPoints; index++) // cycle by MpdTofPoint
   {
      auto mcPoint = (MpdTofPoint *)aPoints->UncheckedAt(index);

      Int_t    tid  = mcPoint->GetTrackID();
      Double_t time = mcPoint->GetTime();

      auto iter = mmT2P.find(tid);
      if (iter != mmT2P.end()) // same trackID already inserted, insert to position(sorting by time)
      {
         int count = mmT2P.count(tid);
         for (int i = 0; i < count; i++, iter++) // cycle by hits with same trackID
         {
            if (time < iter->second->GetTime()) {
               mmT2P.insert(iter, make_pair(tid, mcPoint));
               break;
            }

            if (i == count - 1) mmT2P.insert(++iter, make_pair(tid, mcPoint)); // insert to last
         }
      } else
         mmT2P.insert(make_pair(tid, mcPoint));

   } // cycle by MpdTofPoint
}
//------------------------------------------------------------------------------------------------------------------------
TVector3 MpdTofMatching::EstTrackOnR(const MpdTpcKalmanTrack *tr) const
{
   MpdKalmanHit hEnd;
   hEnd.SetType(MpdKalmanHit::kFixedR);
   MpdTpcKalmanTrack tr1(*tr);
   TObjArray *       hits = tr1.GetHits();

   if (hits->GetEntriesFast() == 0) {
      Int_t nh = tr1.GetNofTrHits();
      for (Int_t j = 0; j < nh; ++j) hits->Add(tr1.GetTrHits()->UncheckedAt(j));
      tr1.SetNofHits(nh);
   }

   tr1.SetParam(*tr1.GetParamAtHit());
   tr1.SetParamNew(*tr1.GetParamAtHit());
   tr1.SetWeight(*tr1.GetWeightAtHit());
   tr1.SetLength(tr1.GetLengAtHit());
   auto h = (MpdKalmanHit *)tr1.GetTrHits()->First();
   tr1.SetPos(tr1.GetPosAtHit());

   if (h->GetType() == MpdKalmanHit::kFixedZ) tr1.SetType(MpdKalmanTrack::kEndcap);
   tr1.SetPosNew(tr1.GetPos());

   hEnd.SetPos(fTofBarrelRadius); // barrel TOF radius, [cm]
   pKF->PropagateToHit(&tr1, &hEnd, kTRUE);

   TVector3 pos(tr1.GetPosNew(), 0., 0.);
   pos.SetPhi(tr1.GetParamNew(0) / tr1.GetPosNew());
   pos.SetZ(tr1.GetParamNew(1));

   return pos;
}
//------------------------------------------------------------------------------------------------------------------------
MpdTpcKalmanTrack MpdTofMatching::RefitTrack(const MpdTpcKalmanTrack *tr) const
{
   MpdTpcKalmanTrack tr1(*tr);
   tr1.SetDirection(MpdKalmanTrack::kOutward);
   TObjArray *hits = tr1.GetHits();

   if (hits->GetEntriesFast() == 0) {
      Int_t nh = tr1.GetNofTrHits();
      for (Int_t j = 0; j < nh; ++j) hits->Add(tr1.GetTrHits()->UncheckedAt(j));
      tr1.SetNofHits(nh);
   }

   tr1.SetParam(*tr1.GetParamAtHit());
   tr1.SetParamNew(*tr1.GetParamAtHit());
   tr1.SetWeight(*tr1.GetWeightAtHit());
   tr1.SetLength(tr1.GetLengAtHit());
   tr1.SetPos(tr1.GetPosAtHit());
   tr1.SetPosNew(tr1.GetPos());

   return tr1;
}
//------------------------------------------------------------------------------------------------------------------------
int MpdTofMatching::EstTrackOnPlate(const MpdTpcKalmanTrack &tr, const TVector3 &point, const TVector3 &perp,
                                    TVector3 &pos, Double_t &length, TVector3 &P) const
{
   MpdTpcKalmanTrack tr1(tr);

   Double_t plane[6] = {point.X(), point.Y(), point.Z(), perp.X(), perp.Y(), perp.Z()};
   if (!pKF->PropagateParamP(&tr1, plane, kTRUE)) return -1;

   Double_t R = tr1.GetPosNew();
   if (TMath::IsNaN(R)) return -2;

   length = tr1.GetLength();
   P      = tr1.Momentum3();

   pos.SetXYZ(R, 0., 0.);
   pos.SetPhi(tr1.GetParamNew(0) / R);
   pos.SetZ(tr1.GetParamNew(1));

   return 0;
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofMatching::EstMarkerTrackOnPlate(size_t mode, const MpdTpcKalmanTrack *tr, const TVector3 &point,
                                           const TVector3 &perp, TVector3 &estPoP, TVector3 &P) const
{
   const static Double_t fError_dZ = 0.03, fError_dPhi = 0.03; // [rads]  msc

   auto      tr1   = RefitTrack(tr);
   TMatrixD &param = *tr1.GetParamNew();

   switch (mode) {
   case 0: param(3, 0) += fError_dZ; break;

   case 1: param(2, 0) -= fError_dPhi; break;
   case 2: param(3, 0) -= fError_dZ; break;

   case 3: param(2, 0) += fError_dPhi; break;

   default:;
   };

   tr1.SetParam(param);

   Double_t plane[6] = {point.X(), point.Y(), point.Z(), perp.X(), perp.Y(), perp.Z()};
   if (!pKF->PropagateParamP(&tr1, plane, kTRUE)) return false;

   Double_t R = tr1.GetPosNew();
   if (TMath::IsNaN(R)) return false;

   P = tr1.Momentum3();

   estPoP.SetXYZ(R, 0., 0.);
   estPoP.SetPhi(tr1.GetParamNew(0) / R);
   estPoP.SetZ(tr1.GetParamNew(1));

   return true;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatching::Print(const MpdTpcKalmanTrack *track, const char *comment, ostream &os) const
{
   if (comment != nullptr) os << comment;

   Double_t rad = track->GetPosNew();
   Double_t phi = track->GetParamNew(2);
   double   X   = rad * TMath::Cos(phi); // X
   double   Y   = rad * TMath::Sin(phi); // Y
   double   Z   = track->GetParamNew(1); // Z
   TVector3 P   = track->Momentum3();
   os << "Mom=(" << P.X() << "," << P.Y() << "," << P.Z() << ";" << P.Perp() << "," << P.Mag() << ") Pos=(" << X << ","
      << Y << "," << Z << "; " << sqrt(X * X + Y * Y) << "," << sqrt(X * X + Y * Y + Z * Z) << ")";
}
//------------------------------------------------------------------------------------------------------------------------
