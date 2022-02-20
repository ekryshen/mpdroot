//------------------------------------------------------------------------------------------------------------------------
/// \class MpdEtofMatching
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <map>

#include <TMath.h>
#include <TFile.h>
#include <TRandom2.h>
#include <TStorage.h>
#include <TEfficiency.h>

#include "FairRootManager.h"
#include "MpdMCTrack.h"
#include "FairRunAna.h"
#include "FairLogger.h"

#include "MpdEctKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdTofPoint.h"
#include "MpdTofHit.h"
#include "MpdTofMatching.h"
#include "MpdTofGeoUtils.h"
#include "MpdEtofGeoUtils.h"
#include "MpdEtof.h"
#include "MpdTofMatchingQA.h"

#include "MpdEtofMatching.h"

using namespace std;
/*
struct less_by_pointer
{
    inline bool operator() (const intervalType& struct1, const intervalType& struct2)
    {
        return (struct1.value < struct2.value);
    }
};
*/
ClassImp(MpdEtofMatching);
//------------------------------------------------------------------------------------------------------------------------
MpdEtofMatching::MpdEtofMatching(const char *name, Int_t verbose, Bool_t test, const char *flnm)
   : FairTask(name, verbose), aMcPoints(nullptr), aMcTracks(nullptr), aTofHits(nullptr), aKFectTracks(nullptr),
     aTofMatchings(nullptr), fMode(kIntervalTree), fDoTest(test), fUseMCData(false), pRandom(new TRandom2),
     fNSmeared(20), fTofZpos(295.2), fTofRmax(140.), fThreshR(15.), fThreshTheta(1.5), // Z- 250.0
     pKF(nullptr)
{
   pMatchingQA = fDoTest ? new MpdTofMatchingQA(flnm, true) : nullptr;
   // pMF = new LMatchingFilter(pMatchingQA, fVerbose);
}
//------------------------------------------------------------------------------------------------------------------------
MpdEtofMatching::~MpdEtofMatching()
{
   //    delete pMF;
   delete pRandom;
   delete pMatchingQA;
   if (aTofMatchings != nullptr) {
      aTofMatchings->Delete();
      delete aTofMatchings;
   }
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdEtofMatching::Init()
{
   aMcPoints    = (TClonesArray *)FairRootManager::Instance()->GetObject("ETOFPoint");
   aMcTracks    = (TClonesArray *)FairRootManager::Instance()->GetObject("MCTrack");
   aTofHits     = (TClonesArray *)FairRootManager::Instance()->GetObject("ETOFHit");
   aKFectTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("EctTrack");

   if (aMcPoints && aMcTracks) fUseMCData = true;

   if (!aTofHits || !aKFectTracks) {
      LOG(FATAL) << "[MpdEtofMatching::Init] Branch not found!";
      return kERROR;
   }

   pKF = MpdKalmanFilter::Instance("KF", "KF");

   // Create and register output array
   aTofMatchings = new TClonesArray("MpdTofMatchingData");
   FairRootManager::Instance()->Register("ETOFMatching", "ETof", aTofMatchings, kTRUE);

   ///	pMF->SetContainer(aTofMatchings);

   //	MpdEtofGeoUtils::Instance()->ParseTGeoManager(fUseMCData, nullptr, false);
   //	MpdEtofGeoUtils::Instance()->FindNeighborStrips(0.8, nullptr, nullptr, false);// 0.8 [cm] <--- thresh. distance
   // between neighbor strips,  (see h1TestDistance histo)

   LOG(INFO) << "MpdEtofMatching initialization finished succesfully.";

   return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofMatching::Exec(Option_t *option)
{
   fDoMCTest = fUseMCData && fDoTest;
   /*
      // Reset event
           aTofMatchings->Clear();
      pMF->Reset();

           Int_t nTofPoints = -1, nMCTracks = -1;
      if(fUseMCData){ nTofPoints = aMcPoints->GetEntriesFast(); nMCTracks = aMcTracks->GetEntriesFast();}
      Int_t nTofHits = aTofHits->GetEntriesFast();
      Int_t nKFTracks = aKFectTracks->GetEntriesFast();

           if(fVerbose) cout<<"\n -I- [MpdEtofMatching::Exec] points= "<<nTofPoints<<", hits= "<<nTofHits<<", mc tracks=
   "<<nMCTracks<<", kf tracks= "<<nKFTracks<<endl;

      // ---------------------------------------------------------------------------------------->>> Sorting & Mapping
   points to MC tracks TmmP2T 			mmMCpoints;	// pair< MCtrackID, MpdTofPoint*> TmmP2T::const_iterator
   mmMCpointCiter;

      if(fUseMCData) MpdTofMatching::MappingMcPoints2McTracks(aMcPoints, mmMCpoints);

      // ---------------------------------------------------------------------------------------->>> Mapping hits to
   stripUIDs typedef std::map<Int_t, std::pair<MpdTofHit*, Int_t> > 		TmS2H; // pair< stripUID, pair<MpdTofHit*,
   hitIndex> > TmS2H 			mHits;

           for(Int_t hitIndex = 0; hitIndex < nTofHits; hitIndex++) // cycle by tof hits
      {
         MpdTofHit *pTofHit = (MpdTofHit*) aTofHits->At(hitIndex);
         mHits.insert(make_pair( pTofHit->GetDetectorID(), make_pair(pTofHit, hitIndex)));
      }
      // ---------------------------------------------------------------------------------------->>> Sorting by Pt  the
   KfEctTracks TsPt		tids; for(Int_t KfIndex = 0; KfIndex < nKFTracks; KfIndex++) 	// cycle by ECT KF tracks
      {
         MpdEctKalmanTrack *pKfTrack = (MpdEctKalmanTrack*) aKFectTracks->UncheckedAt(KfIndex);
         tids.insert(make_pair(pKfTrack, KfIndex));
      }
      // ----------------------------------------------------------------------------------------
      const intervalTreeType*	mDetectorsR = MpdEtofGeoUtils::Instance()->GetDetR();
      const intervalTreeType*	mDetectorsPhi = MpdEtofGeoUtils::Instance()->GetDetPhi();

           vector<intervalType> 	segmentR, segmentPhi, intersect;

           TVector3 	hitPosition, Momentum, estPoint;
           bool 		IsInside;
           Double_t 	trackLength, thR, thPhi;
         Int_t  		charge, selectedTracks  = 0;

      mcInfo 		MCdata; // The MC run variables. VALID ONLY if fUseMCData = true

   //cout<<"\n
   ------------------------------------------------------------------------------------------------------------->>
   EVENT";
   //mDetectorsR->dump("\n\n ----->>>	mDetectorsR INTERVALS");
   //mDetectorsPhi->dump("\n\n ----->>>          mDetectorsPhi INTERVALS");

      TmS2H::const_iterator	hitCiter;
      for(TsPt::const_iterator it = tids.begin(), itEnd = tids.end(); it != itEnd; it++)// cycle by sorted ECT KF tracks
      {
         const MpdEctKalmanTrack *pKfTrack = (MpdEctKalmanTrack*) it->first;
         Int_t KfIndex = it->second;

         if(fUseMCData) MpdTofMatching::GetMcInfo(MCdata, mmMCpoints, pKfTrack, aMcTracks);

         size_t TofSide = (pKfTrack->Momentum3().Pz() > 0.) ? 0 : 1;
         double TofZpos = (TofSide == 0) ? fTofZpos : -fTofZpos;

         if(EstTrackOnPlane(pKfTrack, TofZpos, estPoint, trackLength, Momentum, charge))
         {
            IsInside = (estPoint.Perp() < fTofRmax)  ? true : false;

            if(fDoMCTest) pMatchingQA->FilleRest(estPoint.Perp());
         }

         if(!IsInside) continue;	// KF track out ETof Rmax

         if(fDoMCTest)
         {
            if(MCdata.Npoints == 1) // only SINGLE tof point per MCtrack
               pMatchingQA->FillDevPoint2EstCyl(MCdata.Position, estPoint, pKfTrack->Momentum());

            if(MCdata.TofTouch) // KF track have TOFpoint
            {
               selectedTracks++;
               pMatchingQA->FillSelectedKfTracks(pKfTrack->Momentum3().Eta(), pKfTrack->Momentum());
            }
         }

         // ---------------------------------------------------------------------------------------->>> Looking for
   overlaping of estPoint & detectors double estR = estPoint.Perp(), estPhi = estPoint.Phi();

         const static double Rerror = 10.; // [cm] // FIXME:  should be dependent on the parameters of the KFtrack
         const static double PhiError = 0.1; // [rads] // FIXME:  should be dependent on the parameters of the KFtrack

         segmentR.clear();
         segmentPhi.clear();
         intersect.clear();

         mDetectorsR[TofSide].findOverlapping(estR - Rerror,  estR + Rerror, segmentR);
         mDetectorsPhi[TofSide].findOverlapping(estPhi - PhiError, estPhi + PhiError, segmentPhi);

   //cout<<"\n <<<--R-->>>>> ("<<estR-Rerror<<", "<<estR+Rerror<<")";
   //MpdTof::intervalTreeType::dumpIntervals(&segmentR, "\n\n ----->>>	mDetectorsR findOverlapping");

   //cout<<"\n <<<--Phi-->>>>> ("<<estPhi-PhiError<<", "<<estPhi+PhiError<<")";
   //MpdTof::intervalTreeType::dumpIntervals(&segmentPhi, "\n\n ----->>>	mDetectorsPhi findOverlapping");

         if(!segmentR.empty() && !segmentPhi.empty()) // have overlaped segments both R and Phi
         {
            // calc. intersection
            sort(segmentR.begin(), segmentR.end(), less_by_pointer());
               sort(segmentPhi.begin(), segmentPhi.end(), less_by_pointer());
            set_intersection(segmentR.begin(), segmentR.end(), segmentPhi.begin(), segmentPhi.end(),
   std::back_inserter(intersect), less_by_pointer());

            for(vector<intervalType>::const_iterator cit = intersect.begin(), citEnd = intersect.end(); cit != citEnd;
   cit++) // cycle by the overlaped strips
            {
               Int_t detUID = (*cit).value->volumeUID;

               if(fDoMCTest && MCdata.Npoints == 1) // only one tof point per MCtrack
                  pMatchingQA->FillDevPoint2EstP(MCdata.Position, estPoint, pKfTrack->Momentum());

               hitCiter = mHits.find(detUID);
               if(hitCiter != mHits.end()) // the estimated strip have hit
               {
                  MpdTofHit *TofHit = hitCiter->second.first;
                  Int_t hitIndex = hitCiter->second.second;
                  TofHit->Position(hitPosition);

                  double devHit, devHitZ, devHitR, devHitPhi;
                  MpdTofMatching::GetDelta(hitPosition, estPoint, devHit, devHitZ, devHitR, devHitPhi);

                  bool mcIsSameVolumeUID = (TofHit->GetDetectorID() == MCdata.stripUID);

                  if(fDoMCTest) pMatchingQA->FillHitDev2EstP(mcIsSameVolumeUID, devHitR, devHitPhi);

                  if(devHitR < fThreshR && devHitPhi < fThreshTheta) // inside matching window
                  {
                     pMF->AddCandidate(MpdTofMatchingData(KfIndex, hitIndex, pKfTrack->GetNofTrHits(), TofHit,
   MCdata.pid, TofHit->GetFlag(), trackLength, estPoint, Momentum, charge, devHitR,  devHitPhi));

   //cout<<"\n AddCandidate	KfIndex="<<KfIndex<<" hitIndex="<<hitIndex<<" delta="<<delta<<"
   NmbTrHits="<<pKfTrack->GetNofTrHits()<<"  mcHasCand="<<mcHasCand;
                  }
               }

            } // cycle by the overlaped strips

         } // have overlaped segments both R and Phi

           if(fDoMCTest && MCdata.TofTouch) pMatchingQA->FillCandidates(MCdata.HaveTrueCand, MCdata.HaveCand,
   pKfTrack->Momentum3().Eta(), pKfTrack->Momentum());

      } // cycle by KF tracks

      if(fDoTest) pMatchingQA->FillTrackPerEvent(nKFTracks, selectedTracks);

      double chi2;
      Int_t MatchingOK = pMF->Processing(nKFTracks, chi2);	// accept candidates
      Int_t nEntries = pMF->UpdateContainer(); 		// remove unmatched candidates
   assert(nEntries	== MatchingOK);

      if(fDoMCTest) pMatchingQA->FillMatchingEfficiency(aTofMatchings, aTofHits, aKFectTracks, aMcTracks);

      if(fVerbose)
      {
         cout<<"\n -I- [MpdEtofMatching::Exec] MatchingOK = "<<MatchingOK;
         if(MatchingOK != 0) cout<<" ("<<chi2/MatchingOK<<")";
      }
   */
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofMatching::Finish()
{
   if (pMatchingQA) pMatchingQA->Finish();
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdEtofMatching::EstTrackOnPlane(const MpdEctKalmanTrack *tr, Double_t Zetof, TVector3 &pos, Double_t &length,
                                      TVector3 &Mom, Int_t &charge) const
{
   MpdEctKalmanTrack tr1(*tr);
   TObjArray        *hits = tr1.GetHits();

   if (hits->GetEntriesFast() == 0) {
      Int_t nh = tr1.GetNofTrHits();
      for (Int_t j = 0; j < nh; ++j) hits->Add(tr1.GetTrHits()->UncheckedAt(j));
      tr1.SetNofHits(nh);
   }

   MpdKalmanHit *h = (MpdKalmanHit *)tr1.GetTrHits()->First();
   if (h->GetType() == MpdKalmanHit::kFixedZ) tr1.SetType(MpdKalmanTrack::kEndcap);
   tr1.SetPos(tr1.GetPosAtHit());
   tr1.SetPosNew(tr1.GetPos());
   tr1.SetParamNew(*tr1.GetParamAtHit());
   tr1.SetParam(*tr1.GetParamAtHit());
   tr1.SetWeight(*tr1.GetWeightAtHit());
   tr1.SetLength(tr1.GetLengAtHit());
   tr1.SetDirection(MpdKalmanTrack::kOutward);

   MpdKalmanHit hEnd;
   hEnd.SetType(MpdKalmanHit::kFixedZ);
   hEnd.SetPos(Zetof); // eTOF Z position, [cm]

   if (!pKF->PropagateParamZ(&tr1, &hEnd, kTRUE)) return false;

   Double_t Z = tr1.GetPosNew();
   if (TMath::IsNaN(Z)) return false;

   Double_t Pt_inv = tr1.GetParamNew(4); // 1/Pt
   if (Pt_inv == 0.) return false;

   length = tr1.GetLength();
   Mom    = tr1.Momentum3();
   charge = tr1.Charge();
   pos.SetXYZ(tr1.GetParamNew(0), tr1.GetParamNew(1), Z);

   return true;
}
//------------------------------------------------------------------------------------------------------------------------