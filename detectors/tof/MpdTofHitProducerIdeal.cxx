//------------------------------------------------------------------------------------------------------------------------
#include <assert.h>
#include <iostream>

#include <TClonesArray.h>

#include "FairRootManager.h"
#include "FairLogger.h"

#include "MpdTofUtils.h"
#include "MpdTofHit.h"
#include "MpdTofPoint.h"
#include "MpdTofHitProducerQA.h"

#include "MpdTofHitProducerIdeal.h"

using namespace std;

ClassImp(MpdTofHitProducerIdeal);
//------------------------------------------------------------------------------------------------------------------------
MpdTofHitProducerIdeal::MpdTofHitProducerIdeal(const char *name, Bool_t useMCdata, Int_t verbose, Bool_t test,
                                               Bool_t merge, const char *flnm, bool IsEndcap)
   : FairTask(name, verbose), fDoTest(test), fDoMergeHits(merge), fUseMCData(useMCdata)
{
   pQA = fDoTest ? new MpdTofHitProducerQA(flnm, IsEndcap) : nullptr;
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofHitProducerIdeal::~MpdTofHitProducerIdeal()
{
   delete pQA;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerIdeal::AddParameters(TString &buf) const
{
   buf += "\n Run parameters:";
   buf += "  fDoTest=";
   buf += fDoTest;
   buf += ", fDoMergeHits=";
   buf += fDoMergeHits;
   buf += ", fUseMCData=";
   buf += fUseMCData;
   buf += ", fOnlyPrimary=";
   buf += fOnlyPrimary;
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdTofHitProducerIdeal::Init()
{
   if (fUseMCData) {
      aMcPoints = (TClonesArray *)FairRootManager::Instance()->GetObject("TOFPoint");
      aMcTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("MCTrack");
      assert(aMcPoints);
      assert(aMcTracks);
   } else {
      aExpDigits = (TClonesArray *)FairRootManager::Instance()->GetObject("??????"); // FIXME: NOW unknown name
      assert(aExpDigits);
   }

   // Create and register output array
   aTofHits = new TClonesArray("MpdTofHit");
   FairRootManager::Instance()->Register("TOFHit", "Tof", aTofHits, kTRUE);

   LOG(INFO) << "[MpdTofHitProducerIdeal::Init] Initialization finished succesfully.";

   return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerIdeal::Exec(Option_t *opt)
{
   static const TVector3 XYZ_err(0., 0., 0.);

   aTofHits->Clear();

   Int_t    nSingleHits = 0;
   TVector3 pos;

   if (fUseMCData) {
      for (Int_t pointIndex = 0, nTofPoint = aMcPoints->GetEntriesFast(); pointIndex < nTofPoint;
           pointIndex++) // cycle by TOF points
      {
         MpdTofPoint *pPoint = (MpdTofPoint *)aMcPoints->UncheckedAt(pointIndex);
         pPoint->Position(pos);

         AddHit(MpdTofPoint::ClearGap(pPoint->GetDetectorID()), pos, XYZ_err, pointIndex, pPoint->GetTrackID(),
                pPoint->GetTime(), MpdTofUtils::IsSingle);
         nSingleHits++;
      }
   } else //  exp. data used
   {
      // FIXME: now not realized
      // AddHit(Int_t detUID, const TVector3 &posHit, const TVector3 &posHitErr, Int_t expDigitIndex, Double_t time,
      // Int_t flag)
      assert(false);
   }

   int nFinally;
   if (fDoMergeHits) {
      MergeHitsOnStrip();        // save only the fastest hit in the strip
      nFinally = CompressHits(); // remove blank slotes
   } else
      nFinally = aTofHits->GetEntriesFast();

   LOG(DEBUG1) << "[MpdTofHitProducerIdeal::Exec] single hits= " << nSingleHits << ", final hits= " << nFinally;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerIdeal::Finish()
{
   if (pQA) pQA->Finish();
}
//------------------------------------------------------------------------------------------------------------------------
size_t MpdTofHitProducerIdeal::MergeHitsOnStrip(void)
{
   size_t                  mergedNmb = 0;
   map<Int_t, MpdTofHit *> mHits;  // pair<suid, MpdTofHit*> fastest hits map
   multiset<Int_t>         msUIDs; // suid for Hits
   MpdTofHit              *fastHit, *slowHit;

   for (Int_t hitIndex = 0, nHits = aTofHits->GetEntriesFast(); hitIndex < nHits; hitIndex++) // cycle by hits
   {
      auto pHit = (MpdTofHit *)aTofHits->UncheckedAt(hitIndex);
      assert(nullptr != pHit);
      Int_t suid = pHit->GetDetectorID();

      if (fDoTest) msUIDs.insert(suid);

      auto it = mHits.find(suid);
      if (it != mHits.end()) // hit for this suid already exist
      {
         mergedNmb++;

         if (pHit->GetTime() < it->second->GetTime()) //  faster hit found
         {
            fastHit = pHit;
            slowHit = it->second;
         } else {
            fastHit = it->second;
            slowHit = pHit;
         }

         fastHit->AddLinks(slowHit->GetLinks()); // copy links
         aTofHits->Remove(slowHit);              // remove old hit   --> make blank slote !!
         fastHit->AddFlag(MpdTofUtils::HaveTail);
         it->second = fastHit; // change pair value to current suid

         if (pQA) pQA->FillMergedTimes(fastHit->GetTime(), slowHit->GetTime());
      } else
         mHits.insert(make_pair(suid, pHit)); // insert new suid pair

   } // cycle by hits

   if (pQA) {
      for (auto it = msUIDs.begin(), itEnd = msUIDs.end(); it != itEnd; it = msUIDs.upper_bound(*it))
         pQA->FillOccupancy(msUIDs.count(*it));
   }

   return mergedNmb;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerIdeal::AddHit(Int_t suid, const TVector3 &posHit, const TVector3 &posHitErr, Int_t expDigitIndex,
                                    Double_t time, Int_t flag)
{
   MpdTofHit *pHit =
      new ((*aTofHits)[aTofHits->GetEntriesFast()]) MpdTofHit(suid, posHit, posHitErr, expDigitIndex, time, flag);

   pHit->AddLink(FairLink(MpdTofUtils::expDigitIndex, expDigitIndex));
   pHit->AddLink(FairLink(MpdTofUtils::volumeUID, suid));
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerIdeal::AddHit(Int_t suid, const TVector3 &posHit, const TVector3 &posHitErr, Int_t mcPointIndex,
                                    Int_t mcTrackIndex, Double_t time, Int_t flag)
{
   MpdTofHit *pHit =
      new ((*aTofHits)[aTofHits->GetEntriesFast()]) MpdTofHit(suid, posHit, posHitErr, mcPointIndex, time, flag);

   pHit->AddLink(FairLink(MpdTofUtils::mcPointIndex, mcPointIndex));
   pHit->AddLink(FairLink(MpdTofUtils::mcTrackIndex, mcTrackIndex));
   pHit->AddLink(FairLink(MpdTofUtils::volumeUID, suid));
}
//------------------------------------------------------------------------------------------------------------------------
Int_t MpdTofHitProducerIdeal::CompressHits()
{
   aTofHits->Compress();
   return aTofHits->GetEntriesFast();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerIdeal::Dump(const char *title, ostream &out) const
{
   out << "\n [MpdTofHitProducer::Dump]   ";
   if (title) out << title;
   out << ", size= " << aTofHits->GetEntriesFast();

   MpdTofPoint *point;
   MpdTofHit   *pHit;
   TVector3     hitPos, pointPos;

   TIterator *iter = aTofHits->MakeIterator();
   while ((pHit = (MpdTofHit *)iter->Next())) {
      pHit->Position(hitPos);
      out << "\n    hit suid = " << pHit->GetDetectorID() << ", hit pos(" << hitPos.X() << "," << hitPos.Y() << ","
          << hitPos.Z() << "), flag =" << pHit->GetFlag();

      if (aMcPoints) {
         point = (MpdTofPoint *)aMcPoints->UncheckedAt(pHit->GetRefIndex());
         point->Position(pointPos);
         out << "\n point suid = " << point->GetDetectorID() << ", point pos(" << pointPos.X() << "," << pointPos.Y()
             << "," << pointPos.Z() << "), dev=" << (hitPos - pointPos).Mag();
      }
   }

   delete iter;
}
//------------------------------------------------------------------------------------------------------------------------
