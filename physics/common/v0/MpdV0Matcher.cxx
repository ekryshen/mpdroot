/*
 * MpdV0Matcher.cxx
 *
 *  Created on: 27 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0Matcher.h"

#include <FairRootManager.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TObjArray.h>

#include "FairLogger.h"
#include "MpdMCTrack.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdMiniMcTrack.h"
#include "MpdV0Particle.h"
#include "MpdV0Track.h"

#include <iostream>
using std::cout;
using std::endl;

MpdV0Matcher::MpdV0Matcher(MpdCommonV0::EParticleType pidHypo, EMatchType type)
   : fWrite(kFALSE), fFirst(kFALSE), fPidHipo(pidHypo), fMethod(type), fMomThreshold(0.1), fMpdEvent(nullptr),
     fMiniEvents(nullptr), fMiniTracks(nullptr), fMcTracks(nullptr), fV0s(nullptr), fFormat(EFormat::kDst),
     fLinks(nullptr)
{
}

void MpdV0Matcher::MatchV0ByMomentum()
{
   fMcV0Indexes.clear();
   fMcPdg.clear();
   fMcV0Momenta.clear();
   switch (fFormat) {
   case EFormat::kDst: {
      PrepareMomentumMatchingDst();
   } break;
   case EFormat::kMiniDst: {
      PrepareMomentumMatchinMiniDst();
   } break;
   }
   Int_t matched = 0;
   for (int i = 0; i < fV0s->GetEntriesFast(); i++) {
      if (fLinks->GetLink(i) >= 0) continue;  // this track is already matched !
      if (fLinks->GetLink(i) == -2) continue; // this track is marked as bad !
      MpdV0Track *track = (MpdV0Track *)fV0s->UncheckedAt(i);
      if (track->GetPdg() != fMcV0sPid[i]) continue; // this v0 has different pid hypo then our v0
      if (!IsGoodV0(track->GetPdg())) continue;
      const TVector3 momTrack   = track->GetMomentum();
      Double_t       v0p2scaled = 1.0 / momTrack.Mag2();
      Double_t       momDiff    = 1E+9;
      Int_t          minIndex   = -1;
      for (int j = 0; j < fMcV0Momenta.size(); j++) {
         const TVector3 &momParticle = fMcV0Momenta[j];
         Double_t        diff        = (momParticle - momTrack).Mag2() * v0p2scaled;
         if (diff < momDiff) {
            minIndex = j;
            momDiff  = diff;
         }
      }
      momDiff = TMath::Sqrt(momDiff);
      if (momDiff > fMomThreshold) minIndex = -2; // mark as bad v0
      fLinks->SetLink(i, minIndex);
      if (minIndex >= 0) matched++;
   }
}

InitStatus MpdV0Matcher::Init()
{

   FairRootManager *mngr = FairRootManager::Instance();
   fMpdEvent             = (MpdEvent *)mngr->GetObject("MPDEvent.");
   fMiniEvents           = (TClonesArray *)mngr->GetObject("Event");
   fMiniTracks           = (TClonesArray *)mngr->GetObject("Track");
   fMcTracks             = (TClonesArray *)mngr->GetObject("MCTrack");
   if (mngr->GetObject("MpdV0")) {
      fV0s = (TClonesArray *)mngr->GetObject("MpdV0");
   } else {
      LOG(error) << "Lack of V0 particles!";
      return kFATAL;
   }

   /**
    * sprawdzic czy to jest pierwszy tak - tzn. czy branch istnieje żeby czyścić pamięć i reallokować indexy
    */
   if (mngr->CheckBranch("MpdV0Links") == 0) {
      fFirst = kTRUE;
      fLinks = new MpdSimpleLinks<Int_t>();
      mngr->Register("MpdV0Links", "V0", fLinks, fWrite);
   } else {

      fLinks = (MpdSimpleLinks<Int_t> *)mngr->GetObject("MpdV0Links");
   }

   if (fMpdEvent) {
      fFormat = EFormat::kDst;
   } else if (fMiniEvents) {
      fFormat = EFormat::kMiniDst;
   } else {
      LOG(ERROR) << "No MpdEvent/MiniDstEvent in tree";
      return kFATAL;
   }
   fMcV0sPid = MpdCommonV0::GetPdgs(fPidHipo);
   return kSUCCESS;
}

void MpdV0Matcher::Exec(Option_t *option)
{
   fLinks->Clear();
   if (fFirst) { // set links to -1
      for (int i = 0; i < fV0s->GetEntriesFast(); i++) {
         fLinks->PushBack(-1);
      }
   }
   switch (fMethod) {
   case EMatchType::kMatchByMomentum: {
      MatchV0ByMomentum();
   } break;
   case EMatchType::kMatchByFirstDaughter: {
      switch (fFormat) {
      case EFormat::kDst: MatchDstByDaugher(0); break;
      case EFormat::kMiniDst: MatchMiniDstByDaughter(0); break;
      }
   } break;
   case EMatchType::kMatchBySecondDaughet: {
      switch (fFormat) {
      case EFormat::kDst: MatchDstByDaugher(1); break;
      case EFormat::kMiniDst: MatchMiniDstByDaughter(1); break;
      }
   } break;
   case EMatchType::kMatchByBothDaughters: {
      switch (fFormat) {
      case EFormat::kDst: MatchDstByDaughers(); break;
      case EFormat::kMiniDst: MatchMiniDstByDaughters(); break;
      }
   }
   }
}

void MpdV0Matcher::MatchDstByDaugher(Int_t dau)
{
   TClonesArray *tracks = fMpdEvent->GetGlobalTracks();
   cout << "Matching " << fV0s->GetEntriesFast() << endl;

   for (int i = 0; i < fV0s->GetEntriesFast(); i++) {
      if (fLinks->GetLink(i) >= 0) continue;  // this track is already matched !
      if (fLinks->GetLink(i) == -2) continue; // this track is marked as bad !
      MpdV0Track *track = (MpdV0Track *)fV0s->UncheckedAt(i);
      if (!IsGoodV0(track->GetPdg())) continue; // this track should not be matched by this task
      Int_t dauIndex = -1;
      if (dau == 0) dauIndex = track->GetPositiveDaughterIndex();
      if (dau == 1) dauIndex = track->GetNegativeDaugterIndex();
      MpdTrack *recotrack = (MpdTrack *)tracks->UncheckedAt(dauIndex);
      Int_t     mcId      = recotrack->GetID();
      if (mcId < 0) { // no MC particle assigned to daughter, try another
         if (dau == 0) {
            dauIndex  = track->GetNegativeDaugterIndex();
            recotrack = (MpdTrack *)tracks->UncheckedAt(dauIndex);
            mcId      = recotrack->GetID();
         } else {
            dauIndex  = track->GetPositiveDaughterIndex();
            recotrack = (MpdTrack *)tracks->UncheckedAt(dauIndex);
            mcId      = recotrack->GetID();
         }
         if (mcId < 0) { // still no MC, mark track as bad
            fLinks->SetLink(i, -2);
            continue;
         }
      }
      MpdMCTrack *mctrack = (MpdMCTrack *)fMcTracks->UncheckedAt(mcId);
      if (mctrack->GetMotherId() < 0) { // primary or some "strange" particle cannot match
         fLinks->SetLink(i, -2);
         continue;
      }
      MpdMCTrack *motherTrack = (MpdMCTrack *)fMcTracks->UncheckedAt(mctrack->GetMotherId());
      if (motherTrack->GetPdgCode() == track->GetPdg()) {
         fLinks->SetLink(i, mctrack->GetMotherId()); // yeay this is our V0
      } else {
         fLinks->SetLink(i, -2); // cannot match this is not our V0
      }
   }
}

void MpdV0Matcher::MatchMiniDstByDaughter(Int_t dau)
{
   TClonesArray *tracks = fMpdEvent->GetGlobalTracks();
   cout << "Matching " << fV0s->GetEntriesFast() << endl;

   for (int i = 0; i < fV0s->GetEntriesFast(); i++) {
      if (fLinks->GetLink(i) >= 0) continue;  // this track is already matched !
      if (fLinks->GetLink(i) == -2) continue; // this track is marked as bad !
      MpdV0Track *track = (MpdV0Track *)fV0s->UncheckedAt(i);
      if (!IsGoodV0(track->GetPdg())) continue; // this track should not be matched by this task
      Int_t dauIndex = -1;
      if (dau == 0) dauIndex = track->GetPositiveDaughterIndex();
      if (dau == 1) dauIndex = track->GetNegativeDaugterIndex();
      MpdTrack *recotrack = (MpdTrack *)tracks->UncheckedAt(dauIndex);
      Int_t     mcId      = recotrack->GetID();
      if (mcId < 0) { // no MC particle assigned to daughter, try another
         if (dau == 0) {
            dauIndex  = track->GetNegativeDaugterIndex();
            recotrack = (MpdTrack *)tracks->UncheckedAt(dauIndex);
            mcId      = recotrack->GetID();
         } else {
            dauIndex  = track->GetPositiveDaughterIndex();
            recotrack = (MpdTrack *)tracks->UncheckedAt(dauIndex);
            mcId      = recotrack->GetID();
         }
         if (mcId < 0) { // still no MC, mark track as bad
            fLinks->SetLink(i, -2);
            continue;
         }
      }
      MpdMCTrack *mctrack = (MpdMCTrack *)fMcTracks->UncheckedAt(mcId);
      if (mctrack->GetMotherId() < 0) { // primary or some "strange" particle cannot match
         fLinks->SetLink(i, -2);
         continue;
      }
      MpdMCTrack *motherTrack = (MpdMCTrack *)fMcTracks->UncheckedAt(mctrack->GetMotherId());
      if (motherTrack->GetPdgCode() == track->GetPdg()) {
         fLinks->SetLink(i, mctrack->GetMotherId()); // yeay this is our V0
      } else {
         fLinks->SetLink(i, -2); // cannot match this is not our V0
      }
   }
}

void MpdV0Matcher::PrepareMomentumMatchingDst()
{
   for (int i = 0; i < fMcTracks->GetEntriesFast(); i++) {
      MpdMCTrack *track = (MpdMCTrack *)fMcTracks->UncheckedAt(i);
      if (IsGoodV0(track->GetPdgCode())) {
         fMcV0Indexes.push_back(i);
         fMcV0Momenta.push_back(TVector3(track->GetPx(), track->GetPy(), track->GetPz()));
         fMcV0sPid.push_back(track->GetPdgCode());
      }
   }
}

void MpdV0Matcher::PrepareMomentumMatchinMiniDst()
{
   for (int i = 0; i < fMiniTracks->GetEntriesFast(); i++) {
      MpdMiniMcTrack *track = (MpdMiniMcTrack *)fMiniTracks->UncheckedAt(i);
      if (IsGoodV0(track->pdgId())) {
         fMcV0Indexes.push_back(i);
         fMcV0Momenta.push_back(TVector3(track->px(), track->py(), track->pz()));
         fMcV0sPid.push_back(track->pdgId());
      }
   }
}

Bool_t MpdV0Matcher::IsGoodV0(Int_t pid) const
{
   for (unsigned int i = 0; i < fMcV0sPid[i]; i++) {
      if (pid == fMcV0sPid[i]) return kTRUE;
   }
   return kFALSE;
}

MpdV0Matcher::~MpdV0Matcher() {}

void MpdV0Matcher::MatchDstByDaughers()
{
   TClonesArray *tracks = fMpdEvent->GetGlobalTracks();
   cout << "Matching " << fV0s->GetEntriesFast() << endl;

   for (int i = 0; i < fV0s->GetEntriesFast(); i++) {
      if (fLinks->GetLink(i) >= 0) continue;  // this track is already matched !
      if (fLinks->GetLink(i) == -2) continue; // this track is marked as bad !
      MpdV0Track *track = (MpdV0Track *)fV0s->UncheckedAt(i);
      if (!IsGoodV0(track->GetPdg())) continue; // this track should not be matched by this task
      Int_t     dauIndex1  = track->GetPositiveDaughterIndex();
      Int_t     dauIndex2  = track->GetNegativeDaugterIndex();
      MpdTrack *recotrack1 = (MpdTrack *)tracks->UncheckedAt(dauIndex1);
      MpdTrack *recotrack2 = (MpdTrack *)tracks->UncheckedAt(dauIndex2);
      Int_t     mcId1      = recotrack1->GetID();
      Int_t     mcId2      = recotrack1->GetID();
      if (mcId1 < 0 || mcId2 < 0) {
         fLinks->SetLink(i, -2);
         continue;
      }
      MpdMCTrack *mctrack1 = (MpdMCTrack *)fMcTracks->UncheckedAt(mcId1);
      MpdMCTrack *mctrack2 = (MpdMCTrack *)fMcTracks->UncheckedAt(mcId2);
      if (mctrack1->GetMotherId() != mctrack2->GetMotherId()) { // different mothers
         fLinks->SetLink(i, -2);
         continue;
      }
      if (mctrack1->GetMotherId() < 0) { // primary or some "strange" particle cannot match
         fLinks->SetLink(i, -2);
         continue;
      }
      if (mctrack2->GetMotherId() < 0) { // primary or some "strange" particle cannot match
         fLinks->SetLink(i, -2);
         continue;
      }
      MpdMCTrack *motherTrack = (MpdMCTrack *)fMcTracks->UncheckedAt(mctrack1->GetMotherId());
      if (motherTrack->GetPdgCode() == track->GetPdg()) {
         fLinks->SetLink(i, mctrack1->GetMotherId()); // yeay this is our V0
      } else {
         fLinks->SetLink(i, -2); // cannot match this is not our V0
      }
   }
}

void MpdV0Matcher::MatchMiniDstByDaughters() {}
