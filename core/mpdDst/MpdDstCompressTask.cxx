/*
 * MpdDstWriteTask.cxx
 *
 *  Created on: 23 lut 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdDstCompressTask.h"

#include "FairMCEventHeader.h"
#include "FairRunAna.h"
MpdDstCompressTask::MpdDstCompressTask() : MpdDstCompressTask("dstwirite", 1) {}

InitStatus MpdDstCompressTask::CheckBranches()
{
   FairRootManager *mngr = FairRootManager::Instance();
   fMpdEvent             = (MpdEvent *)mngr->GetObject("MPDEvent.");
   if (fMpdEvent == nullptr) return kFATAL;
   mngr->Register("MPDEvent.", "MPD", fMpdEvent, kTRUE);
   if (fUseMC) {
      fMCTracks = (TClonesArray *)mngr->GetObject("MCTrack");
      if (fMCTracks == nullptr) {
         LOG(WARNING) << "MC tracks requested but not found!";
         fUseMC = kFALSE;
      } else {
         mngr->Register("MCTrack", "MC", fMCTracks, kTRUE);
      }
   }
   if (fUseFreezouts) {
      fFreezouts = (TClonesArray *)mngr->GetObject("Freezouts.");
      if (fFreezouts == nullptr) {
         LOG(WARNING) << "Freeouts tracks requested but not found!";
         fUseFreezouts = kFALSE;
      } else {
         mngr->Register("Freezouts.", "Freezouts", fFreezouts, kTRUE);
      }
   }
   if (fUseTpcKalmans) {
      fTpcKalmans = (TClonesArray *)mngr->GetObject("TpcKalmanTrack");
      if (fTpcKalmans == nullptr) {
         LOG(WARNING) << "Kalman TPC tracks requested but not found!";
         fUseTpcKalmans = kFALSE;
      } else {
         mngr->Register("TpcKalmanTrack", "TPC", fTpcKalmans, kTRUE);
      }
   }
   if (fUseTpcHits) {
      fTpcHits = (TClonesArray *)mngr->GetObject("TpcHit");
      if (fTpcHits == nullptr) {
         LOG(WARNING) << "TPC hits requested but not found!";
         fUseTpcHits = kFALSE;
      } else {
         mngr->Register("TpcHi.", "TPC", fTpcHits, kTRUE);
      }
   }
   if (!fUseHeader) {
      /*	 TODO when FairROOT will be ugpraded
         FairRunAna::Instance()->SetEventHeaderPersistence(kFALSE);
         */
      TNamed *EventHeader = (TNamed *)mngr->GetObject("MCEventHeader.");
      if (EventHeader) {
         mngr->Register("MCEventHeader.", "MC", EventHeader, kTRUE);
      }
      if (EventHeader == nullptr) {
         EventHeader = (FairMCEventHeader *)mngr->GetObject("EventHeader.");
         mngr->Register("EventHeader.", "MC", EventHeader, kTRUE);
      }
   }
   return kSUCCESS;
}

InitStatus MpdDstCompressTask::Init()
{
   if (CheckBranches() == kFATAL) {
      LOG(FATAL) << "End of macro MPDEvent not found";
      return kFATAL;
   }
   fMCMapSize  = 1000;
   fMCIndexMap = new Int_t[fMCMapSize];
   return kSUCCESS;
}

MpdDstCompressTask::MpdDstCompressTask(const char *name, Int_t Verbose)
   : FairTask(name, Verbose), fUseMC(kFALSE), fUseFreezouts(kFALSE), fUseTpcKalmans(kFALSE), fUseTpcHits(kFALSE),
     fUseHeader(kFALSE), fMCCompression(kFALSE), fMpdEvent(nullptr), fFreezouts(nullptr), fMCTracks(nullptr),
     fTpcKalmans(nullptr), fTpcHits(nullptr), fMCMapSize(0), fMCIndexMap(nullptr)
{
}

MpdDstCompressTask::~MpdDstCompressTask()
{
   if (fMCIndexMap) delete[] fMCIndexMap;
}

void MpdDstCompressTask::Exec(Option_t *optio)
{
   if (fMCCompression) {
      TClonesArray *glob_tracks = fMpdEvent->GetGlobalTracks();
      TClonesArray *prim_tracks = fMpdEvent->GetPrimaryTracks();
      if (fMCTracks->GetEntriesFast() > fMCMapSize) {
         delete[] fMCIndexMap;
         fMCMapSize  = fMCTracks->GetSize() * 2;
         fMCIndexMap = new Int_t[fMCMapSize];
      }
      for (int iMCTrack = 0; iMCTrack < fMCTracks->GetEntriesFast(); iMCTrack++) {
         fMCIndexMap[iMCTrack] = -1;
      }
      for (int iTrack = 0; iTrack < glob_tracks->GetEntriesFast(); iTrack++) {
         MpdTrack *track   = (MpdTrack *)glob_tracks->UncheckedAt(iTrack);
         Int_t     matched = track->GetID();
         if (matched >= 0) fMCIndexMap[matched] = 0;
      }
      for (int iTrack = 0; iTrack < prim_tracks->GetEntriesFast(); iTrack++) {
         MpdTrack *track   = (MpdTrack *)prim_tracks->UncheckedAt(iTrack);
         Int_t     matched = track->GetID();
         if (matched >= 0) fMCIndexMap[matched] = 0;
      }
      Int_t index = 0;
      for (int iMCTrack = 0; iMCTrack < fMCTracks->GetEntriesFast(); iMCTrack++) {
         if (fMCIndexMap[iMCTrack] > -1) { // this track is matched
            fMCIndexMap[iMCTrack] = index++;
         } else {
            fMCTracks->RemoveAt(iMCTrack);
         }
      }
      fMCTracks->Compress();
      // set new matching ID's
      for (int iTrack = 0; iTrack < glob_tracks->GetEntriesFast(); iTrack++) {
         MpdTrack *track   = (MpdTrack *)glob_tracks->UncheckedAt(iTrack);
         Int_t     matched = track->GetID();
         if (matched > -1) track->SetID(fMCIndexMap[matched]);
      }
      for (int iTrack = 0; iTrack < prim_tracks->GetEntriesFast(); iTrack++) {
         MpdTrack *track   = (MpdTrack *)prim_tracks->UncheckedAt(iTrack);
         Int_t     matched = track->GetID();
         if (matched > -1) track->SetID(fMCIndexMap[matched]);
      }
   }
}
