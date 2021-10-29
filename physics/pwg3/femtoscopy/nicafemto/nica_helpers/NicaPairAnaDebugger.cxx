/*
 * NicaPairAnaDebugger.cxx
 *
 *  Created on: 21 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "NicaPairAnaDebugger.h"

#include <iostream>

#include "NicaMemoryMapManager.h"
#include "NicaCutCollection.h"
#include "NicaCutContainer.h"
#include "NicaComplexTrack.h"
#include "NicaEvent.h"
#include "NicaTrack.h"

NicaPairAnaDebugger::NicaPairAnaDebugger() {}

void NicaPairAnaDebugger::ProcessEvent() {
  fMemoryMap->PrepareMaps(fCurrentEventCollectionID);
  NicaCutCollection *cont =
      fCutContainer->GetEventCollection(fCurrentEventCollectionID);
  for (int i = 0; i < fMemoryMap->GetTemporaryTotalTracksNo(); i++) {
    fCurrentTrack = fCurrentEvent->GetTrack(i);
    for (int j = 0; j < cont->GetNextNo(); j++) {
      fCurrentTrackCollectionID = cont->GetNextAddr(j);
      if (fCutContainer->PassTrack(fCurrentTrack, fCurrentTrackCollectionID)) {
        fMemoryMap->AddTrackToMapTrack(
            fCurrentEventCollectionID, fCurrentTrackCollectionID,
            i);  // load track into memory map - may be usefull at finish event
        ProcessTrack();
        NicaComplexTrack *t = (NicaComplexTrack *)fCurrentTrack;
        //			std::cout<<"TR"<<t->GetRealTrack()->GetPx()<<"
        //"<<t->GetImgTrack()->GetPx()<<std::endl;
        //			std::cout<<"TR"<<t->GetRealTrack()<<"
        //"<<t->GetImgTrack()<<std::endl;
      }
    }
  }
  fMemoryMap->BufferEvent(fCurrentEventCollectionID);
  NicaTrack *track, *track2;
  if (fMemoryMap->GetTracksNo(fCurrentEventCollectionID,
                              fCurrentTrackCollectionID) == 0) {
    if (fBackgroundMode == kMixedPairsID) {
      fMemoryMap->RejectLastEvent(fCurrentEventCollectionID);
    }
    return;
  }
  Int_t nTrack = fMemoryMap->GetTracksNo(fCurrentEventCollectionID,
                                         fCurrentTrackCollectionID);
  for (int i = 0; i < nTrack; i++) {
    track = fMemoryMap->GetTrack(fCurrentEventCollectionID,
                                 fCurrentTrackCollectionID, i);
    NicaComplexTrack *cpz = (NicaComplexTrack *)track;
    std::cout << " NMX " << std::endl;
    std::cout << "\t" << cpz->GetRealTrack()->GetPx() << " "
              << cpz->GetImgTrack()->GetPx() << std::endl;
    std::cout << "\t" << cpz->GetRealTrack() << " " << cpz->GetImgTrack()
              << std::endl;
  }
  // return;
  NicaCutCollection *track_cuts =
      fCutContainer->GetTrackCollection(fCurrentTrackCollectionID);
  Int_t tt_cut_no = track_cuts->GetNextNoBackround();
  Int_t nTrackA = fMemoryMap->GetTracksNo(fCurrentEventCollectionID,
                                          fCurrentTrackCollectionID);
  NicaComplexTrack *tr1, *tr2;
  if (fMemoryMap->IsReadyToMixing(fCurrentEventCollectionID)) {
    for (int l = 0; l < fMixSize; l++) {
      if (l == fMemoryMap->GetCounter(fCurrentEventCollectionID)) continue;
      Int_t nTrackB = fMemoryMap->GetTracksNo(fCurrentEventCollectionID,
                                              fCurrentTrackCollectionID, l);
      for (int i = 0; i < nTrackA; i++) {
        tr1 = (NicaComplexTrack *)fMemoryMap->GetTrack(
            fCurrentEventCollectionID, fCurrentTrackCollectionID, i);
        for (int j = 0; j < nTrackB; j++) {
          tr2 = (NicaComplexTrack *)fMemoryMap->GetTrack(
              fCurrentEventCollectionID, fCurrentTrackCollectionID, l, j);
          tr1->ClassName();
          tr2->ClassName();
        }
      }
    }
  }
}

NicaPairAnaDebugger::~NicaPairAnaDebugger() {}
