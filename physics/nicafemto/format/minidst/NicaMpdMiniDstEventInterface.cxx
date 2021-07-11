/*
 * NicaMpdMiniDstEventInterface.cxx
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "NicaMpdMiniDstEventInterface.h"

#include <FairRootManager.h>
#include <TClonesArray.h>
#include <TObject.h>

#include "MpdMiniBTofPidTraits.h"
#include "MpdMiniTrack.h"
#include "NicaMpdMiniDstTrack.h"
#include "NicaTrackInterface.h"
#include "NicaV0CandidateHelix.h"

NicaMpdMiniDstEventInterface::NicaMpdMiniDstEventInterface() : fEvent(nullptr), fTracks(nullptr) {
  fEvent = new NicaTrackClones("MpdMiniEvent", "Event", "event");
  fEvent->GetArray()->ExpandCreateFast(1);
  fTracks  = new NicaTrackClones("MpdMiniTrack", "Track", "track");
  fTofInfo = new NicaTrackClones("MpdMiniBTofPidTraits", "BTofPidTraits", "tof");
  fEmcInfo = nullptr;
}

void NicaMpdMiniDstEventInterface::CopyData(NicaEventInterface* s) {
  NicaMpdMiniDstEventInterface* interface = static_cast<NicaMpdMiniDstEventInterface*>(s);
  fEvent->CopyFrom<MpdMiniEvent>(interface->fEvent->GetArray());
  fTracks->CopyFrom<MpdMiniTrack>(interface->fTracks->GetArray());
  fTofInfo->CopyFrom<MpdMiniBTofPidTraits>(interface->fTofInfo->GetArray());
}

void NicaMpdMiniDstEventInterface::Compress(Int_t* map, Int_t map_size) {
  Int_t track_pos = 0;
  for (int i = 0; i < map_size; i++) {
    Int_t good_track = map[i];
    for (int j = track_pos; j < good_track; j++) {
      fTracks->GetArray()->RemoveAt(j);
      track_pos = good_track + 1;
    }
  }
  fTracks->GetArray()->Compress();
}

void NicaMpdMiniDstEventInterface::CopyAndCompress(NicaEventInterface* s, Int_t* map, Int_t map_size) {
  NicaMpdMiniDstEventInterface* interface = static_cast<NicaMpdMiniDstEventInterface*>(s);
  fTracks->Clear();
  fTracks->CopyCompress<MpdMiniTrack>(interface->fTracks->GetArray(), map, map_size);
  fTofInfo->CopyFrom<MpdMiniBTofPidTraits>(interface->fTofInfo->GetArray());
}

void NicaMpdMiniDstEventInterface::FillTrackInterface(NicaTrackInterface* track, Int_t index) {
  track->SetRawTrack(fTracks->UncheckedAt(index));
}

Int_t NicaMpdMiniDstEventInterface::GetTotalTrackNo() const { return fTracks->GetEntriesFast(); }

NicaTrackInterface* NicaMpdMiniDstEventInterface::GetTrackInterface() const { return new NicaMpdMiniDstTrackInterface(); }

TObject* NicaMpdMiniDstEventInterface::GetRawTrackPointer(Int_t index) const { return fTracks->UncheckedAt(index); }

void NicaMpdMiniDstEventInterface::SetRunInfoId(Int_t i) { GetMiniEvent()->setRunId(i); }

Int_t NicaMpdMiniDstEventInterface::GetRunInfoId() const { return static_cast<MpdMiniEvent*>(fEvent->UncheckedAt(0))->runId(); }

void NicaMpdMiniDstEventInterface::SetMagneticField(TVector3 mag) { GetMiniEvent()->setMagneticField(mag.Z()); }

TVector3 NicaMpdMiniDstEventInterface::GetMagneticField() const {
  return TVector3(0, 0, static_cast<MpdMiniEvent*>(fEvent->UncheckedAt(0))->bField());
}

TLorentzVector NicaMpdMiniDstEventInterface::GetVertexError() const {
  TLorentzVector vec;
  vec.SetVect(static_cast<MpdMiniEvent*>(fEvent->UncheckedAt(0))->primaryVertexError());
  return vec;
}

TLorentzVector NicaMpdMiniDstEventInterface::GetVertex() const {
  TLorentzVector vec;
  vec.SetVect(static_cast<MpdMiniEvent*>(fEvent->UncheckedAt(0))->primaryVertex());
  return vec;
}

void NicaMpdMiniDstEventInterface::ConnectToTree() {
  FairRootManager* manager = FairRootManager::Instance();
  if (CanDeleteEvent()) {
    fEvent->DeleteClones();
    fTracks->DeleteClones();
    fTofInfo->DeleteClones();
  }
  fEvent->GetFromTree();
  fTracks->GetFromTree();
  fTofInfo->GetFromTree();
}

void NicaMpdMiniDstEventInterface::Register(Bool_t write) {
  fTracks->Register(write);
  fTofInfo->Register(write);
  fEvent->Register(write);
}

NicaMpdMiniDstEventInterface::~NicaMpdMiniDstEventInterface() {
  if (CanDeleteEvent()) {
    delete fEvent;
    delete fTracks;
    delete fTofInfo;
  }
}
