/*
 * NicaMpdDstEventV0.cxx
 *
 *  Created on: 25 lut 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdMiniDstEventV0.h"
#include "NicaCout.h"
#include "NicaMpdMiniDstEventV0Interface.h"
#include "NicaMpdMiniDstTrack.h"
#include "NicaV0CandidateHelix.h"

#include <FairRootManager.h>


NicaMpdMiniDstEventV0::NicaMpdMiniDstEventV0(TString trackname) : NicaMpdMiniDstEvent(trackname) {}

NicaMpdMiniDstEventV0::NicaMpdMiniDstEventV0(NicaMpdMiniDstEvent::eMode Mode) : NicaMpdMiniDstEventV0("NicaMpdMiniDstTrack") {
  fMode = Mode;
}

void NicaMpdMiniDstEventV0::CreateSource() { fSource = new NicaMpdMiniDstEventV0Interface(); }

void NicaMpdMiniDstEventV0::ShallowCopyEvent(NicaEvent* event) { NicaMpdMiniDstEvent::ShallowCopyEvent(event); }

void NicaMpdMiniDstEventV0::Update() {
  NicaMpdMiniDstEvent::Update();
  fV0sHiddenInfo->Clear();
  NicaMpdMiniDstEventV0Interface* interface = (NicaMpdMiniDstEventV0Interface*) fSource;
  Int_t normalTracks                        = fTotalTracksNo;
  Int_t nV0s                                = interface->fV0Tracks->GetEntriesFast();
  fTotalTracksNo                            = normalTracks + nV0s;
  fTracks->ExpandCreateFast(fTotalTracksNo);
  Int_t posV0 = 0;
  for (int i = normalTracks; i < fTotalTracksNo; i++) {
    NicaMpdMiniDstTrack* track = (NicaMpdMiniDstTrack*) fTracks->UncheckedAt(i);
    track->ResetTrack(i, this);
    NicaV0Track* v0 = (NicaV0Track*) interface->fV0Tracks->UncheckedAt(posV0);
    track->SetV0(kTRUE, kTRUE);
    track->GetV0Info()->CopyData(v0);
    track->SetPrimary();
    posV0++;
  }
}

Bool_t NicaMpdMiniDstEventV0::ExistInTree() const {
  if (NicaMpdMiniDstEvent::ExistInTree() == kFALSE) return kFALSE;
  FairRootManager* manager = FairRootManager::Instance();
  manager->Print();
  if (manager->GetObject("NicaV0Tracks") != nullptr) { return kTRUE; }
  return kFALSE;
}

TString NicaMpdMiniDstEventV0::GetFormatName() const {
  switch (fMode) {
    case kGlobalTrack: return "NicaMpdMiniDstEventGlobalV0"; break;
    case kPrimaryTrack: return "NicaMpdMiniDstEventPrimaryV0"; break;
  }
  return "";
}

NicaMpdMiniDstEventV0::~NicaMpdMiniDstEventV0() {}
