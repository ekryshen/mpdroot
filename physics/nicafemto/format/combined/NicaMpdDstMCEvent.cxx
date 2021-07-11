/*
 * NicaMpdDstMCEvent.cxx
 *
 *  Created on: 1 kwi 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "NicaMpdDstMCEvent.h"

#include <TClonesArray.h>
#include <TObjArray.h>
#include <stddef.h>

#include "MpdTrack.h"
#include "NicaComplexTrack.h"
#include "NicaEvent.h"
#include "NicaMpdEvent.h"
#include "NicaMpdMcEvent.h"
#include "NicaTrack.h"

NicaMpdDstMCEvent::NicaMpdDstMCEvent() : NicaComplexEvent(new NicaMpdEvent(), new NicaMpdMcEvent()) {}

void NicaMpdDstMCEvent::OnlyPrimary() { ((NicaMpdEvent*) fRealEvent)->OnlyPrimary(); }

void NicaMpdDstMCEvent::OnlyGlobal() { ((NicaMpdEvent*) fRealEvent)->OnlyGlobal(); }

void NicaMpdDstMCEvent::Update() {
  fImgEvent->Update();
  fRealEvent->Update();
  NicaEvent::ShallowCopyEvent(fRealEvent);
  fTracks->Clear();
  fTotalTracksNo = fRealEvent->GetTotalTrackNo();
  fTracks->ExpandCreateFast(fTotalTracksNo);
  for (int i = 0; i < fTotalTracksNo; i++) {
    NicaComplexTrack* track = (NicaComplexTrack*) fTracks->UncheckedAt(i);
    track->NicaTrack::ResetTrack(i, this);
    track->SetRealTrack(fRealEvent->GetTrack(i));
    MpdTrack* mpd_track = (MpdTrack*) track->GetRealTrack()->GetTrackPointer();
    Int_t parent_id     = mpd_track->GetID();
    track->NicaTrack::CopyData(fRealEvent->GetTrack(i));
    track->SetMatchID(parent_id);
    if (parent_id >= 0) { track->SetImgTrack(fImgEvent->GetTrack(parent_id)); }
  }
}

NicaMpdDstMCEvent::~NicaMpdDstMCEvent() {}
