/*
 * NicaMiniDstFullV0Event.cxx
 *
 *  Created on: 26 lut 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMiniDstFullV0Event.h"
#include "MpdMiniTrack.h"
#include "NicaComplexTrack.h"
#include "NicaMpdMiniDstEventV0.h"
#include "NicaMpdMiniDstEventV0Interface.h"
#include "NicaMpdMiniDstMcEvent.h"

#include <iostream>

NicaMiniDstFullV0Event::NicaMiniDstFullV0Event() :
  NicaMpdMiniDstFullEvent(new NicaMpdMiniDstEventV0(), new NicaMpdMiniDstMcEvent()) {}

void NicaMiniDstFullV0Event::Update() {
  NicaMpdMiniDstEventV0Interface* interface = (NicaMpdMiniDstEventV0Interface*) fSource;
  fImgEvent->Update();
  fRealEvent->Update();
  NicaEvent::ShallowCopyEvent(fRealEvent);
  fTracks->Clear();
  fTotalTracksNo = fRealEvent->GetTotalTrackNo();
  fTracks->ExpandCreateFast(fTotalTracksNo);
  Int_t nV0s         = fRealEvent->GetTotalV0No();
  Int_t normalTracks = fTotalTracksNo - nV0s;
  for (int i = 0; i < fTotalTracksNo; i++) {
    NicaComplexTrack* track = (NicaComplexTrack*) fTracks->UncheckedAt(i);
    track->ResetTrack(i, this);
    track->SetRealTrack(fRealEvent->GetTrack(i));
    track->NicaTrack::CopyData(fRealEvent->GetTrack(i));
    if (track->GetRealTrack()->IsGoodV0()) {
      NicaTrack* reTrack = track->GetRealTrack();
      TClonesArray* v0   = interface->GetV0Links();
      Int_t v0Index      = reTrack->GetHiddenInfoIndex();
      NicaLink* link     = (NicaLink*) v0->UncheckedAt(v0Index);
      Int_t index        = link->GetLink(0);
      track->SetMatchID(index);
      if (index >= 0) track->SetImgTrack(fImgEvent->GetTrack(index));
    } else {
      MpdMiniTrack* mpd_track = (MpdMiniTrack*) track->GetRealTrack()->GetTrackPointer();
      Int_t parent_id         = mpd_track->mcTrackIndex();
      track->SetMatchID(parent_id);
      if (parent_id >= 0) track->SetImgTrack(fImgEvent->GetTrack(parent_id));
    }
  }
}

NicaMiniDstFullV0Event::~NicaMiniDstFullV0Event() {}
