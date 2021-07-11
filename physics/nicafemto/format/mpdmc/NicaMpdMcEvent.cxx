/*
 * NicaEventFair.cxx
 *
 *  Created on: 05-07-2014
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdMcEvent.h"

#include "FairRootManager.h"
#include "NicaMpdMcEventInterface.h"
#include <iostream>

NicaMpdMcEvent::NicaMpdMcEvent() : NicaMCEvent("NicaMpdMcTrack") {}

NicaMpdMcEvent::NicaMpdMcEvent(TString trackname) : NicaMCEvent(trackname) {}

NicaMpdMcEvent::NicaMpdMcEvent(const NicaMpdMcEvent& other) : NicaMCEvent(other) {}

void NicaMpdMcEvent::Update() {
  FairMCEventHeader* event = (FairMCEventHeader*) ((NicaMpdMcEventInterface*) fSource)->fEvent;
  TClonesArray* tracks     = (TClonesArray*) ((NicaMpdMcEventInterface*) fSource)->fMcracks->GetArray();
  fB                       = event->GetB();
  fVertex->SetXYZT(event->GetX(), event->GetY(), event->GetZ(), event->GetT());
  fTracks->Clear();
  fTotalTracksNo = tracks->GetEntriesFast();
  fTracks->ExpandCreateFast(fTotalTracksNo);
  for (int i = 0; i < tracks->GetEntriesFast(); i++) {
    MpdMCTrack* track = (MpdMCTrack*) tracks->UncheckedAt(i);
    NicaMCTrack* mc   = (NicaMCTrack*) fTracks->UncheckedAt(i);
    mc->ResetTrack(i, this);
    mc->SetMomentum(track->GetPx(), track->GetPy(), track->GetPz(), track->GetEnergy());
    mc->SetCharge(CalculateCharge(track->GetPdgCode()));
    mc->SetPdg(track->GetPdgCode());
    mc->SetStartPosition(track->GetStartX(), track->GetStartY(), track->GetStartZ(), track->GetStartT());
    mc->SetID(i);
    if (track->GetMotherId() == -1) {
      mc->SetPrimary();
    } else {
      if (track->GetMotherId() == -2)
        mc->SetSecondary(kFALSE);
      else {
        mc->SetMotherIndex(track->GetMotherId());
      }
    }
  }
}

void NicaMpdMcEvent::Clear(Option_t* opt) { NicaMCEvent::Clear(opt); }

void NicaMpdMcEvent::Print() {}

void NicaMpdMcEvent::CreateSource() {
  std::cout << "Create source" << std::endl;
  fSource = new NicaMpdMcEventInterface();
}

NicaMpdMcEvent::~NicaMpdMcEvent() {}

TString NicaMpdMcEvent::GetFormatName() const { return "NicaMpdMcEvent"; }

Bool_t NicaMpdMcEvent::ExistInTree() const {
  FairRootManager* manager = FairRootManager::Instance();
  Int_t header             = manager->CheckBranch("MCEventHeader.") + manager->CheckBranch("EventHeader.");
  if (header > 1) header = 1;
  Int_t tracks = manager->CheckBranch("MCTrack");
  if ((header + tracks) == 2) { return kTRUE; }
  return kFALSE;
}
