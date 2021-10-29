/*
 * NicaMpdMiniDstEventV0Interface.cxx
 *
 *  Created on: 25 lut 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdMiniDstEventV0Interface.h"
#include "NicaCout.h"
#include "NicaLink.h"
#include "NicaMpdMiniDstTrack.h"
#include "NicaV0CandidateHelix.h"
#include <FairLogger.h>


NicaMpdMiniDstEventV0Interface::NicaMpdMiniDstEventV0Interface() : NicaMpdMiniDstEventInterface() {
  fV0Tracks = new NicaTrackClones("NicaV0Track", "NicaV0Tracks", "NicaV0");
  fV0Links  = new NicaTrackClones("NicaLink", "NicaV0Matches", "NicaV0");
}

void NicaMpdMiniDstEventV0Interface::ConnectToTree() {
  NicaMpdMiniDstEventInterface::ConnectToTree();
  if (CanDeleteEvent()) {
    fV0Tracks->DeleteClones();
    fV0Links->DeleteClones();
  }
  fV0Tracks->GetFromTree();
  fV0Links->GetFromTree();
}

void NicaMpdMiniDstEventV0Interface::Register(Bool_t write) {
  NicaMpdMiniDstEventInterface::Register(write);
  fV0Tracks->Register(write);
  fV0Links->Register(write);
}

void NicaMpdMiniDstEventV0Interface::CopyData(NicaEventInterface* s) { LOG(ERROR) << "Function not supported"; }

void NicaMpdMiniDstEventV0Interface::Compress(Int_t* map, Int_t map_size) { LOG(ERROR) << "Function not supproted !"; }

void NicaMpdMiniDstEventV0Interface::CopyAndCompress(NicaEventInterface* s, Int_t* map, Int_t map_size) {
  NicaMpdMiniDstEventInterface::CopyAndCompress(s, map, map_size);
  NicaMpdMiniDstEventV0Interface* interface = (NicaMpdMiniDstEventV0Interface*) s;
  fV0Tracks->CopyFrom<NicaV0CandidateHelix>(interface->fV0Tracks->GetArray());
  fV0Links->CopyFrom<NicaLink>(interface->fV0Links->GetArray());
}

NicaTrackInterface* NicaMpdMiniDstEventV0Interface::GetTrackInterface() const { return new NicaMpdMiniDstTrackInterface(); }

Int_t NicaMpdMiniDstEventV0Interface::GetTotalTrackNo() const {
  return NicaMpdMiniDstEventInterface::GetTotalTrackNo() + fV0Tracks->GetEntriesFast();
}

NicaMpdMiniDstEventV0Interface::~NicaMpdMiniDstEventV0Interface() {
  if (CanDeleteEvent()) {
    delete fV0Links;
    delete fV0Tracks;
  }
}
