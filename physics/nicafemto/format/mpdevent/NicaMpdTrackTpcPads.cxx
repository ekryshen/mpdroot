/*
 * NicaMpdTrackTpcPads.cxx
 *
 *  Created on: 2 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdTrackTpcPads.h"
#include "NicaEvent.h"
#include "TLorentzVector.h"

NicaMpdTrackTpcPads::NicaMpdTrackTpcPads() : fPads(new NicaTrackTpcPads()) {}

void NicaMpdTrackTpcPads::Update(MpdTrack* track) {
  NicaMpdTrack::Update(track);
  fPads->Reset();
}

NicaMpdTrackTpcPads::~NicaMpdTrackTpcPads() { delete fPads; }

NicaMpdTrackTpcPads::NicaMpdTrackTpcPads(const NicaMpdTrackTpcPads& other)
    : NicaMpdTrack(other), fPads(new NicaTrackTpcPads(*other.fPads)) {}

NicaMpdTrackTpcPads& NicaMpdTrackTpcPads::operator=(
    const NicaMpdTrackTpcPads& other) {
  if (this != &other) {
    NicaMpdTrack::operator=(other);
    *fPads = *other.fPads;
  }
  return *this;
}

void NicaMpdTrackTpcPads::CopyData(NicaTrack* other) {
  NicaMpdTrack::CopyData(other);
  NicaMpdTrackTpcPads* track = (NicaMpdTrackTpcPads*)other;
  *fPads = *track->fPads;
  fPads->Reset();
}

void NicaMpdTrackTpcPads::CalculatePads(Bool_t shift) {
  if (shift) {
    fPads->Calculate(GetHelix(), GetEvent()->GetVertex());
  } else {
    fPads->Calculate(GetHelix(), nullptr);
  }
}

Bool_t NicaMpdTrackTpcPads::PadsCalculated() const {
  return fPads->AreCalculated();
  return kTRUE;
}
