/*
 * NicaTrackIsGlobalCut.cxx
 *
 *  Created on: 22 gru 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdTrackIsGlobalCut.h"

#include "NicaMpdMiniDstTrack.h"

MpdTrackIsGlobalCut::MpdTrackIsGlobalCut() : NicaTrackCut(1) { SetUnitName("IsGlobal [AU]"); }

Bool_t MpdTrackIsGlobalCut::Pass(NicaTrack* track) {
  NicaMpdMiniDstTrack* tr = static_cast<NicaMpdMiniDstTrack*>(track);
  if (tr->IsGlobal()) {
    SetValue(1);
  } else {
    SetValue(0);
  }
  return Validate();
}

MpdTrackIsGlobalCut::~MpdTrackIsGlobalCut() {}

Bool_t MpdTrackIsGlobalCut::Init(Int_t format_id) {
  return FormatInhertis("NicaMpdMiniDstEvent", format_id, ENicaFormatDepth::kNonBuffered);
}
