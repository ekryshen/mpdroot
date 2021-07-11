/*
 * MpdRejectSPlittedPairsCut.cxx
 *
 *  Created on: 28 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdSplittedPairsCut.h"
#include "NicaTwoTrack.h"

MpdSplittedPairsCut::MpdSplittedPairsCut() : NicaTwoTrackCut(1) {
  fReject = kTRUE;
}

Bool_t MpdSplittedPairsCut::Pass(NicaTwoTrack *pair) {
  NicaComplexTrack *track1 = (NicaComplexTrack *)pair->GetTrack1();
  NicaComplexTrack *track2 = (NicaComplexTrack *)pair->GetTrack2();
  if (track1->GetMatchID() == track2->GetMatchID()) {
    return ForcedUpdate(!fReject);
  } else {
    return ForcedUpdate(fReject);
  }
}

Bool_t MpdSplittedPairsCut::Init(Int_t task_id) {
  return FormatInhertis("NicaComplexEvent", task_id,
                        ENicaFormatDepth::kBuffered);
}

MpdSplittedPairsCut::~MpdSplittedPairsCut() {}
