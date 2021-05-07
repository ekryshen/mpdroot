/*
 * MpdSharedPadsCut.cxx
 *
 *  Created on: 30 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairSharedPadsCut.h"

#include "NicaMpdTrackTpcPads.h"
#include "NicaTpcTrack.h"
#include "NicaTwoTrack.h"

namespace MpdPadsFormat {
MpdPairSharedPadsCut::MpdPairSharedPadsCut(Bool_t shift)
    : MpdModularTpcPairCut(2), fShift(shift) {
  SetUnitName("Overlap Pads [%]", OverlappedPads());
  SetUnitName("SharedPads [%]", SharedPads());
  SetMinMax(0, 100, OverlappedPads());
  SetMinMax(0, 100, SharedPads());
}

Bool_t MpdPairSharedPadsCut::Pass(NicaTwoTrack *pair) {
  NicaMpdTrackTpcPads *track1 = (NicaMpdTrackTpcPads *)pair->GetTrack1();
  NicaMpdTrackTpcPads *track2 = (NicaMpdTrackTpcPads *)pair->GetTrack2();
  track1->CalculatePads(fShift);
  track2->CalculatePads(fShift);
  const Double_t pads1 = track1->GetTpcPadsInfo()->GetPadsNo();
  const Double_t pads2 = track2->GetTpcPadsInfo()->GetPadsNo();
  const Int_t start1_pad = track1->GetTpcPadsInfo()->GetFirstGoodPad();
  const Int_t start2_pad = track2->GetTpcPadsInfo()->GetFirstGoodPad();
  const Int_t end1_pad = track1->GetTpcPadsInfo()->GetLastGoodPad() + 1;
  const Int_t end2_pad = track2->GetTpcPadsInfo()->GetLastGoodPad() + 1;
  Double_t maxPads =
      TMath::Max(end1_pad, end2_pad) - TMath::Min(start1_pad, start2_pad);
  const Int_t first_common_pad = TMath::Max(start1_pad, start2_pad);
  const Int_t last_common_pad = TMath::Min(end1_pad, end2_pad);
  Double_t overlapPads = 0;
  Double_t sharedPads = 0;
  ULong64_t map1 = track1->GetHitMap();
  ULong64_t map2 = track2->GetHitMap();
  for (int iLay = first_common_pad; iLay < last_common_pad; iLay++) {
    Int_t pad1_id = track1->GetTpcPadsInfo()->GetPadID(iLay);
    Int_t pad2_id = track2->GetTpcPadsInfo()->GetPadID(iLay);
    if (pad1_id == pad2_id) {
      overlapPads++;
      if (TESTBIT(map1, iLay) &&
          TESTBIT(map2, iLay)) {  // have also hit in this rows
        sharedPads++;
      }
    }
  }
  SetValue(100.0 * overlapPads / maxPads, OverlappedPads());
  SetValue(100.0 * sharedPads / maxPads, SharedPads());
  return Validate();
}

MpdPairSharedPadsCut::~MpdPairSharedPadsCut() {
  // TODO Auto-generated destructor stub
}
}  // namespace MpdPadsFormat
