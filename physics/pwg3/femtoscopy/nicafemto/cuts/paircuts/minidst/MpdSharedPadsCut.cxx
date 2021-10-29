/*
 * MpdSharedPadsCut.cxx
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdSharedPadsCut.h"
#include "NicaMpdHbtTrack.h"
#include "NicaTrackTpcPads.h"
#include "NicaTwoTrack.h"

namespace MpdHbtDst {
MpdSharedPadsCut::MpdSharedPadsCut() : MpdFemtoPairCut(2) {
  SetUnitName("SharedPads [N]", SharedPadsNo());
  SetUnitName("SharedPadsFrac [%]", SharedPadsFraction());
}

Bool_t MpdSharedPadsCut::Pass(NicaTwoTrack *pair) {
  NicaMpdHbtTrack *track1 = (NicaMpdHbtTrack *)pair->GetTrack1();
  NicaMpdHbtTrack *track2 = (NicaMpdHbtTrack *)pair->GetTrack2();
  NicaTrackTpcPads *pads1 = track1->GetPadsInfo();
  NicaTrackTpcPads *pads2 = track2->GetPadsInfo();
  Int_t start_pad =
      TMath::Max(pads1->GetFirstGoodPad(), pads2->GetFirstGoodPad());
  Int_t end_pad = TMath::Min(pads1->GetLastGoodPad(), pads2->GetLastGoodPad());
  Double_t padsTotal = end_pad - start_pad + 1;
  Double_t shared_pads = 0;
  for (int i = start_pad; i <= end_pad; i++) {
    if (pads1->GetPadID(i) == pads2->GetPadID(i)) shared_pads++;
  }
  if (padsTotal == 0) {
    SetValue(100, SharedPadsFraction());
  } else {
    SetValue(100.0 * shared_pads / padsTotal, SharedPadsFraction());
  }
  SetValue(shared_pads, SharedPadsNo());
  return Validate();
}

MpdSharedPadsCut::~MpdSharedPadsCut() {}

}  // namespace MpdHbtDst
