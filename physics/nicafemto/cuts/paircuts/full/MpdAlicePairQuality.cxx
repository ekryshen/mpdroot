/*
 * MpdAlicePairQuality.cxx
 *
 *  Created on: 26 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdAlicePairQuality.h"

#include "NicaMpdTrack.h"
#include "NicaTwoTrack.h"
namespace MpdPadsFormat {
MpdAlicePairQuality::MpdAlicePairQuality() : MpdTpcPadsPairCut(1) {
  SetUnitName("Q_{ALICE}[%]");
  SetMinMax(-1, 1);
}

Bool_t MpdAlicePairQuality::Pass(NicaTwoTrack *pair) {
  NicaMpdTrack *tr1 = (NicaMpdTrack *)pair->GetTrack1();
  NicaMpdTrack *tr2 = (NicaMpdTrack *)pair->GetTrack2();
  Double_t nhits = tr1->GetNHitsTpc() + tr2->GetNHitsTpc();
  ULong64_t shared_map1 = tr1->GetSharedHitMap();
  ULong64_t shared_map2 = tr2->GetSharedHitMap();
  ULong64_t sum = shared_map1 & shared_map2;
  ULong64_t hits_map1 = tr1->GetHitMap();
  ULong64_t hits_map2 = tr2->GetHitMap();
  Double_t q = 0;
  for (int i = 0; i < 53; i++) {
    Int_t hit1 = TESTBIT(hits_map1, i);
    Int_t hit2 = TESTBIT(hits_map2, i);
    Int_t sum_hit = hit1 + hit2;
    switch (sum_hit) {
      case 2:                   // both have hits
        if (TESTBIT(sum, i)) {  // and shared
          q++;
        } else {
          q--;
        }
        break;
      case 1: {  // one hit in padrow
        q++;
      } break;
      default:
        break;
    }
  }
  SetValue(q / nhits);
  return Validate();
}

MpdAlicePairQuality::~MpdAlicePairQuality() {}
}  // namespace MpdPadsFormat
