/*
 * MpdPairSharedHitsCut.cxx
 *
 *  Created on: 26 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairSharedHitsCut.h"

#include "NicaMpdTrack.h"
#include "NicaTwoTrack.h"

namespace MpdPadsFormat {
MpdPairSharedHitsCut::MpdPairSharedHitsCut() : MpdTpcPadsPairCut(1) {
  SetUnitName("SharedHits[%]");
}

Bool_t MpdPairSharedHitsCut::Pass(NicaTwoTrack* pair) {
  ULong64_t tr1 = ((NicaMpdTrack*)pair->GetTrack1())->GetSharedHitMap();
  ULong64_t tr2 = ((NicaMpdTrack*)pair->GetTrack2())->GetSharedHitMap();
  ULong64_t sum = tr1 & tr2;
  Float_t shared = 0;
  for (int i = 0; i < 53; i++) {
    if (TESTBIT(sum, i)) shared++;
  }
  SetValue(shared / 53.0);
  return Validate();
}

MpdPairSharedHitsCut::~MpdPairSharedHitsCut() {}
}  // namespace MpdPadsFormat
