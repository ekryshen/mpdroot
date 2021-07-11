/*
 * MpdStarQualityPairCut.cxx
 *
 *  Created on: 30 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdStarPairQualityCut.h"

#include "NicaDataFormat.h"
#include "NicaMpdTrack.h"
#include "NicaTwoTrack.h"

namespace MpdPadsFormat {
  MpdStarPairQualityCut::MpdStarPairQualityCut() : NicaTwoTrackCut(1) {
    SetUnitName("Q_{STAR} [AU]");
    SetMinMax(-0.5, 1);
  }

  Bool_t MpdStarPairQualityCut::Pass(NicaTwoTrack* pair) {
    NicaMpdTrack* track1 = (NicaMpdTrack*) pair->GetTrack1();
    NicaMpdTrack* track2 = (NicaMpdTrack*) pair->GetTrack2();
    Double_t hits        = ((NicaTpcTrack*) track1->GetDetTrack(NicaDetectorID::kTPC))->GetNHits()
                    + ((NicaTpcTrack*) track2->GetDetTrack(NicaDetectorID::kTPC))->GetNHits();
    ULong64_t map1 = track1->GetHitMap();
    ULong64_t map2 = track2->GetHitMap();
    Double_t q     = 0;
    for (int iLay = 0; iLay < 53; iLay++) {
      Int_t hit = TESTBIT(map1, iLay) + TESTBIT(map2, iLay);
      switch (hit) {
        case 1: q++; break;
        case 2: q--; break;
        default: break;
      }
    }
    SetValue(q / hits);
    return Validate();
  }

  MpdStarPairQualityCut::~MpdStarPairQualityCut() {}
}  // namespace MpdPadsFormat
