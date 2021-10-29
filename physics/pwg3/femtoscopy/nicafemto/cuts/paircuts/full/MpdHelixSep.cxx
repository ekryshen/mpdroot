/*
 * MpdHelixSep.cxx
 *
 *  Created on: 18 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdHelixSep.h"

#include "NicaMpdTrackTpcPads.h"
#include "NicaTpcTrack.h"
#include "NicaTwoTrack.h"
namespace MpdPadsFormat {
  MpdHelixSep::MpdHelixSep() : MpdTpcPadsPairCut(1) { SetUnitName("DistQuality[%"); }

  Bool_t MpdHelixSep::Pass(NicaTwoTrack* pair) {
    NicaMpdTrackTpcPads* track1 = (NicaMpdTrackTpcPads*) pair->GetTrack1();
    NicaMpdTrackTpcPads* track2 = (NicaMpdTrackTpcPads*) pair->GetTrack2();
    track1->CalculatePads();
    track2->CalculatePads();
    const Double_t hits1 = ((NicaTpcTrack*) track1->GetDetTrack(NicaDetectorID::kTPC))->GetNHits();
    const Double_t hits2 = ((NicaTpcTrack*) track2->GetDetTrack(NicaDetectorID::kTPC))->GetNHits();
    Int_t padmin         = TMath::Min(track1->GetTpcPadsInfo()->GetFirstGoodPad(), track2->GetTpcPadsInfo()->GetFirstGoodPad());
    if (padmin >= 0) {
      TVector3 ent1        = track1->GetTpcPadsInfo()->GetMomentumAtLayer(padmin);
      TVector3 ent2        = track1->GetTpcPadsInfo()->GetMomentumAtLayer(padmin);
      Double_t entry_dist  = (ent1 - ent2).Mag();
      const TVector3& dca1 = track1->GetDCA();
      const TVector3& dca2 = track2->GetDCA();
      Double_t dca_dist    = (dca1 - dca2).Mag();
      if (dca_dist == 0) dca_dist = 1E-6;
      SetValue(entry_dist / dca_dist);
    } else {
      SetValue(100);
    }
    return Validate();
  }

  MpdHelixSep::~MpdHelixSep() {
    // TODO Auto-generated destructor stub
  }
}  // namespace MpdPadsFormat
