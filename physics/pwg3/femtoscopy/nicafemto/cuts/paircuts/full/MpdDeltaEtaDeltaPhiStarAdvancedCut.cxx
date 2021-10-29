/*
 * MpdDeltaEtaDeltaPhiStarCut.cxx
 *
 *  Created on: 4 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdDeltaEtaDeltaPhiStarAdvancedCut.h"

#include "NicaMpdTrackTpcPads.h"
#include "NicaTwoTrack.h"

namespace MpdPadsFormat {
  MpdDeltaEtaDeltaPhiStarAdvancedCut::MpdDeltaEtaDeltaPhiStarAdvancedCut() : MpdTpcPadsPairCut(3) {
    SetUnitName("#Delta#eta", DeltaEta());
    SetUnitName("#Delta#phi^{*}", DeltaPhiStar());
    SetUnitName("#Delta#phi^{*}_{min}", DeltaPhiStarMin());
  }

  Bool_t MpdDeltaEtaDeltaPhiStarAdvancedCut::Pass(NicaTwoTrack* pair) {
    NicaMpdTrackTpcPads* track1 = (NicaMpdTrackTpcPads*) pair->GetTrack1();
    NicaMpdTrackTpcPads* track2 = (NicaMpdTrackTpcPads*) pair->GetTrack2();
    const Double_t hits1        = track1->GetTpcPadsInfo()->GetPadsNo();
    const Double_t hits2        = track2->GetTpcPadsInfo()->GetPadsNo();
    Double_t minHits            = TMath::Min(hits1, hits2);
    Double_t maxHits            = TMath::Max(hits1, hits2);
    Double_t phi1, phi2;
    Double_t phi_star, phi_star_min;
    SetValue(track1->GetMomentum().Eta() - track2->GetMomentum().Eta(), DeltaEta());
    phi_star_min     = 1E+9;
    phi_star         = 0;
    Double_t Phi1    = track1->GetMomentum().Phi();
    Double_t Phi2    = track2->GetMomentum().Phi();
    Double_t ptv1    = track1->GetMomentum().Pt();
    Double_t ptv2    = track2->GetMomentum().Pt();
    Double_t charge1 = track1->GetCharge();
    Double_t charge2 = track2->GetCharge();
    Double_t magSign = 1;  // FIXME
    for (int i = 0; i < minHits; i++) {
      if (track1->GetTpcPadsInfo()->GetPadID(i) == -1) continue;
      if (track2->GetTpcPadsInfo()->GetPadID(i) == -1) continue;
      Double_t rad = track1->GetTpcPadsInfo()->GetR(i) * 0.01;

      // phi1 = track1->GetPhi(i);
      // phi2 = track2->GetPhi(i);
      Double_t dPhi = (Phi2 - Phi1 + (TMath::ASin(-0.075 * charge2 * magSign * rad / ptv2))
                       - (TMath::ASin(-0.075 * charge1 * magSign * rad / ptv1)));
      Double_t dphi = TVector2::Phi_mpi_pi(dPhi);
      if (TMath::Abs(dphi) < TMath::Abs(phi_star_min)) phi_star_min = dphi;
      phi_star += dphi;
    }
    SetValue(phi_star / minHits, DeltaPhiStar());
    SetValue(phi_star_min, DeltaPhiStarMin());
    if (InLimits(DeltaEta()) && InLimits(DeltaPhiStar()) && InLimits(DeltaPhiStarMin())) { return ForcedUpdate(kFALSE); }
    return ForcedUpdate(kTRUE);
  }

  MpdDeltaEtaDeltaPhiStarAdvancedCut::~MpdDeltaEtaDeltaPhiStarAdvancedCut() {
    // TODO Auto-generated destructor stub
  }
}  // namespace MpdPadsFormat
