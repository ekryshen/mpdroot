/*
 * MpdPairDeltaPhiStarDeltaEtaCut.cxx
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "MpdPairDeltaPhiStarDeltaEtaCut.h"

#include "NicaDataFormatManager.h"
#include "NicaExpEvent.h"
#include "NicaHelix.h"
#include "NicaMpdConst.h"
#include "NicaPackage.h"
#include "NicaParameter.h"
#include "NicaTwoTrack.h"
#include <RtypesCore.h>
#include <iostream>

MpdPairDeltaPhiStarDeltaEtaCut::MpdPairDeltaPhiStarDeltaEtaCut() : NicaTwoTrackCut(2), fR(0.27) {
  SetUnitName("#Delta#phi^{*}[rad]");
  SetUnitName("#Delta#eta []", 1);
}

MpdPairDeltaPhiStarDeltaEtaCut::MpdPairDeltaPhiStarDeltaEtaCut(const MpdPairDeltaPhiStarDeltaEtaCut& other) :
  NicaTwoTrackCut(other), fR(other.fR) {}

Bool_t MpdPairDeltaPhiStarDeltaEtaCut::Pass(NicaTwoTrack* pair) {
  Double_t MagField = static_cast<NicaExpEvent*>(pair->GetTrack1()->GetEvent())->GetMagField()->Z();
  Double_t phi1     = pair->GetTrack1()->GetMomentum().Phi();
  Double_t phi2     = pair->GetTrack2()->GetMomentum().Phi();
  Double_t chg1     = pair->GetTrack1()->GetCharge();
  Double_t chg2     = pair->GetTrack2()->GetCharge();
  Double_t ptv1     = pair->GetTrack1()->GetMomentum().Pt();
  Double_t ptv2     = pair->GetTrack2()->GetMomentum().Pt();
  Double_t eta1     = pair->GetTrack1()->GetMomentum().Eta();
  Double_t eta2     = pair->GetTrack2()->GetMomentum().Eta();
  SetValue(eta1 - eta2, 1);
  // B*R
  const Double_t scale = -0.149896229 * MagField * fR;
  Double_t afsi0b      = scale * chg1 / ptv1;
  Double_t afsi1b      = scale * chg2 / ptv2;
  Double_t dps         = phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
  dps                  = TVector2::Phi_mpi_pi(dps);
  SetValue(dps, 0);
  if (TMath::Abs(afsi0b) >= 0 || TMath::Abs(afsi1b) >= 0) { return ForcedUpdate(kTRUE); }
  if (InLimits(0) && InLimits(1)) return ForcedUpdate(kFALSE);
  return ForcedUpdate(kTRUE);
}

NicaPackage* MpdPairDeltaPhiStarDeltaEtaCut::Report() const {
  NicaPackage* report = NicaTwoTrackCut::Report();
  report->AddObject(new NicaParameterDouble("R [cm]", fR * 100.0));
  return report;
}

void MpdPairDeltaPhiStarDeltaEtaCut::SetDeltaEtaCut(Double_t min, Double_t max) { SetMinMax(min, max, 0); }

void MpdPairDeltaPhiStarDeltaEtaCut::SetDeltaPhiStarCut(Double_t min, Double_t max) { SetMinMax(min, max, 1); }

Bool_t MpdPairDeltaPhiStarDeltaEtaCut::Init(Int_t task_id) {
  Bool_t fine = FormatInhertis("NicaExpEvent", task_id, ENicaFormatDepth::kBuffered);
  const NicaExpEvent* ev =
    (const NicaExpEvent*) NicaDataFormatManager::Instance()->GetFormat(task_id, ENicaFormatDepth::kBuffered);
  return fine;
}

MpdPairDeltaPhiStarDeltaEtaCut::~MpdPairDeltaPhiStarDeltaEtaCut() {}

Bool_t MpdPairDeltaPhiStarDeltaEtaMinCut::Pass(NicaTwoTrack* pair) {
  Double_t MagField = static_cast<NicaExpEvent*>(pair->GetTrack1()->GetEvent())->GetMagField()->Z();
  Double_t phi1     = pair->GetTrack1()->GetMomentum().Phi();
  Double_t phi2     = pair->GetTrack2()->GetMomentum().Phi();
  Double_t chg1     = pair->GetTrack1()->GetCharge();
  Double_t chg2     = pair->GetTrack2()->GetCharge();
  Double_t ptv1     = pair->GetTrack1()->GetMomentum().Pt();
  Double_t ptv2     = pair->GetTrack2()->GetMomentum().Pt();
  Double_t eta1     = pair->GetTrack1()->GetMomentum().Eta();
  Double_t eta2     = pair->GetTrack2()->GetMomentum().Eta();
  SetValue(eta1 - eta2, 1);
  Double_t RphiMin = 1;
  for (Double_t R = NicaMpdConst::TpcInnerDriftRadius * 0.01; R < NicaMpdConst::TPcOuterDriftRadius * 0.01; R += 0.01) {
    Double_t deps =
      (phi2 - phi1 + (TMath::ASin(-0.15 * MagField * chg2 * R / ptv2)) - (TMath::ASin(-0.15 * MagField * chg1 * R / ptv1)));
    deps = TVector2::Phi_mpi_pi(deps);
    if (TMath::Abs(deps) < RphiMin) { RphiMin = deps; }
  }
  SetValue(RphiMin, 0);

  if (InLimits(0) && InLimits(1)) return ForcedUpdate(kFALSE);
  return ForcedUpdate(kTRUE);
}

NicaPackage* MpdPairDeltaPhiStarDeltaEtaMinCut::Report() const { return NicaTwoTrackCut::Report(); }
