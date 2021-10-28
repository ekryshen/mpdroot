/*
 * MpdPairDeltaPhistarMinCut.cxx
 *
 *  Created on: 22 lip 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairDeltaPhistarMinCut.h"
#include "NicaExpEvent.h"
#include "NicaMpdConst.h"
#include "NicaPackage.h"
#include "NicaParameter.h"
#include "NicaTwoTrack.h"
#include "NicaTwoTrackCut.h"


MpdPairDeltaPhistarMinCut::MpdPairDeltaPhistarMinCut() : NicaTwoTrackCut(1) {
  SetUnitName("#Delta phi^{*}_{min} [rad]", 0);
  SetNSteps(10);
}

Bool_t MpdPairDeltaPhistarMinCut::Pass(NicaTwoTrack* pair) {
  Double_t MagField = static_cast<NicaExpEvent*>(pair->GetTrack1()->GetEvent())->GetMagField()->Z();
  Double_t phi1     = pair->GetTrack1()->GetMomentum().Phi();
  Double_t phi2     = pair->GetTrack2()->GetMomentum().Phi();
  Double_t chg1     = pair->GetTrack1()->GetCharge();
  Double_t chg2     = pair->GetTrack2()->GetCharge();
  Double_t ptv1     = pair->GetTrack1()->GetMomentum().Pt();
  Double_t ptv2     = pair->GetTrack2()->GetMomentum().Pt();
  Double_t eta1     = pair->GetTrack1()->GetMomentum().Eta();
  Double_t eta2     = pair->GetTrack2()->GetMomentum().Eta();
  // B*R
  const Double_t scale1 = -0.00299792458 * MagField * chg1 / ptv1;  // 100 to convert units
  const Double_t scale2 = -0.00299792458 * MagField * chg2 / ptv2;

  Double_t min = 1E+6;
  for (double R = NicaMpdConst::TpcInnerDriftRadius; R < NicaMpdConst::TPcOuterDriftRadius; R += fStep) {
    //   const Double_t scale = -0.149896229 * MagField * R * 2.0 * 0.01;
    Double_t afsi0b = scale1 * R;
    Double_t afsi1b = scale2 * R;
    Double_t tphi1  = phi1 + TMath::ASin(afsi0b);
    Double_t tphi2  = phi2 + TMath::ASin(afsi1b);
    Double_t dps    = tphi1 - tphi2;
    dps             = TVector2::Phi_mpi_pi(dps);
    if (TMath::Abs(min) > TMath::Abs(dps)) min = dps;
  }
  SetValue(min);
  if (InLimits(0)) return ForcedUpdate(kFALSE);
  return ForcedUpdate(kTRUE);
}

void MpdPairDeltaPhistarMinCut::SetNSteps(Int_t steps) {
  Double_t rMin = NicaMpdConst::TpcInnerDriftRadius;
  Double_t rMax = NicaMpdConst::TPcOuterDriftRadius;
  Double_t dx   = steps;
  fStep         = (rMax - rMin) / dx;
}

NicaPackage* MpdPairDeltaPhistarMinCut::Report() const {
  NicaPackage* pack = NicaTwoTrackCut::Report();
  Double_t rMin     = NicaMpdConst::TpcInnerDriftRadius;
  Double_t rMax     = NicaMpdConst::TPcOuterDriftRadius;
  Int_t steps       = (rMax - rMin) / fStep;
  pack->AddObject(new NicaParameterInt("NSteps", steps));
  return pack;
}

MpdPairDeltaPhistarMinCut::~MpdPairDeltaPhistarMinCut() {
  // TODO Auto-generated destructor stub
}
