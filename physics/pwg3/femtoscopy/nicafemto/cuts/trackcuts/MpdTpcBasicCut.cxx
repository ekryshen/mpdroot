/*
 * NicaTrackTpcSigmaCut.cpp
 *
 *  Created on: 28 mar 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdTpcBasicCut.h"
#include "NicaComplexTrack.h"
#include "NicaDataFormatManager.h"
#include "NicaMpdConst.h"
#include "NicaMpdTrack.h"
#include <fstream>

Double_t MpdTcpCalibratedCut::BetheBlochFunction(Double_t X, Double_t* par) {
  Float_t x     = X;
  Float_t p     = x;       // kaon momentum
  Float_t mass  = par[0];  // Mass of particle
  Float_t beta  = p / TMath::Sqrt(p * p + mass * mass);
  Float_t gamma = 1.0 / TMath::Sqrt(1 - beta * beta);
  Float_t bg    = beta * gamma;

  Float_t kp1 = par[1];
  Float_t kp2 = par[2];
  Float_t kp3 = par[3];
  Float_t kp4 = par[4];
  Float_t kp5 = par[5];

  beta       = bg / TMath::Sqrt(1. + bg * bg);
  Float_t aa = TMath::Power((double) beta, (double) kp4);
  Float_t bb = TMath::Power(1. / bg, kp5);
  // cout<<1./bg<<" "<<kp5<<endl;
  // cout<<bb<<endl;
  bb = TMath::Log(kp3 + bb);
  return (kp2 - aa - bb) * kp1 / aa / 2E+9;
}

MpdTcpCalibratedCut::MpdTcpCalibratedCut(TString calib_file) : NicaTrackTpcCut(), fUseDst(kFALSE) {
  if (calib_file.Length() == 0) {
    fUseDst = kTRUE;
    return;
  }
  std::ifstream calib;
  calib.open(calib_file.Data());
  TString trash;
  calib >> trash >> trash;
  for (int i = 0; i < 6; i++)
    calib >> fPionBB[i];
  calib >> trash >> trash;
  for (int i = 0; i < 6; i++)
    calib >> fKaonBB[i];
  calib >> trash >> trash;
  for (int i = 0; i < 6; i++)
    calib >> fProtonBB[i];
  calib >> trash >> trash;
  for (int i = 0; i < 6; i++)
    calib >> fElectronBB[i];
  for (int i = 0; i < 7; i++)
    calib >> trash;
  for (int i = 0; i < 72; i++)
    calib >> fSigmasPi[i];
  for (int i = 0; i < 7; i++)
    calib >> trash;
  for (int i = 0; i < 72; i++)
    calib >> fSigmasK[i];
  for (int i = 0; i < 7; i++)
    calib >> trash;
  for (int i = 0; i < 72; i++)
    calib >> fSigmasP[i];
  for (int i = 0; i < 7; i++)
    calib >> trash;
  for (int i = 0; i < 72; i++)
    calib >> fSigmasE[i];
  for (int i = 0; i < 7; i++)
    calib >> trash;
  for (int i = 0; i < 72; i++)
    calib >> fPInt[i];
  calib.close();
}

Bool_t MpdTcpCalibratedCut::Pass(NicaTrack* track) {
  NicaTpcTrack* tpc = (NicaTpcTrack*) ((NicaMpdTrack*) track)->GetDetTrack(MpdDetectorID::kTPC);
  SetValue(track->GetCharge(), fgChargeId);
  SetValue(tpc->GetDeDx(), fgDeDxId);
  SetValue(tpc->GetNHits(), fgTpcHitsId);
  if (fUseDst) {
    SetValue(tpc->GetSigmaElectron(), fgSigmaElectronId);
    SetValue(tpc->GetSigmaPion(), fgSigmaPionId);
    SetValue(tpc->GetSigmaKaon(), fgSigmaKaonId);
    SetValue(tpc->GetSigmaProton(), fgSigmaProtonId);
  } else {
    Double_t mom = track->GetMomentum().P();
    Int_t i_p    = 0;
    Int_t k      = 0;
    while (mom > fPInt[k]) {
      k++;
    }
    // std::cout<<">" <<k<<std::endl;
    if (k != 0 && k < 71) {
      i_p = k - 1;
    } else {
      if (k < 0)
        i_p = 0;
      else
        i_p = 71;
    }
    Double_t f_pion = BetheBlochFunction(mom, fPionBB);
    Double_t dedx   = tpc->GetDeDx();
    SetValue((dedx - BetheBlochFunction(mom, fPionBB)) / fSigmasPi[i_p], fgSigmaPionId);
    SetValue((dedx - BetheBlochFunction(mom, fKaonBB)) / fSigmasK[i_p], fgSigmaKaonId);
    SetValue((dedx - BetheBlochFunction(mom, fProtonBB)) / fSigmasP[i_p], fgSigmaProtonId);
    SetValue((dedx - BetheBlochFunction(mom, fElectronBB)) / fSigmasE[i_p], fgSigmaElectronId);
  }
  return ForcedUpdate(Verify());
}
