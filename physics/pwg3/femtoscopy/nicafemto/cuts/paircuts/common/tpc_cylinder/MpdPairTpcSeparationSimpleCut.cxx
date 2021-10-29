/*
 * MpdPairTpcSeparatioSimpleCut.cxx
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairTpcSeparationSimpleCut.h"

#include "NicaTpcSectorGeo.h"
#include "NicaTwoTrack.h"

MpdPairTpcSeparationSimpleCut::MpdPairTpcSeparationSimpleCut(
    Int_t averaging_steps)
    : MpdPairTpcSimpleCut(5),
      fNSteps(averaging_steps),
      fSteps(averaging_steps),
      fRMin(0),
      fRMax(0),
      fZMax(0),
      fRadii(NULL) {
  SetUnitName("TPC Entry Sep. [cm]", Entry());
  SetUnitName("TPC Exit Sep [cm]", Exit());
  SetUnitName("TPC Average Sep. [cm]", Average());
  SetUnitName("TPC Min Sep [cm]", Min());
  SetUnitName("TPC Max Sep [cm]", Max());
  for (int i = 0; i < 5; i++) SetMinMax(0, 1E+5, i);
}

MpdPairTpcSeparationSimpleCut::~MpdPairTpcSeparationSimpleCut() {
  if (fRadii) delete[] fRadii;
}

Bool_t MpdPairTpcSeparationSimpleCut::Pass(NicaTwoTrack *pair) {
  NicaExpTrackHelix *track1 = (NicaExpTrackHelix *)pair->GetTrack1();
  NicaExpTrackHelix *track2 = (NicaExpTrackHelix *)pair->GetTrack2();
  TVector3 pos1 = PosNominal(track1, fRMin);
  TVector3 pos2 = PosNominal(track2, fRMin);
  pos1 -= pos2;
  SetValue(pos1.Mag(), Entry());
  pos1 = PosNominal(track1, fRMax);
  pos2 = PosNominal(track2, fRMax);
  pos1 -= pos2;
  SetValue(pos1.Mag(), Exit());
  Double_t min = 1E+5;
  Double_t max = 0;
  Double_t steps = 0;
  Double_t sum = 0;
  for (int i = 0; i < fNSteps; i++) {
    pos1 = PosNominal(track1, fRadii[i]);
    if (pos1.Mag2() == 0) break;
    if (TMath::Abs(pos1.Z()) > fZMax) break;
    pos2 = PosNominal(track2, fRadii[i]);
    if (TMath::Abs(pos2.Z()) > fZMax) break;
    if (pos2.Mag2() == 0) break;

    pos1 -= pos2;
    Double_t dt = pos1.Mag();
    min = TMath::Min(dt, min);
    max = TMath::Max(dt, max);
    sum += dt;
    steps++;
  }
  if (steps > 0)
    SetValue(sum / steps, Average());
  else
    SetValue(0, Average());
  SetValue(max, Max());
  SetValue(min, Min());
  return Validate();
}

Bool_t MpdPairTpcSeparationSimpleCut::Init(Int_t format_id) {
  if (fSteps < 0) fSteps = 10;
  if (fSteps > 60) fSteps = 50;
  fSteps = (Int_t)fSteps;
  fNSteps = fSteps;
  NicaTpcSectorGeo *geo = NicaTpcSectorGeo::Instance();
  fRMin = geo->GetMinY();
  fRMax = geo->GetMaxY();
  Double_t dR = (fRMax - fRMin) / fSteps;
  fRadii = new Float_t[fNSteps];
  fZMax = geo->GetZmax();
  for (int i = 0; i < fNSteps; i++) {
    fRadii[i] = fRMin + dR * i;
  }
  return MpdPairTpcSimpleCut::Init(format_id);
}

MpdPairTcpSeparationSimpleNegCut::MpdPairTcpSeparationSimpleNegCut() {
  for (int i = 0; i < 5; i++) SetMinMax(-1E+9, 1E+9, i);
}

Bool_t MpdPairTcpSeparationSimpleNegCut::Pass(NicaTwoTrack *pair) {
  NicaExpTrackHelix *track1 = (NicaExpTrackHelix *)pair->GetTrack1();
  NicaExpTrackHelix *track2 = (NicaExpTrackHelix *)pair->GetTrack2();
  TVector3 pos1 = PosNominal(track1, fRMin);
  TVector3 pos2 = PosNominal(track2, fRMin);
  pos1 -= pos2;
  SetValue(pos1.Mag(), Entry());
  pos1 = PosNominal(track1, fRMax);
  pos2 = PosNominal(track2, fRMax);
  pos1 -= pos2;
  SetValue(pos1.Mag(), Exit());
  Double_t min = 1E+5;
  Double_t max = 0;
  Double_t steps = 0;
  Double_t sum = 0;
  for (int i = 0; i < fSteps; i++) {
    pos1 = PosNominal(track1, fRadii[i]);
    if (pos1.Mag2() == 0) break;
    if (TMath::Abs(pos1.Z()) > fZMax) break;
    pos2 = PosNominal(track2, fRadii[i]);
    if (TMath::Abs(pos2.Z()) > fZMax) break;
    if (pos2.Mag2() == 0) break;
    pos1 -= pos2;
    Double_t dt = pos1.Mag();
    min = TMath::Min(dt, min);
    max = TMath::Max(dt, max);
    sum += dt;
    steps++;
  }
  if (steps > 0)
    SetValue(sum / steps, Average());
  else
    SetValue(0, Average());
  Double_t maxu = TMath::Max(min, TMath::Max(sum / steps, GetValue(Entry())));
  SetValue(max, Max());
  SetValue(min, Min());
  Double_t gMin = TMath::Min(GetValue(Entry()),
                             TMath::Min(GetValue(Max()), GetValue(Average())));
  Double_t gMax = TMath::Max(GetValue(Entry()),
                             TMath::Max(GetValue(Max()), GetValue(Average())));
  // SetValue(gMin,5);
  // SetValue(gMax,6);
  if (InLimits(Average())) return ForcedUpdate(kTRUE);
  if (InLimits(Max())) return ForcedUpdate(kTRUE);
  if (InLimits(Min())) return ForcedUpdate(kTRUE);
  if (InLimits(Exit())) return ForcedUpdate(kTRUE);
  if (InLimits(Entry())) return ForcedUpdate(kTRUE);
  return Validate();
}

MpdPairTcpSeparationSimpleNegCut::~MpdPairTcpSeparationSimpleNegCut() {}
