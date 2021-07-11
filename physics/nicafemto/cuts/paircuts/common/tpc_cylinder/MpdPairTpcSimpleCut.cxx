/*
 * MpdPairTpcSimpleCut..cxx
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairTpcSimpleCut.h"
#include "NicaEvent.h"
#include "NicaHelix.h"

TVector3 MpdPairTpcSimpleCut::PosNominal(NicaExpTrackHelix* track, Double_t R) const {
  const NicaHelix& helix = track->GetHelix();
  Double_t s1, s2;
  NicaEvent* ev = track->GetEvent();
  helix.PathLength(R, ev->GetVertex()->X(), ev->GetVertex()->Y(), s1, s2);
  if (s1 == s2) { return TVector3(0, 0, 0); }
  while (s1 < 0) {
    s1 += helix.GetPeriod();
  }
  while (s2 < 0) {
    s2 += helix.GetPeriod();
  }
  TVector3 p = helix.EvalPos(TMath::Min(s1, s2));
  // shift by vertex pos
  p.SetXYZ(p.X() - ev->GetVertex()->X(), p.Y() - ev->GetVertex()->Y(), p.Z() - ev->GetVertex()->Z());
  return p;
}

TVector3 MpdPairTpcSimpleCut::PosReal(NicaExpTrackHelix* track, Double_t R) const {
  const NicaHelix& helix = track->GetHelix();
  Double_t s1, s2;
  NicaEvent* ev = track->GetEvent();
  helix.PathLength(R, s1, s2);
  if (s1 == s2) return TVector3(0, 0, 0);
  while (s1 < 0) {
    s1 += helix.GetPeriod();
  }
  while (s2 < 0) {
    s2 += helix.GetPeriod();
  }
  return helix.EvalPos(TMath::Min(s1, s2));
}

MpdPairTpcSimpleCut::MpdPairTpcSimpleCut(Int_t size) : MpdCylinderTpcPairCut(size) {}

MpdPairTpcSimpleCut::~MpdPairTpcSimpleCut() {}
