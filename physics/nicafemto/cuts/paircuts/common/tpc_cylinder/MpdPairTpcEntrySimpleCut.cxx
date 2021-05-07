/*
 * MpdPairTpcEntrySimpleCut.cxx
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairTpcEntrySimpleCut.h"

#include "NicaTpcSectorGeo.h"
#include "NicaTwoTrack.h"

MpdPairTpcEntrySimpleCut::MpdPairTpcEntrySimpleCut()
    : MpdPairTpcSimpleCut(3), fRmin(40.3) {
  SetUnitName("Tpc Entry Sep R [cm]", XYZ());
  SetUnitName("Tpc Entry Sep XY [cm]", XY());
  SetUnitName("Tpc Entry Sep Z [cm]", Z());
}

MpdPairTpcEntrySimpleCut::~MpdPairTpcEntrySimpleCut() {}

Bool_t MpdPairTpcEntrySimpleCut::Pass(NicaTwoTrack *pair) {
  NicaExpTrackHelix *track1 = (NicaExpTrackHelix *)pair->GetTrack1();
  NicaExpTrackHelix *track2 = (NicaExpTrackHelix *)pair->GetTrack2();
  TVector3 pos1 = PosNominal(track1, fRmin);
  TVector3 pos2 = PosNominal(track2, fRmin);
  pos1 -= pos2;
  SetValue(pos1.Pt(), XY());
  SetValue(TMath::Abs(pos1.Z()), Z());
  SetValue(pos1.Mag(), XYZ());
  return Validate();
}

Bool_t MpdPairTpcEntrySimpleCut::Init(Int_t format_id) {
  NicaTpcSectorGeo *geo = NicaTpcSectorGeo::Instance();
  fRmin = geo->GetMinY();
  return NicaTwoTrackCut::Init(format_id);
}
