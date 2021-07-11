/*
 * MpdPairDeltaDCA.cxx
 *
 *  Created on: 9 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairDeltaDCA.h"

#include "NicaMpdTrack.h"
#include "NicaTwoTrack.h"

MpdPairDeltaDCA::MpdPairDeltaDCA() : NicaTwoTrackCut(3) {
  SetUnitName("#DeltaDCA_{xy} [cm]", DeltaDCAxy());
  SetUnitName("#DeltaDCA_{z} [cm]", DeltaDCAz());
  SetUnitName("#DeltaDCA [cm]", DeltaDCA());
  SetMinMax(0, 1E+6, DeltaDCAxy());
  SetMinMax(0, 1E+6, DeltaDCAz());
  SetMinMax(0, 1E+6, DeltaDCA());
}

Bool_t MpdPairDeltaDCA::Pass(NicaTwoTrack* pair) {
  NicaExpTrack* tr1    = (NicaExpTrack*) pair->GetTrack1();
  NicaExpTrack* tr2    = (NicaExpTrack*) pair->GetTrack2();
  const TVector3& dca1 = tr1->GetDCA();
  const TVector3& dca2 = tr2->GetDCA();
  TVector3 deltaDCA    = dca1 - dca2;
  SetValue(deltaDCA.Pt(), DeltaDCAxy());
  SetValue(TMath::Abs(deltaDCA.Z()), DeltaDCAz());
  SetValue(deltaDCA.Mag(), DeltaDCA());
  return Validate();
}

MpdPairDeltaDCA::~MpdPairDeltaDCA() {}
