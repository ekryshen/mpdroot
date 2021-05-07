/*
 * NicaTpcPads.cxx
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaTrackTpcPads.h"

#include "NicaTpcSectorGeo.h"
#include <TLorentzVector.h>
#include <iostream>

NicaTrackTpcPads::NicaTrackTpcPads() { fPadsNo[0] = fPadsNo[1] = -2; }

NicaTrackTpcPads::~NicaTrackTpcPads() {}

NicaTrackTpcPads::NicaTrackTpcPads(const NicaTrackTpcPads& other) : TObject(other) {
  fNominalHelix = other.fNominalHelix;
  for (int i = 0; i < NicaMpdConst::TpcLayers; i++) {
    fPadID[i] = other.fPadID[i];
    fPaths[i] = other.fPaths[i];
  }
  fPadsNo[0] = other.fPadsNo[0];
  fPadsNo[1] = other.fPadsNo[1];
}

void NicaTrackTpcPads::Calculate(const NicaHelix& helix, TLorentzVector* vector) {
  if (AreCalculated()) return;  // pads are calculated
  NicaTpcSectorGeo* sec = NicaTpcSectorGeo::Instance();
  fNominalHelix         = helix;
  if (vector) {
    fNominalHelix.Shift(-vector->X(), -vector->Y(), -vector->Z());
    sec->CalculatePads(fNominalHelix, fPaths, fPadsNo);
  } else {
    sec->CalculatePads(fNominalHelix, fPaths, fPadsNo);
  }
  for (int i = fPadsNo[0]; i < fPadsNo[1]; i++) {
    TVector3 glob = fNominalHelix.EvalPos(fPaths[i]);
    TVector3 loc;
    fPadID[i] = sec->Global2Local(glob, loc, -1);
  }
  for (int i = 0; i < fPadsNo[0]; i++) {  // fill other layers
    fPadID[i] = -1;
  }
  for (int i = fPadsNo[1]; i < 53; i++) {
    fPadID[i] = -1;
  }
}

Double_t NicaTrackTpcPads::GetR(Int_t lay) const {
  NicaTpcSectorGeo* sec = NicaTpcSectorGeo::Instance();
  Double_t R            = 0;
  Double_t L            = lay;
  if (lay < sec->NofRowsReg(0)) {
    return sec->GetMinY() + sec->PadHeight(0) * (L + 0.5);
  } else {
    return sec->GetRocY(1) + sec->PadHeight(1) * (L - sec->NofRowsReg(0) + 0.5);
  }
}

Bool_t NicaTrackTpcPads::AreCalculated() const {
  if (fPadsNo[0] != -2) return kTRUE;
  return kFALSE;
}

NicaTrackTpcPads& NicaTrackTpcPads::operator=(const NicaTrackTpcPads& other) {
  if (this == &other) return *this;
  fNominalHelix = other.fNominalHelix;
  for (int i = 0; i < NicaMpdConst::TpcLayers; i++) {
    fPadID[i] = other.fPadID[i];
    fPaths[i] = other.fPaths[i];
  }
  fPadsNo[0] = other.fPadsNo[0];
  fPadsNo[1] = other.fPadsNo[1];
  return *this;
}
