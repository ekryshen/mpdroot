/*
 * NicaMpdTrack.cxx
 *
 *  Created on: 28 mar 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdTrack.h"
#include "NicaDataFormat.h"
#include <Rtypes.h>
#include <iostream>
NicaMpdTrack::NicaMpdTrack() : NicaExpTrackHelix(), fHitsMap(0), fSharedHitsMap(0) {
  fFirstPoint = new TVector3();
  fLastPoint  = new TVector3();
  fTpcTrack   = new NicaTpcTrack();
  fToFTrack   = new NicaToFTrack();
}

NicaMpdTrack::~NicaMpdTrack() {
  delete fFirstPoint;
  delete fLastPoint;
  delete fTpcTrack;
  delete fToFTrack;
}

void NicaMpdTrack::Update(MpdTrack* track) {
  SetMomentum(track->GetPx(), track->GetPy(), track->GetPz(), 0);
  SetID(track->GetID());
  if (track->GetPt() > 0)
    SetCharge(-1);
  else
    SetCharge(1);
  SetPrimary();
  SetNHits(track->GetNofHits());
  fTpcTrack->SetNHits(track->GetNofHits());  //! FIXME
  fTpcTrack->SetPidProb(
    track->GetTPCPidProbPion(), track->GetTPCPidProbKaon(), track->GetTPCPidProbProton(), track->GetTPCPidProbElectron());
  fToFTrack->SetBeta(track->GetTofBeta());
  fToFTrack->SetMass2(track->GetTofMass2());
  fTpcTrack->SetDeDx(track->GetdEdXTPC());
  fToFTrack->SetFlag(track->GetTofFlag());
  fTpcTrack->SetSigma(track->GetNSigmaPion(), track->GetNSigmaKaon(), track->GetNSigmaProton(), track->GetNSigmaProton());
  SetChi2(track->GetChi2());
  SetDCA(track->GetDCAX(), track->GetDCAY(), track->GetDCAZ());
  fFirstPoint->SetXYZ(track->GetFirstPointX(), track->GetFirstPointY(), track->GetFirstPointZ());
  fLastPoint->SetXYZ(track->GetLastPointX(), track->GetLastPointY(), track->GetLastPointZ());
  fSharedHitsMap = track->GetSharedHitMap();
  fHitsMap       = track->GetLayerHitMap();
  UpdateHelix();
}

void NicaMpdTrack::CopyData(NicaTrack* other) {
  NicaMpdTrack* track = (NicaMpdTrack*) other;
  *fFirstPoint        = *track->fFirstPoint;
  *fLastPoint         = *track->fLastPoint;
  *fTpcTrack          = *track->fTpcTrack;
  *fToFTrack          = *track->fToFTrack;
  fHitsMap            = track->fHitsMap;
  fSharedHitsMap      = track->fSharedHitsMap;
  NicaExpTrackHelix::CopyData(track);
}

NicaMpdTrack::NicaMpdTrack(const NicaMpdTrack& other) : NicaExpTrackHelix(other) {
  fFirstPoint    = new TVector3(*other.fFirstPoint);
  fLastPoint     = new TVector3(*other.fLastPoint);
  fTpcTrack      = new NicaTpcTrack(*other.fTpcTrack);
  fToFTrack      = new NicaToFTrack(*other.fToFTrack);
  fHitsMap       = other.fHitsMap;
  fSharedHitsMap = other.fSharedHitsMap;
}

NicaMpdTrack& NicaMpdTrack::operator=(const NicaMpdTrack& other) {
  if (this == &other) return *this;
  NicaExpTrackHelix::operator=(other);
  *fFirstPoint               = *other.fFirstPoint;
  *fLastPoint                = *other.fLastPoint;
  *fTpcTrack                 = *other.fTpcTrack;
  *fToFTrack                 = *other.fToFTrack;
  fHitsMap                   = other.fHitsMap;
  fSharedHitsMap             = other.fSharedHitsMap;
  return *this;
}

void NicaMpdTrack::PrintHitMaps() const {
  std::cout << "HIT MAP:\t";
  for (int i = 0; i < 53; i++) {
    Bool_t bit = TESTBIT(fHitsMap, i);
    if (bit) {
      std::cout << '1';
    } else {
      std::cout << '0';
    }
  }
  std::cout << std::endl;
  std::cout << "SHARED HIT MAP\t";
  for (int i = 0; i < 53; i++) {
    Bool_t bit = TESTBIT(fSharedHitsMap, i);
    if (bit) {
      std::cout << '1';
    } else {
      std::cout << '0';
    }
  }
  std::cout << std::endl;
}

NicaDetectorTrack* NicaMpdTrack::GetDetTrack(const UInt_t detID) const {
  switch (detID) {
    case NicaDetectorID::kTPC: {
      return fTpcTrack;
    } break;
    case NicaDetectorID::kTOF: {
      return fToFTrack;
    } break;
    default: return nullptr; break;
  }
}
