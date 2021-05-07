/*
 * NicaMpdMiniDstTrack.cxx
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "NicaMpdMiniDstTrack.h"


#include <TLorentzVector.h>
#include <TVector3.h>

#include "MpdMiniTrack.h"
#include "NicaEvent.h"
#include "NicaExpTrack.h"
#include "NicaHelix.h"
#include "NicaToFTrack.h"
#include "NicaTpcTrack.h"
#include "NicaTrack.h"

NicaMpdMiniDstTrack::NicaMpdMiniDstTrack() {}

void NicaMpdMiniDstTrack::Update(MpdMiniTrack* track, NicaMpdMiniDstEvent::eMode mode) {
  TVector3 gMom, mom;
  TVector3 dca = track->origin();
  gMom         = track->gMom();
  fType        = 0;
  switch (mode) {
    case NicaMpdMiniDstEvent::eMode::kGlobalTrack: {
      mom = track->gMom();
      SetGlobal(kTRUE);
    } break;
    case NicaMpdMiniDstEvent::eMode::kPrimaryTrack: {
      SetGlobal(kFALSE);
      if (track->isPrimary()) {
        mom = track->pMom();
      } else {
        mom = track->gMom();
      }

    } break;
  }
  SetMomentum(mom.X(), mom.Y(), mom.Z(), 0);
  SetDCA(dca.X(), dca.Y(), dca.Z());

  SetID(track->id());
  SetCharge(track->charge());
  SetNHits(track->nHits());
  SetHitMap(track->hitMap());

  if (!track->isPrimary()) { SetGlobal(kTRUE); }
  SetChi2(track->chi2());
  fTpcTrack->SetNHits(track->nHits());  //! FIXME
  fTpcTrack->SetDeDx(track->dEdx());
  fTpcTrack->SetSigma(track->nSigmaPion(), track->nSigmaKaon(), track->nSigmaProton(), track->nSigmaProton());

  fToFTrack->SetBeta(NicaToFTrack::DummyVal());
  fToFTrack->SetMass2(NicaToFTrack::DummyVal());
  fToFTrack->SetFlag(0);
  fHitsMap = track->hitMap();

  TVector3 origin = track->origin();
  fFirstPoint->SetXYZ(origin.X(), origin.Y(), origin.Z());
  Double_t vx = GetEvent()->GetVertex()->X();
  Double_t vy = GetEvent()->GetVertex()->Y();
  Double_t vz = GetEvent()->GetVertex()->Z();
  SetDCA(origin.X() - vx, origin.Y() - vy, origin.Z());
  switch (mode) {
    case NicaMpdMiniDstEvent::eMode::kGlobalTrack: {
      SetDCA(origin.X() - vx, origin.Y() - vy, origin.Z() - vz);
    } break;
    case NicaMpdMiniDstEvent::eMode::kPrimaryTrack: {
      SetDCA(0, 0, 0);
      // GetDCA()->SetXYZ(origin.X()-vx, origin.Y()-vy,origin.Z());
    } break;
  }
  // fSharedHitsMap = track->topologyMap(0);
  // fHitsMap = track->topologyMap(1);
  UpdateHelix();
}

NicaMpdMiniDstTrack::~NicaMpdMiniDstTrack() {
  // TODO Auto-generated destructor stub
}

//==================================================
NicaMpdMiniDstTrackInterface::NicaMpdMiniDstTrackInterface() {
  // TODO Auto-generated constructor stub
}

NicaMpdMiniDstTrackInterface::~NicaMpdMiniDstTrackInterface() {
  // TODO Auto-generated destructor stub
}
