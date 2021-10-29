/*
 * MpdKinMonitor.cxx
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdKinMonitor.h"
#include "NicaTrack.h"
#include <TH1.h>
#include <TLorentzVector.h>

MpdKinMonitor::MpdKinMonitor() : NicaPropertyMonitorXY("#eta []", "p_{T} [GeV/c]", ENicaCutUpdate::kTrackUpdate) {
  SetXaxis(200, -2, 2);
  SetYaxis(200, 0, 4);
}

void MpdKinMonitor::Update(Bool_t passed, TObject* obj) {
  NicaTrack* tr = (NicaTrack*) obj;
  if (passed) {
    fHistoPassed->Fill(tr->GetMomentum().Eta(), tr->GetMomentum().Pt());
  } else {
    fHistoFailed->Fill(tr->GetMomentum().Eta(), tr->GetMomentum().Pt());
  }
}

MpdKinMonitor::~MpdKinMonitor() {
  // TODO Auto-generated destructor stub
}
