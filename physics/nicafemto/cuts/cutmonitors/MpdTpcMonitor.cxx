/*
 * MpdTpcMonitor.cxx
 *
 *  Created on: 20 sie 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdTpcMonitor.h"
#include "NicaComplexTrack.h"
#include "NicaDataFormatManager.h"
#include "NicaMpdTrack.h"
#include <TH1.h>

MpdTpcMonitor::MpdTpcMonitor() : NicaPropertyMonitorXY("p [GeV/c]", "dEdX [AU]", ENicaCutUpdate::kTrackUpdate) {
  SetXaxis(100, 0, 4);
  SetYaxis(100, 0, 1E+5);
}

Bool_t MpdTpcMonitor::Init(Int_t task_id) { return NicaPropertyMonitorXY::Init(task_id); }

void MpdTpcMonitor::Update(Bool_t passed, TObject* obj) {
  switch (fFormatType) {
    case ENicaFormatType::kReco: {
      NicaMpdTrack* track = static_cast<NicaMpdTrack*>(obj);
      if (passed) {
        fHistoPassed->Fill(track->GetMomentum().P(), track->GetTpcTrack()->GetDeDx());
      } else {
        fHistoFailed->Fill(track->GetMomentum().P(), track->GetTpcTrack()->GetDeDx());
      }
    } break;
    case ENicaFormatType::kComplexReco: {
      NicaComplexTrack* complex_track = static_cast<NicaComplexTrack*>(obj);
      NicaMpdTrack* track             = static_cast<NicaMpdTrack*>(complex_track->GetRealTrack());
      if (passed) {
        fHistoPassed->Fill(track->GetMomentum().P(), track->GetTpcTrack()->GetDeDx());
      } else {
        fHistoFailed->Fill(track->GetMomentum().P(), track->GetTpcTrack()->GetDeDx());
      }
    } break;
    default: break;
  }
}

MpdTpcMonitor::MpdTpcMonitor(const MpdTpcMonitor& other) : NicaPropertyMonitorXY(other) {}

MpdTpcMonitor::~MpdTpcMonitor() {
  // TODO Auto-generated destructor stub
}
