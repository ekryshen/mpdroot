/*
 * MpdTofMonitor.cxx
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdTofMonitor.h"
#include "NicaComplexTrack.h"
#include "NicaMpdTrack.h"
#include <TH1.h>

MpdTofMonitor::MpdTofMonitor() : NicaPropertyMonitorXY("p [GeV/c]", "m^{2}_{TOF}", ENicaCutUpdate::kTrackUpdate) {
  SetXaxis(100, 0, 4);
  SetYaxis(150, -0.25, 1.25);
}

void MpdTofMonitor::Update(Bool_t passed, TObject* obj) {
  switch (fFormatType) {
    case ENicaFormatType::kReco: {
      NicaMpdTrack* track = static_cast<NicaMpdTrack*>(obj);
      if (passed) {
        fHistoPassed->Fill(track->GetMomentum().P(), track->GetToFTrack()->GetMass2());
      } else {
        fHistoFailed->Fill(track->GetMomentum().P(), track->GetToFTrack()->GetMass2());
      }
    } break;
    case ENicaFormatType::kComplexReco: {
      NicaComplexTrack* complex_track = static_cast<NicaComplexTrack*>(obj);
      NicaMpdTrack* track             = static_cast<NicaMpdTrack*>(complex_track->GetRealTrack());
      if (passed) {
        fHistoPassed->Fill(track->GetMomentum().P(), track->GetToFTrack()->GetMass2());
      } else {
        fHistoFailed->Fill(track->GetMomentum().P(), track->GetToFTrack()->GetMass2());
      }
    } break;
    default: break;
  }
}

MpdTofMonitor::~MpdTofMonitor() {}
