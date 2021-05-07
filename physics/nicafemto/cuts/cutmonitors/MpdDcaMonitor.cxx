/*
 * MpdDCAMonitor.cxx
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdDcaMonitor.h"

#include "NicaComplexTrack.h"
#include "NicaExpTrack.h"

#include <TH1.h>

MpdDcaMonitor::MpdDcaMonitor() : NicaPropertyMonitorXY("DCA_{xy} [cm]", "DCA_{z} [cm]", ENicaCutUpdate::kTrackUpdate) {
  SetXaxis(100, 0, 10);
  SetYaxis(200, -10, 10);
}

void MpdDcaMonitor::Update(Bool_t passed, TObject* obj) {
  switch (fFormatType) {
    case ENicaFormatType::kReco: {
      NicaExpTrack* track = (NicaExpTrack*) obj;
      Double_t dcaxy      = track->GetDCA().Pt();
      if (passed) {
        fHistoPassed->Fill(track->GetDCA().Pt(), track->GetDCA().Z());
      } else {
        fHistoFailed->Fill(track->GetDCA().Pt(), track->GetDCA().Z());
      }
    } break;
    case ENicaFormatType::kComplexReco: {
      NicaComplexTrack* z_track = (NicaComplexTrack*) obj;
      NicaExpTrack* track       = (NicaExpTrack*) z_track->GetRealTrack();
      Double_t dcaxy            = track->GetDCA().Pt();
      if (passed) {
        fHistoPassed->Fill(track->GetDCA().Pt(), track->GetDCA().Z());
      } else {
        fHistoFailed->Fill(track->GetDCA().Pt(), track->GetDCA().Z());
      }
    } break;
    default: break;
  }
}

MpdDcaMonitor::~MpdDcaMonitor() {}
