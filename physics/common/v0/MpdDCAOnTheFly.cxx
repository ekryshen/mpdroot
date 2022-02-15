/*
 * MpdDCAOnTheFly.cxx
 *
 *  Created on: 14 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdDCAOnTheFly.h"

#include <FairRootManager.h>

#include "MpdEvent.h"
#include "MpdHelix.h"

MpdDCAOnTheFly::MpdDCAOnTheFly() : fEvent(nullptr), fBz(0.0) {}

InitStatus MpdDCAOnTheFly::Init()
{
   fEvent = (MpdEvent *)FairRootManager::Instance()->GetObject("MPDEvent.");
   if (fEvent) return kSUCCESS;
   return kERROR;
}

void MpdDCAOnTheFly::Exec(Option_t *option)
{
   TClonesArray *tracks = fEvent->GetGlobalTracks();
   TVector3      vertex(fEvent->GetPrimaryVerticesX(), fEvent->GetPrimaryVerticesY(), fEvent->GetPrimaryVerticesZ());
   for (int i = 0; i < tracks->GetEntriesFast(); i++) {
      MpdTrack *track = (MpdTrack *)tracks->UncheckedAt(i);
      TVector3  mom(track->GetPx(), track->GetPy(), track->GetPz());
      TVector3  pos(track->GetFirstPointX(), track->GetFirstPointY(), track->GetFirstPointZ());
      MpdHelix  h(mom, pos, track->GetCharge(), fBz);
      Double_t  s   = h.pathLength(vertex, kTRUE);
      TVector3  dca = h.at(s);
      track->SetDCAX(dca.X() - vertex.X());
      track->SetDCAX(dca.Y() - vertex.Y());
      track->SetDCAX(dca.Z() - vertex.Z());
   }
}

MpdDCAOnTheFly::~MpdDCAOnTheFly() {}
