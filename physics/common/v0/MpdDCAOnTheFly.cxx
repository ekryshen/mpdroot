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

#include <iostream>
#include "MpdEvent.h"
#include "MpdMiniPhysicalHelix.h"
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
      MpdTrack *           track = (MpdTrack *)tracks->UncheckedAt(i);
      TVector3             mom(track->GetPx(), track->GetPy(), track->GetPz());
      TVector3             pos(track->GetFirstPointX(), track->GetFirstPointY(), track->GetFirstPointZ());
      MpdMiniPhysicalHelix helix2(mom, pos, 5 * 1e-14, track->GetCharge());
      MpdHelix             h(mom, pos, track->GetCharge(), 0.5);
      Double_t             s   = h.pathLength(vertex, kTRUE);
      TVector3             dca = h.cat(s);
      // std::cout << "DCA before" << std::endl;
      // std::cout << pos.X() << " " << pos.Y() << " " << pos.Z() << std::endl;
      // std::cout << track->GetDCAX() << " " << track->GetDCAY() << " " << track->GetDCAZ() << " " << track <<
      // std::endl;
      /* track->SetDCAX(dca.X() - vertex.X());
       track->SetDCAY(dca.Y() - vertex.Y());
       track->SetDCAZ(dca.Z() - vertex.Z());
       // std::cout << "CurvHelix = " << h.curvature() << std::endl;
       std::cout << track->GetDCAX() << " " << track->GetDCAY() << " " << track->GetDCAZ() << std::endl;

       s   = helix2.pathLength(vertex, kTRUE);
       dca = helix2.at(s);
       track->SetDCAX(dca.X() - vertex.X());
       track->SetDCAY(dca.Y() - vertex.Y());
       track->SetDCAZ(dca.Z() - vertex.Z());
       // std::cout << "CurvmINIHelix = " << helix2.curvature() << std::endl;
       std::cout << track->GetDCAX() << " " << track->GetDCAY() << " " << track->GetDCAZ() << std::endl;
 */
      /*  std::cout << "Vertex" << s << std::endl;
        std::cout << fEvent->GetPrimaryVerticesX() << " " << fEvent->GetPrimaryVerticesY() << " "
                  << fEvent->GetPrimaryVerticesZ() << std::endl;*/

      track->SetPhi(dca.Phi());
      track->SetTheta(dca.Theta());
   }
}

MpdDCAOnTheFly::~MpdDCAOnTheFly() {}
