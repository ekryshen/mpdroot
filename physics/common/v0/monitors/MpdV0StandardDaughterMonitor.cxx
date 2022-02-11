/*
 * MpdV0StandardDaughetrMonitor.cxx
 *
 *  Created on: 11 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0StandardDaughterMonitor.h"

#include <TH2.h>

#include "MpdTrack.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"

MpdV0StandardDaughterMonitor::MpdV0StandardDaughterMonitor() : MpdV0DaughterMonitor(1, 1)
{
   SetDCAAxis(100, 0, 10);
   SetTpcDeDxXaxis(100, 0, 5);
   SetTpcDeDxYaxis(100, 0, 10000);
}

void MpdV0StandardDaughterMonitor::Init()
{
   if (fInit == kTRUE) return;
   fInit = kTRUE;
   MakeHistogram1d("DCA", "DCA [cm]", 0);
   MakeHistogram2d("TpcDeDx", "p [GeV/c]", "dEdX [AU]", 0);
}

MpdV0StandardDaughterMonitor::~MpdV0StandardDaughterMonitor() {}

void MpdV0StandardDaughterMonitor::FillDstTrack(const MpdTrack &track, Bool_t status)
{
   TVector3 dca(track.GetDCAX(), track.GetDCAY(), track.GetDCAZ());
   TVector3 mom(track.GetPx(), track.GetPy(), track.GetPz());
   Fill1D(0, dca.Mag(), status);
   Fill2D(0, mom.Mag(), track.GetdEdXTPC(), status);
}

void MpdV0StandardDaughterMonitor::FillMiniDstTrack(const MpdMiniTrack &track, Bool_t status)
{
   TVector3 vertex = GetMiniEvent()->primaryVertex();
   TVector3 dca(track.gDCA(vertex));
   Fill1D(0, dca.Mag(), status);
   Fill2D(0, track.gPtot(), track.dEdx(), status);
}

MpdV0DaughterMonitor *MpdV0StandardDaughterMonitor::MakeCopy() const
{
   return new MpdV0StandardDaughterMonitor(*this);
}
