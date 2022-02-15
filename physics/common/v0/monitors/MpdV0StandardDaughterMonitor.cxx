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

#include "MpdTrack.h"
#include "MpdMiniBTofPidTraits.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"
#include "MpdV0Monitor.h"

MpdV0StandardDaughterMonitor::MpdV0StandardDaughterMonitor() : MpdV0DaughterMonitor(1, 2)
{
   SetDCAAxis(100, 0, 10);
   SetTpcDeDxXaxis(200, 0, 5);
   SetTpcDeDxYaxis(200, 0, 10000);
   SetTofM2Axis(1000, -1, 1);
   SetTofPAxis(100, 0, 5);
}

void MpdV0StandardDaughterMonitor::Init()
{
   if (fInit == kTRUE) return;
   fInit = kTRUE;
   MakeHistogram1d("DCA", "DCA [cm]", 0);
   MakeHistogram2d("TpcDeDx", "p [GeV/c]", "dEdX [AU]", 0);
   MakeHistogram2d("MTof", "p [GeV/c]", "m^{2}_{ToF}", 1);
}

MpdV0StandardDaughterMonitor::~MpdV0StandardDaughterMonitor() {}

void MpdV0StandardDaughterMonitor::FillDstTrack(const MpdTrack &track, Bool_t status)
{
   TVector3 dca(track.GetDCAX(), track.GetDCAY(), track.GetDCAZ());
   TVector3 mom(track.GetPx(), track.GetPy(), track.GetPz());
   Fill1D(0, dca.Mag(), status);
   Fill2D(0, mom.Mag(), track.GetdEdXTPC(), status);
   Fill2D(1, mom.Mag(), track.GetTofMass2(), status);
}

void MpdV0StandardDaughterMonitor::FillMiniDstTrack(const MpdMiniTrack &track, Bool_t status)
{
   TVector3 vertex = GetMiniEvent()->primaryVertex();
   TVector3 dca(track.gDCA(vertex));
   Fill1D(0, dca.Mag(), status);
   Fill2D(0, track.gPtot(), track.dEdx(), status);
   Int_t tof = track.bTofPidTraitsIndex();
   if (tof < 0) {
      Fill2D(1, track.gPtot(), -1E+4, status);
   } else {
      MpdMiniBTofPidTraits *tofInfo = GetTofPidTraits(tof);
      Fill2D(1, track.gPtot(), tofInfo->massSqr(), status);
   }
}

MpdV0DaughterMonitor *MpdV0StandardDaughterMonitor::MakeCopy() const
{
   return new MpdV0StandardDaughterMonitor(*this);
}
