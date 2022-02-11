/*
 * MpdV0StandardDaughetrMonitor.h
 *
 *  Created on: 11 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0STANDARDDAUGHTERMONITOR_H_
#define MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0STANDARDDAUGHTERMONITOR_H_

#include <RtypesCore.h>
#include <TH1.h>
#include <TString.h>
#include <TVector3.h>
#include <vector>

#include "MpdV0DaughterMonitor.h"

class MpdMiniTrack;

class TH2D;

class MpdV0StandardDaughterMonitor : public MpdV0DaughterMonitor {

public:
   MpdV0StandardDaughterMonitor();
   MpdV0StandardDaughterMonitor(const MpdV0StandardDaughterMonitor &other) = default;
   void                  SetDCAAxis(Int_t bins, Double_t min, Double_t max) { SetXaxis1d(0, bins, min, max); };
   void                  SetTpcDeDxXaxis(Int_t bins, Double_t min, Double_t max) { SetXaxis2d(0, bins, min, max); };
   void                  SetTpcDeDxYaxis(Int_t bins, Double_t min, Double_t max) { SetYaxis2d(0, bins, min, max); };
   virtual void          Init();
   virtual void          FillDstTrack(const MpdTrack &track, Bool_t status);
   virtual void          FillMiniDstTrack(const MpdMiniTrack &track, Bool_t status);
   MpdV0DaughterMonitor *MakeCopy() const;
   virtual ~MpdV0StandardDaughterMonitor();
   ClassDef(MpdV0StandardDaughterMonitor, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0STANDARDDAUGHTERMONITOR_H_ */
