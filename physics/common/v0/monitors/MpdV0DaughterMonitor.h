/*
 * MpdV0DaughterMonitor.h
 *
 *  Created on: 11 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0DAUGHTERMONITOR_H_
#define MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0DAUGHTERMONITOR_H_

#include "MpdV0Monitor.h"

#include <RtypesCore.h>
#include <TH1.h>
#include <TString.h>
#include <vector>

class MpdEvent;
class MpdMiniEvent;

class MpdMiniTrack;

class MpdTrack;

class MpdV0DaughterMonitor : public MpdV0Monitor {
protected:
   MpdV0DaughterMonitor(Int_t d1, Int_t d2) : MpdV0Monitor(d1, d2){};

public:
   MpdV0DaughterMonitor();
   MpdV0DaughterMonitor(const MpdV0DaughterMonitor &other)                                  = default;
   virtual void                  FillDstTrack(const MpdTrack &track, Bool_t status)         = 0;
   virtual void                  FillMiniDstTrack(const MpdMiniTrack &track, Bool_t status) = 0;
   virtual void                  Init()                                                     = 0;
   virtual MpdV0DaughterMonitor *MakeCopy() const                                           = 0;
   virtual ~MpdV0DaughterMonitor();
   ClassDef(MpdV0DaughterMonitor, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0DAUGHTERMONITOR_H_ */
