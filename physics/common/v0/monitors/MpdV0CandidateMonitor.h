/*
 * MpdV0CandidateMonitor.h
 *
 *  Created on: 11 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0CANDIDATEMONITOR_H_
#define MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0CANDIDATEMONITOR_H_

#include "MpdV0Monitor.h"

class MpdV0Particle;

class MpdV0CandidateMonitor : public MpdV0Monitor {
protected:
   MpdV0CandidateMonitor(Int_t d1, Int_t d2) : MpdV0Monitor(d1, d2){};

public:
   MpdV0CandidateMonitor();
   virtual void                   Fill(const MpdV0Particle &particle, Bool_t status) = 0;
   virtual MpdV0CandidateMonitor *MakeCopy() const                                   = 0;
   virtual ~MpdV0CandidateMonitor();
   ClassDef(MpdV0CandidateMonitor, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0CANDIDATEMONITOR_H_ */
