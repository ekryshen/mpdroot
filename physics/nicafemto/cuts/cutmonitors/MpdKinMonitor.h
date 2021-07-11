/*
 * MpdKinMonitor.h
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDKINMONITOR_H_
#define MPDKINMONITOR_H_
#include "NicaPropertyMonitorXY.h"
class MpdKinMonitor : public NicaPropertyMonitorXY {
 public:
  MpdKinMonitor();
  virtual void Update(Bool_t passed, TObject *obj);
  NicaCutMonitor *MakeCopy() const { return new MpdKinMonitor(*this); };
  virtual ~MpdKinMonitor();
  ClassDef(MpdKinMonitor, 1)
};

#endif /* MPDROOT_NICA_MPD_CUTMONITORS_MPDKINMONITOR_H_ */
