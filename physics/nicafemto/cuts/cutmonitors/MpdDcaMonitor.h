/*
 * MpdDCAMonitor.h
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDDCAMONITOR_H_
#define MPDDCAMONITOR_H_

#include "NicaPropertyMonitorXY.h"
class MpdDcaMonitor : public NicaPropertyMonitorXY {
 public:
  MpdDcaMonitor();
  virtual void Update(Bool_t passed, TObject *obj);
  NicaCutMonitor *MakeCopy() const { return new MpdDcaMonitor(*this); };
  virtual ~MpdDcaMonitor();
  ClassDef(MpdDcaMonitor, 1)
};

#endif /* MPDDCAMONITOR_H_ */
