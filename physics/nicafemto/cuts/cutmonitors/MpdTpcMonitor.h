/*
 * MpdTpcMonitor.h
 *
 *  Created on: 20 sie 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTMONITORS_MPDTPCMONITOR_H_
#define MPDROOT_NICA_MPD_CUTMONITORS_MPDTPCMONITOR_H_

#include "NicaPropertyMonitorXY.h"
#include "NicaStd.h"

class MpdTpcMonitor : public NicaPropertyMonitorXY {
 public:
  MpdTpcMonitor();
  MpdTpcMonitor(const MpdTpcMonitor &other);
  virtual Bool_t Init(Int_t task_id);
  void SetAxisP(Int_t bins, Double_t min, Double_t max) {
    SetXaxis(bins, min, max);
  }
  void SetAxisDeDx(Int_t bins, Double_t min, Double_t max) {
    SetYaxis(bins, min, max);
  }
  virtual void Update(Bool_t passed, TObject *obj);
  NicaCutMonitor *MakeCopy() const { return new MpdTpcMonitor(*this); };
  virtual ~MpdTpcMonitor();
  ClassDef(MpdTpcMonitor, 1)
};

#endif /* MPDROOT_NICA_MPD_CUTMONITORS_MPDTPCMONITOR_H_ */
