/*
 * MpdTofMonitor.h
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTMONITORS_MPDTOFMONITOR_H_
#define MPDROOT_NICA_MPD_CUTMONITORS_MPDTOFMONITOR_H_
#include "NicaPropertyMonitorXY.h"
class MpdTofMonitor : public NicaPropertyMonitorXY {
 public:
  MpdTofMonitor();
  virtual void Update(Bool_t passed, TObject *obj);
  void SetAxisP(Int_t bins, Double_t min, Double_t max) {
    SetXaxis(bins, min, max);
  }
  void SetAxisM2(Int_t bins, Double_t min, Double_t max) {
    SetYaxis(bins, min, max);
  }
  NicaCutMonitor *MakeCopy() const { return new MpdTofMonitor(*this); };
  virtual ~MpdTofMonitor();
  ClassDef(MpdTofMonitor, 1)
};

#endif /* MPDROOT_NICA_MPD_CUTMONITORS_MPDTOFMONITOR_H_ */
