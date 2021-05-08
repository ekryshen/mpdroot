/*
 * MpdBasicTrackCut.h
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDBASICTRACKCUT_H_
#define MPDBASICTRACKCUT_H_

#include "NicaCutsAndMonitors.h"
#include <TVector3.h>
class MpdKinMonitor;
class MpdTofMonitor;
class MpdTpcMonitor;
class NicaTrackDCACut;
class NicaTrackTpcToFCut;
class MpdDcaMonitor;
class NicaTrackEtaCut;
class NicaTrackPtCut;
class NicaTrackChargeCut;

class MpdBasicTrackCut : public NicaCutsAndMonitors {
  MpdTofMonitor* fToFMonitor;  //!
  MpdKinMonitor* fKinMonitor;  //!
  MpdTpcMonitor* fTpcMonitor;  //!
  MpdDcaMonitor* fDCAMonitor;  //!
protected:
  /**
   * create monitors
   * @param opt dca for dca monitor, kin for kinetic monitor
   * tpc for tpc monitor tof for tof monitor
   */
  virtual void AddAllCutMonitorRequests(Option_t* opt = "dca+kin+tpc+tof");

public:
  MpdBasicTrackCut();
  MpdTofMonitor* GetTofMonitor() const;
  MpdKinMonitor* GetKinMonitor() const;
  MpdTpcMonitor* GetTpcMonitor() const;
  MpdDcaMonitor* GetDCAMonitor() const;
  NicaTrackDCACut* GetDCACut() const;
  NicaTrackEtaCut* GetEtaCut() const;
  NicaTrackPtCut* GetPtCut() const;
  NicaTrackChargeCut* GetChargeCut() const;
  NicaTrackTpcToFCut* GetTpcTofCut() const;
  virtual ~MpdBasicTrackCut();
  ClassDef(MpdBasicTrackCut, 1)
};

#endif /* MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDBASICTRACKCUT_H_ */
