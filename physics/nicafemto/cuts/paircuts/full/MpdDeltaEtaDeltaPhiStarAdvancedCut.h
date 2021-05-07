/*
 * MpdDeltaEtaDeltaPhiStarCut.h
 *
 *  Created on: 4 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_MPDDELTAETADELTAPHISTARCUT_H_
#define MPDROOT_NICA_MPD_CUTS_MPDDELTAETADELTAPHISTARCUT_H_

#include "MpdTpcPadsPairCut.h"

/**
 * calculate delta-eta delta-phi*
 */
namespace MpdPadsFormat {
class MpdDeltaEtaDeltaPhiStarAdvancedCut : public MpdTpcPadsPairCut {
 public:
  MpdDeltaEtaDeltaPhiStarAdvancedCut();
  virtual Bool_t Pass(NicaTwoTrack *pair);
  static Int_t DeltaPhiStar() { return 0; };
  static Int_t DeltaEta() { return 1; };
  static Int_t DeltaPhiStarMin() { return 2; };
  virtual ~MpdDeltaEtaDeltaPhiStarAdvancedCut();
  ClassDef(MpdDeltaEtaDeltaPhiStarAdvancedCut, 1)
};
}  // namespace MpdPadsFormat
#endif /* MPDROOT_NICA_MPD_CUTS_MPDDELTAETADELTAPHISTARCUT_H_ */
