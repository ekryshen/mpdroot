/*
 * NicaTrackIsGlobalCut.h
 *
 *  Created on: 22 gru 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDTRACKISGLOBALCUT_H_
#define MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDTRACKISGLOBALCUT_H_

#include "NicaTrackCut.h"

class MpdTrackIsGlobalCut : public NicaTrackCut {
public:
  MpdTrackIsGlobalCut();
  void AcceptAll() { SetMinMax(0, 1); };
  void AcceptOnlyGlobals() { SetMinAndMax(1); };
  void AcceptOnlyNonGlobal() { SetMinAndMax(0); };
  virtual Bool_t Init(Int_t task_id);
  virtual Bool_t Pass(NicaTrack* track);
  virtual ~MpdTrackIsGlobalCut();
  ClassDef(MpdTrackIsGlobalCut, 1)
};


#endif /* MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDTRACKISGLOBALCUT_H_ */
