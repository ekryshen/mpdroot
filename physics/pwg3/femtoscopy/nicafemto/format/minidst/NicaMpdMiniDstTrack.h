/*
 * NicaMpdMiniDstTrack.h
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_FORMAT_MINIDST_NICAMPDMINIDSTTRACK_H_
#define MPDROOT_NICA_MPD_FORMAT_MINIDST_NICAMPDMINIDSTTRACK_H_

#include "NicaMpdMiniDstEvent.h"
#include "NicaMpdTrack.h"
#include "NicaTrackInterface.h"

class MpdMiniTrack;

class NicaMpdMiniDstTrack : public NicaMpdTrack {

public:
  NicaMpdMiniDstTrack();
  virtual void Update(MpdMiniTrack* track, NicaMpdMiniDstEvent::eMode mode);
  virtual ~NicaMpdMiniDstTrack();
  ClassDef(NicaMpdMiniDstTrack, 1)
};

class NicaMpdMiniDstTrackInterface : public NicaTrackInterface {
public:
  NicaMpdMiniDstTrackInterface();
  virtual ~NicaMpdMiniDstTrackInterface();
  ClassDef(NicaMpdMiniDstTrackInterface, 1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_MINIDST_NICAMPDMINIDSTTRACK_H_ */
