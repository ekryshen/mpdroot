/*
 * NicaMiniDstFullV0Event.h
 *
 *  Created on: 26 lut 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_FORMAT_COMBINED_NICAMINIDSTFULLV0EVENT_H_
#define MPDROOT_NICA_MPD_FORMAT_COMBINED_NICAMINIDSTFULLV0EVENT_H_

#include "NicaMpdMiniDstFullEvent.h"


class NicaMiniDstFullV0Event : public NicaMpdMiniDstFullEvent {
public:
  NicaMiniDstFullV0Event();
  virtual void Update();
  virtual TString GetFormatName() const { return "NicaMpdMiniDstFullV0Event"; };
  virtual ~NicaMiniDstFullV0Event();
  ClassDef(NicaMiniDstFullV0Event, 1)
};


#endif /* MPDROOT_NICA_MPD_FORMAT_COMBINED_NICAMINIDSTFULLV0EVENT_H_ */
