/*
 * NicaMpdDstEventV0.h
 *
 *  Created on: 25 lut 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_FORMAT_MINIDSTV0_NICAMPDMINIDSTEVENTV0_H_
#define MPDROOT_NICA_MPD_FORMAT_MINIDSTV0_NICAMPDMINIDSTEVENTV0_H_

#include "NicaMpdMiniDstEvent.h"

class NicaMpdMiniDstEventV0 : public NicaMpdMiniDstEvent {

protected:
  virtual void ShallowCopyEvent(NicaEvent* event);
  NicaMpdMiniDstEventV0(TString trackname);

public:
  NicaMpdMiniDstEventV0(NicaMpdMiniDstEvent::eMode = kGlobalTrack);
  virtual void CreateSource();
  virtual void Update();
  virtual Bool_t ExistInTree() const;
  virtual TString GetFormatName() const;
  virtual ~NicaMpdMiniDstEventV0();
  ClassDef(NicaMpdMiniDstEventV0, 1)
};


#endif /* MPDROOT_NICA_MPD_FORMAT_MINIDSTV0_NICAMPDMINIDSTEVENTV0_H_ */
