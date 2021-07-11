/*
 * NicaMpdMiniDstEvent.h
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_FORMAT_MINIDST_NICAMPDMINIDSTEVENT_H_
#define MPDROOT_NICA_MPD_FORMAT_MINIDST_NICAMPDMINIDSTEVENT_H_

#include "NicaExpEvent.h"

class NicaMpdMiniDstEvent : public NicaExpEventHelix {
public:
  enum eMode { kGlobalTrack, kPrimaryTrack };

protected:
  eMode fMode;
  NicaMpdMiniDstEvent(TString trackname);

public:
  NicaMpdMiniDstEvent(eMode = kGlobalTrack);
  virtual void CreateSource();
  virtual void Update();
  void OnlyPrimary() { fMode = kPrimaryTrack; };
  void OnlyGlobal() { fMode = kGlobalTrack; };
  virtual Bool_t ExistInTree() const;
  virtual TString GetFormatName() const;
  virtual ~NicaMpdMiniDstEvent();
  ClassDef(NicaMpdMiniDstEvent, 1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_MINIDST_NICAMPDMINIDSTEVENT_H_ */
