/*
 * NicaMpdEvent.h
 *
 *  Created on: 28 mar 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_NICAMPDEVENT_H_
#define INTERFACES_MPDROOT_NICAMPDEVENT_H_
#include "NicaMpdEvent.h"

#include "MpdEvent.h"
#include "NicaExpEvent.h"
#include "NicaMpdTrack.h"
class NicaMpdEvent : public NicaExpEventHelix {
  enum kTrackType { kAllTracks, kPrimaryTracks, kGlobalTracks };
  kTrackType fMode;

protected:
  NicaMpdEvent(TString trackname);

public:
  NicaMpdEvent();
  NicaMpdEvent(const NicaMpdEvent& other);
  void CreateSource();
  void Update();
  void OnlyPrimary();
  void OnlyGlobal();
  virtual Bool_t ExistInTree() const;
  virtual TString GetFormatName() const;
  virtual ~NicaMpdEvent();
  ClassDef(NicaMpdEvent, 1)
};

#endif /* INTERFACES_MPDROOT_NICAMPDEVENT_H_ */
