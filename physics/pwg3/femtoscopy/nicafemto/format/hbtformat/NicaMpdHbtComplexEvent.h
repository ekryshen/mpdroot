/*
 * NicaMpdHbtComplexEvent.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDHBTCOMPLEXEVENT_H_
#define NICAMPDHBTCOMPLEXEVENT_H_

#include "NicaComplexEvent.h"

/**
 * complex event that represent minidst + femto data
 */
class NicaMpdHbtComplexEvent : public NicaComplexEvent {
 public:
  NicaMpdHbtComplexEvent();
  void Update(){};
  virtual ~NicaMpdHbtComplexEvent();
  NicaMpdHbtComplexEvent(const NicaMpdHbtComplexEvent &other) = default;
  NicaMpdHbtComplexEvent &operator=(const NicaMpdHbtComplexEvent &other) =
      default;
  ClassDef(NicaMpdHbtComplexEvent, 1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_HBTFORMAT_NICAMPDHBTCOMPLEXEVENT_H_ */
