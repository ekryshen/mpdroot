/*
 * NicaMpdMiniDstFullEvent.h
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDMINIDSTFULLEVENT_H_
#define NICAMPDMINIDSTFULLEVENT_H_

#include <TString.h>

#include "NicaComplexEvent.h"

class NicaMpdMiniDstFullEvent : public NicaComplexEvent {
protected:
  NicaMpdMiniDstFullEvent(NicaEvent* re, NicaEvent* im);

public:
  NicaMpdMiniDstFullEvent();
  void OnlyPrimary();
  void OnlyGlobal();
  virtual void Update();
  virtual TString GetFormatName() const { return "NicaMpdMiniDstFullEvent"; };
  virtual ~NicaMpdMiniDstFullEvent();
  ClassDef(NicaMpdMiniDstFullEvent, 1)
};

#endif /* NICAMPDMINIDSTFULLEVENT_H_ */
