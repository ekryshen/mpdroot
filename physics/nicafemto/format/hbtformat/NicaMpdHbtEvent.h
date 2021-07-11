/*
 * MpdHbtEvent.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDHBTEVENT_H_
#define NICAMPDHBTEVENT_H_

#include "NicaExpEvent.h"
class NicaPackage;
class NicaMpdHbtEvent : public NicaExpEvent {
 public:
  NicaMpdHbtEvent();
  NicaMpdHbtEvent(const NicaMpdHbtEvent &other);
  NicaMpdHbtEvent &operator=(const NicaMpdHbtEvent &other);
  virtual void CreateSource();
  virtual Bool_t IsCompatible(const NicaEvent *buffered) const;
  virtual NicaPackage *Report() const;
  virtual ~NicaMpdHbtEvent();
  ClassDef(NicaMpdHbtEvent, 1)
};

#endif /* NICAMPDHBTEVENT_H_ */
