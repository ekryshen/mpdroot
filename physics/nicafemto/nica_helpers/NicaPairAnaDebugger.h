/*
 * NicaPairAnaDebugger.h
 *
 *  Created on: 21 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAPAIRANADEBUGGER_H_
#define NICAPAIRANADEBUGGER_H_
#include "FairTask.h"
#include "NicaTwoTrackAna.h"
class NicaPairAnaDebugger : public NicaTwoTrackAna {
 public:
  NicaPairAnaDebugger();
  void ProcessEvent();
  virtual ~NicaPairAnaDebugger();
  ClassDef(NicaPairAnaDebugger, 1)
};

#endif /* NICAPAIRANADEBUGGER_H_ */
