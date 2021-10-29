/*
 * NicaMpdHbtComplexEvent.cxx
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdHbtComplexEvent.h"
#include "NicaMpdHbtEvent.h"
#include "NicaMpdMiniDstMcEvent.h"

NicaMpdHbtComplexEvent::NicaMpdHbtComplexEvent()
    : NicaComplexEvent(new NicaMpdHbtEvent(), new NicaMpdMiniDstMcEvent()) {}

NicaMpdHbtComplexEvent::~NicaMpdHbtComplexEvent() {
  // TODO Auto-generated destructor stub
}
