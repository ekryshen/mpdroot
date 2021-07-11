/*
 * NicaMpdDstMCEvent.h
 *
 *  Created on: 1 kwi 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_NICAMPDDSTMCEVENT_H_
#define INTERFACES_MPDROOT_NICAMPDDSTMCEVENT_H_


#include <TString.h>

#include "NicaComplexEvent.h"

class NicaMpdDstMCEvent : public NicaComplexEvent{
public:
	NicaMpdDstMCEvent();
	void OnlyPrimary();
	void OnlyGlobal();
	void Update();
	virtual TString GetFormatName()const{return "NicaMpdDstMCEvent";};
	virtual ~NicaMpdDstMCEvent();
	ClassDef(NicaMpdDstMCEvent,1)
};

#endif /* INTERFACES_MPDROOT_NICAMPDDSTMCEVENT_H_ */
