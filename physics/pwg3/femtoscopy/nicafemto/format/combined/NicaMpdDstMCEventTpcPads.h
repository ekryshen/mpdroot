/*
 * NicaMpdDstMCEventTpcPads.h
 *
 *  Created on: 1 kwi 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_NICAMPDDSTMCEVENTTPCPADS_H_
#define INTERFACES_MPDROOT_NICAMPDDSTMCEVENTTPCPADS_H_

#include <TString.h>

#include "NicaComplexEvent.h"

class NicaMpdDstMCEventTpcPads : public NicaComplexEvent{
public:
	NicaMpdDstMCEventTpcPads();
	void OnlyPrimary();
	void OnlyGlobal();
	void Update();
	virtual TString GetFormatName()const{return "NicaMpdDstMCEventTpcPads";};
	virtual ~NicaMpdDstMCEventTpcPads();
	ClassDef(NicaMpdDstMCEventTpcPads,1)
};

#endif /* INTERFACES_MPDROOT_NICAMPDDSTMCEVENTTPCPADS_H_ */
