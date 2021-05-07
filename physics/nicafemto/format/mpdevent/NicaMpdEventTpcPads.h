/*
 * NicaMpdEventTpcPads.h
 *
 *  Created on: 2 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_FORMAT_MPDEVENT_NICAMPDEVENTTPCPADS_H_
#define MPDROOT_NICA_MPD_FORMAT_MPDEVENT_NICAMPDEVENTTPCPADS_H_
#include "NicaMpdEvent.h"

class NicaMpdEventTpcPads :public NicaMpdEvent{
protected:
	virtual void ShallowCopyTracks(NicaEvent *event);
public:
	NicaMpdEventTpcPads();
	virtual ~NicaMpdEventTpcPads();
	ClassDef(NicaMpdEventTpcPads,1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_MPDEVENT_NICAMPDEVENTTPCPADS_H_ */

