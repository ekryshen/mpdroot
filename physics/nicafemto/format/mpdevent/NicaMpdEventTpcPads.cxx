/*
 * NicaMpdEventTpcPads.cxx
 *
 *  Created on: 2 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdEventTpcPads.h"
#include "NicaMpdTrackTpcPads.h"

NicaMpdEventTpcPads::NicaMpdEventTpcPads() :NicaMpdEvent("NicaMpdTrackTpcPads"){

}

NicaMpdEventTpcPads::~NicaMpdEventTpcPads() {
	// TODO Auto-generated destructor stub
}

void NicaMpdEventTpcPads::ShallowCopyTracks(NicaEvent* event) {
	fTracks->Clear();
	fTotalTracksNo = event->GetTotalTrackNo();
	fTracks->ExpandCreateFast(fTotalTracksNo);
	NicaMpdEventTpcPads *Event = (NicaMpdEventTpcPads*)event;
	for(int i=0;i<fTotalTracksNo;i++){
		NicaMpdTrackTpcPads *from = (NicaMpdTrackTpcPads*)Event->fTracks->UncheckedAt(i);
		NicaMpdTrackTpcPads *to = (NicaMpdTrackTpcPads*)fTracks->UncheckedAt(i);
		*to = *from;
	}
}
