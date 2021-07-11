/*
 * MpdTrackChargeCut.h
 *
 *  Created on: 13 kwi 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDTRACKCHARGECUT_H_
#define INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDTRACKCHARGECUT_H_

#include "NicaTrackCut.h"

class MpdTrackChargeCut : public NicaTrackCut{
public:
	MpdTrackChargeCut();
	Bool_t Pass(NicaTrack *tr);
	virtual ~MpdTrackChargeCut();
	ClassDef(MpdTrackChargeCut,1)
};

#endif /* INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDTRACKCHARGECUT_H_ */
