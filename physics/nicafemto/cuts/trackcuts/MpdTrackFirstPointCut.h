/*
 * MpdTrackFirstPointCut.h
 *
 *  Created on: 21 cze 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDTRACKFIRSTPOINTCUT_H_
#define INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDTRACKFIRSTPOINTCUT_H_

#include "NicaTrackCut.h"

class MpdTrackFirstPointCut : public NicaTrackCut{
public:
	MpdTrackFirstPointCut();
	Bool_t Pass(NicaTrack *track);
	virtual ~MpdTrackFirstPointCut();
	ClassDef(MpdTrackFirstPointCut,1)
};

#endif /* INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDTRACKFIRSTPOINTCUT_H_ */
