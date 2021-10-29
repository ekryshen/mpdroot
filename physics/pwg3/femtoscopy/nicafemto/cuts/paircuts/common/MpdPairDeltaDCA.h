/*
 * MpdPairDeltaDCA.h
 *
 *  Created on: 9 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_COMMON_MPDPAIRDELTADCA_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_COMMON_MPDPAIRDELTADCA_H_

#include "NicaTwoTrackCut.h"
#include "NicaHelix.h"

class MpdPairDeltaDCA : public NicaTwoTrackCut{
public:
	MpdPairDeltaDCA();
	virtual Bool_t Pass(NicaTwoTrack *pair);
	static Int_t DeltaDCAxy(){return 0;};
	static Int_t DeltaDCAz(){return 1;}
	static Int_t DeltaDCA(){return 2;};
	virtual ~MpdPairDeltaDCA();
	ClassDef(MpdPairDeltaDCA,1)
};

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_COMMON_MPDPAIRDELTADCA_H_ */
