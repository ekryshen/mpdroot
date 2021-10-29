/*
 * MpdCylinderTpcPairCut.h
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_TPC_CYLINDER_MPDCYLINDERTPCPAIRCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_TPC_CYLINDER_MPDCYLINDERTPCPAIRCUT_H_

#include "NicaTwoTrackCut.h"

class MpdCylinderTpcPairCut : public NicaTwoTrackCut {
public:
	MpdCylinderTpcPairCut(Int_t size);
	virtual ~MpdCylinderTpcPairCut();
	ClassDef(MpdCylinderTpcPairCut,1)
};

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_TPC_CYLINDER_MPDCYLINDERTPCPAIRCUT_H_ */
