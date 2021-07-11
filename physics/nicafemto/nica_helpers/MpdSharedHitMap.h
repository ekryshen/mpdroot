/*
 * MpdSharedHitMap.h
 *
 *  Created on: 23 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_NICA_HELPERS_MPDSHAREDHITMAP_H_
#define MPDROOT_NICA_MPD_NICA_HELPERS_MPDSHAREDHITMAP_H_
#include "FairTask.h"
#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdTrack.h"
#include "TVector3.h"


class MpdSharedHitMap : public FairTask{
	MpdEvent *fEvent;
	TClonesArray *fKalmanTracks;
	Int_t fShared;
	Int_t fSharedPads;
	Int_t fRealShared;
	ULong64_t CheckTracks(MpdTpcKalmanTrack *ktr1, MpdTpcKalmanTrack *ktr2);
public:
	MpdSharedHitMap();
    virtual InitStatus Init();
    virtual void Exec(Option_t *opt);
	virtual ~MpdSharedHitMap();
	ClassDef(MpdSharedHitMap,1)
};

#endif /* MPDROOT_NICA_MPD_NICA_HELPERS_MPDSHAREDHITMAP_H_ */
