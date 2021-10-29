/*
 * MpdPIDOnTheFly.h
 *
 *  Created on: 19 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_FORMAT_MPDEVENT_MPDPIDONTHEFLY_H_
#define MPDROOT_NICA_MPD_FORMAT_MPDEVENT_MPDPIDONTHEFLY_H_
#include "FairTask.h"
#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdTrack.h"
#include "MpdMiniTrack.h"
#include "TVector3.h"
class MpdPIDOnTheFly : public FairTask{
	MpdEvent *fEvent;
	MpdPid *fPID;
	TVector3 *fEventVector, *fMCVector;
	TClonesArray *fMiniTrack;
    void FillTrackDCA(MpdTrack *track, TVector3 *recoVertex, TVector3 *mcVertex);
	void FillTrackPID(MpdTrack *track);
    void FillTrackPID(MpdMiniTrack *track);
public:
	MpdPIDOnTheFly();
    virtual InitStatus Init();
    virtual void Exec(Option_t *opt);
	virtual ~MpdPIDOnTheFly();
	ClassDef(MpdPIDOnTheFly,1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_MPDEVENT_MPDPIDONTHEFLY_H_ */
