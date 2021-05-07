/*
 * MpdPadsEdgeSectorCut.h
 *
 *  Created on: 7 sie 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDPADSEDGESECTORCUT_H_
#define MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDPADSEDGESECTORCUT_H_


#include "NicaTrackCut.h"

class NicaTpcSectorGeo;

class MpdPadsEdgeSectorCut : public NicaTrackCut{
	Double_t fEdge;
	NicaTpcSectorGeo *fSec;
public:
	MpdPadsEdgeSectorCut();
	MpdPadsEdgeSectorCut(const MpdPadsEdgeSectorCut &other);
	void SetEdge(Double_t edge){fEdge = edge;};
	Bool_t Init(Int_t format_id);
	virtual Bool_t Pass(NicaTrack *track);
	NicaPackage *Report()const;
	virtual ~MpdPadsEdgeSectorCut();
	ClassDef(MpdPadsEdgeSectorCut,1)
};

#endif /* MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDPADSEDGESECTORCUT_H_ */
