/*
 * MpdHitsEdgeSectorCut.h
 *
 *  Created on: 14 sie 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDHITSEDGESECTORCUT_H_
#define MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDHITSEDGESECTORCUT_H_


#include "NicaTrackCut.h"

class NicaTpcSectorGeo;

class MpdHitsEdgeSectorCut  : public NicaTrackCut{
	Double_t fEdge;
	NicaTpcSectorGeo *fSec;
public:
	MpdHitsEdgeSectorCut();
	MpdHitsEdgeSectorCut(const MpdHitsEdgeSectorCut &other);
	void SetEdge(Double_t edge){fEdge = edge;};
	Bool_t Init(Int_t format_id);
	virtual Bool_t Pass(NicaTrack *track);
	NicaPackage *Report()const;
	virtual ~MpdHitsEdgeSectorCut();
	ClassDef(MpdHitsEdgeSectorCut,1)
};

#endif /* MPDROOT_NICA_MPD_CUTS_TRACKCUTS_MPDHITSEDGESECTORCUT_H_ */
