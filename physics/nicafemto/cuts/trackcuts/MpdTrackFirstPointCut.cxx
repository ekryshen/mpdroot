/*
 * MpdTrackFirstPointCut.cxx
 *
 *  Created on: 21 cze 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdTrackFirstPointCut.h"
#include "NicaMpdTrack.h"
#include "MpdTrack.h"
MpdTrackFirstPointCut::MpdTrackFirstPointCut() : NicaTrackCut(3){
	SetUnitName("FirstPointX[cm]",0);
	SetUnitName("FirstPointY[cm]",1);
	SetUnitName("FirstPointZ[cm]",2);
}

Bool_t MpdTrackFirstPointCut::Pass(NicaTrack* track) {
	NicaMpdTrack *mpd = (NicaMpdTrack*)track;
	SetValue(mpd->GetFirstPoint()->X(),0);
	SetValue(mpd->GetFirstPoint()->Y(),1);
	SetValue(mpd->GetFirstPoint()->Z(),2);
	return Validate();
}

MpdTrackFirstPointCut::~MpdTrackFirstPointCut() {
}

