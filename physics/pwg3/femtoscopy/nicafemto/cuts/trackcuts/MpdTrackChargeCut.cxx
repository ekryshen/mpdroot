/*
 * MpdTrackChargeCut.cxx
 *
 *  Created on: 13 kwi 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdTrackChargeCut.h"
#include "NicaMpdTrack.h"
MpdTrackChargeCut::MpdTrackChargeCut():NicaTrackCut(1) {
	SetUnitName("charge");

}

Bool_t MpdTrackChargeCut::Pass(NicaTrack *tr) {
	SetValue(tr->GetCharge());
	return Validate();
}

MpdTrackChargeCut::~MpdTrackChargeCut() {
	// TODO Auto-generated destructor stub
}

