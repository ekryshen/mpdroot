/*
 * MpdFloorPar.cxx
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdFloorGeoPar.h"


#include "FairParamList.h"              // for FairParamList

#include "TObjArray.h"

ClassImp(MpdFloorGeoPar)

MpdFloorGeoPar::~MpdFloorGeoPar() {
	// TODO Auto-generated destructor stub
}

MpdFloorGeoPar::MpdFloorGeoPar(const char* name, const char* title,
		const char* context):
		FairParGenericSet(name,title,context),
		    fGeoSensNodes(new TObjArray()),
		fGeoPassNodes(new TObjArray())
			{
}

void MpdFloorGeoPar::clear(void) {
	delete fGeoSensNodes;
	delete fGeoPassNodes;
}

void MpdFloorGeoPar::putParams(FairParamList* l) {
	  if (!l) { return; }
	  l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
	l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}

Bool_t MpdFloorGeoPar::getParams(FairParamList* l) {
	  if (!l) { return kFALSE; }
	  if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) { return kFALSE; }
	  if (!l->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) { return kFALSE; }
	return kTRUE;
}
