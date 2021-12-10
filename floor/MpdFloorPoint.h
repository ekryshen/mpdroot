/*
 * MpdFloorPoint.h
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef Floor_Floor_MPDFloorPOINT_H_
#define Floor_Floor_MPDFloorPOINT_H_
#include "TObject.h"
#include "TLorentzVector.h"
#include "FairMCPoint.h"
class MpdFloorPoint : public FairMCPoint{
private:
public:
	MpdFloorPoint();
	MpdFloorPoint(Int_t trackId, Int_t detId, TVector3 pos, TVector3 mom,
			Double_t time, Double_t length, Double_t eloss,Int_t eventId = 0);
	virtual ~MpdFloorPoint();
	ClassDef(MpdFloorPoint,1)
};

#endif /* Floor_Floor_MPDFloorPOINT_H_ */
