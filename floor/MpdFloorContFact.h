/*
 * MpdFloorContFact.h
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef Floor_Floor_MPDFloorCONTFACT_H_
#define Floor_Floor_MPDFloorCONTFACT_H_

#include "FairContFact.h"

class FairContainer;

class MpdFloorContFact : public FairContFact{
private:
	void setAllContainers();
public:
	MpdFloorContFact();
	virtual ~MpdFloorContFact();
    FairParSet*	createContainer(FairContainer*);

ClassDef(MpdFloorContFact,0)
};

#endif /* Floor_Floor_MPDFloorCONTFACT_H_ */
