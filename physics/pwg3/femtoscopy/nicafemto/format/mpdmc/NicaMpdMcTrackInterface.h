/*
 * NicaFairTrackInterface.h
 *
 *  Created on: 3 sie 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAFAIRTRACKINTERFACE_H_
#define NICAFAIRTRACKINTERFACE_H_

#include "NicaMCTrackInterface.h"

class NicaMpdMcTrackInterface : public NicaMCTrackInterface{
public:
	NicaMpdMcTrackInterface();
	virtual ~NicaMpdMcTrackInterface();
	ClassDef(NicaMpdMcTrackInterface,1)
};

#endif /* NICAFAIRTRACKINTERFACE_H_ */
