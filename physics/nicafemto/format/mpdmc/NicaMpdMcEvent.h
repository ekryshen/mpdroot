/*
 * NicaFairEvent.h
 *
 *  Created on: 05-07-2014
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAFAIREVENT_H_
#define NICAFAIREVENT_H_
#include "FairMCEventHeader.h"
#include "MpdMCTrack.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TRegexp.h>

#include "NicaMCEvent.h"
#include "NicaMpdMcTrack.h"
#include "NicaTrack.h"
/**
 * basic  class for reading FairMCEventHeader and FairMCTracks
 */
class NicaMpdMcEvent : public NicaMCEvent{
protected:
	NicaMpdMcEvent(TString trackname);
public:
	/**
	 * default constructor
	 */
	NicaMpdMcEvent();
	NicaMpdMcEvent(const NicaMpdMcEvent &other);
	virtual void CreateSource();
	virtual void Update();
	virtual void Clear(Option_t *opt);
	virtual void Print();
	virtual Bool_t ExistInTree()const;
	virtual TString GetFormatName()const;
	virtual ~NicaMpdMcEvent();
	ClassDef(NicaMpdMcEvent,1)
};

#endif /* NICAFAIREVENT_H_ */
