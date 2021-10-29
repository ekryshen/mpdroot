/*
 * NicaFairSource.h
 *
 *  Created on: 2 sie 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAFAIREVENTINTERFACE_H_
#define NICAFAIREVENTINTERFACE_H_

#include "FairMCEventHeader.h"
#include "NicaTrack.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TRegexp.h>
#include "NicaTrackClones.h"

#include "NicaMCEvent.h"
#include "NicaMCEventInterface.h"
#include "NicaMpdMcTrack.h"
/**
 * interface to FairMC data event.
 */
class NicaMpdMcEvent;
class NicaMpdMcEventInterface : public NicaMCEventInterface {
	friend class NicaMpdMcEvent;
	FairMCEventHeader *fEvent;
	NicaTrackClones *fMcracks;
protected:
	virtual void ConnectToTree();
public:
	NicaMpdMcEventInterface();
	virtual Int_t GetTotalTrackNo() const {return fMcracks->GetEntriesFast();};
	virtual void Register(Bool_t write);
	virtual void Clear(Option_t *opt="");
	virtual void Compress(Int_t *map, Int_t map_size);
	virtual void CopyData(NicaEventInterface *s);
	virtual void CopyAndCompress(NicaEventInterface *s, Int_t *map, Int_t map_size);
	virtual void FillTrackInterface(NicaTrackInterface *track, Int_t index);
	virtual TObject *GetRawEventPointer()const{return fEvent;};
	virtual TObject * GetRawTrackPointer(Int_t index)const{return fMcracks->UncheckedAt(index);}
	virtual NicaTrackInterface *GetTrackInterface()const;
	/** GETTERS SETTERS */
	virtual void SetVertex(Double_t x, Double_t y, Double_t z){fEvent->SetVertex(x,y,z);};
	virtual void SetPhi(Double_t phi, Double_t phi_error){};
	virtual void SetB(Double_t b){};
	virtual Double_t GetB()const{return 0;};
	virtual Double_t GetPhi()const{return 0;};
	virtual Double_t GetPhiError()const{return 0;};
	virtual TLorentzVector GetVertex()const;
	virtual ~NicaMpdMcEventInterface();
	ClassDef(NicaMpdMcEventInterface,1)
};

#endif /* NICAFAIREVENTINTERFACE_H_ */
