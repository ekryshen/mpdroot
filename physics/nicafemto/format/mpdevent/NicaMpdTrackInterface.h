/*
 * NicaMpdTrackInterface.h
 *
 *  Created on: 3 sie 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_FORMAT_NICAMPDTRACKINTERFACE_H_
#define INTERFACES_MPDROOT_FORMAT_NICAMPDTRACKINTERFACE_H_

#include "NicaTrackInterface.h"

class NicaMpdTrackInterface : public NicaTrackInterface{
public:
	NicaMpdTrackInterface();
	void SetCharge(Double_t charge);
	void SetID(Int_t id);
	void SetPxPyPzE(Double_t px, Double_t py, Double_t pz, Double_t e);
	void SetMotherIndex(Int_t index);
	void SetStatus(Int_t stat);
	void SetPrimary(Int_t prim);
	Double_t GetPx()const;
	Double_t GetPy()const;
	Double_t GetPz()const;
	Double_t GetE()const;
	Double_t GetCharge()const;
	Int_t GetID()const;
	Int_t GetMotherIndex()const;
	Int_t GetStatus()const;
	Bool_t IsPrimary();
	virtual ~NicaMpdTrackInterface();
	ClassDef(NicaMpdTrackInterface,1)
};

#endif /* INTERFACES_MPDROOT_FORMAT_NICAMPDTRACKINTERFACE_H_ */
