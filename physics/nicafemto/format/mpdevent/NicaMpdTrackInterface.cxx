/*
 * NicaMpdTrackInterface.cxx
 *
 *  Created on: 3 sie 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdTrackInterface.h"
#include "MpdTrack.h"

NicaMpdTrackInterface::NicaMpdTrackInterface() {
	// TODO Auto-generated constructor stub

}

void NicaMpdTrackInterface::SetCharge(Double_t charge) {

}

void NicaMpdTrackInterface::SetID(Int_t id) {
	((MpdTrack*)fRawObject)->SetID(id);
}

void NicaMpdTrackInterface::SetPxPyPzE(Double_t px, Double_t py, Double_t pz,
		Double_t e) {
}

void NicaMpdTrackInterface::SetMotherIndex(Int_t index) {
}

void NicaMpdTrackInterface::SetStatus(Int_t stat) {
}

void NicaMpdTrackInterface::SetPrimary(Int_t prim) {
}

Double_t NicaMpdTrackInterface::GetPx() const {
	return ((MpdTrack*)fRawObject)->GetPx();
}

Double_t NicaMpdTrackInterface::GetPy() const {
	return ((MpdTrack*)fRawObject)->GetPy();
}

Double_t NicaMpdTrackInterface::GetPz() const {
	return ((MpdTrack*)fRawObject)->GetPz();
}

Double_t NicaMpdTrackInterface::GetE() const {
	return 0;
}

Double_t NicaMpdTrackInterface::GetCharge() const {
	return 0;
}

Int_t NicaMpdTrackInterface::GetID() const {
	return ((MpdTrack*)fRawObject)->GetID();
}

Int_t NicaMpdTrackInterface::GetMotherIndex() const {
	return 0;
}

Int_t NicaMpdTrackInterface::GetStatus() const {
	return 0;
}

Bool_t NicaMpdTrackInterface::IsPrimary() {
	return kTRUE;
}

NicaMpdTrackInterface::~NicaMpdTrackInterface() {
	// TODO Auto-generated destructor stub
}

