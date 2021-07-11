/*
 * NicaMpdEventInterface.cxx
 *
 *  Created on: 3 sie 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdEventInterface.h"
#include "NicaMpdTrackInterface.h"
#include "FairRootManager.h"

NicaMpdEventInterface::NicaMpdEventInterface() :fEvent(NULL),fMode(kAllTracks){
	fEvent = new MpdEvent();
}

NicaMpdEventInterface::~NicaMpdEventInterface() {
	if(CanDeleteEvent())
	if(fEvent) delete fEvent;
}

void NicaMpdEventInterface::Compress(Int_t* map, Int_t map_size) {
	TClonesArray *prim= fEvent->GetPrimaryTracks();
	TClonesArray *glob = fEvent->GetGlobalTracks();
	switch(fMode){
		case kAllTracks:{
			Int_t glob_count = 0;
			Int_t prim_count = 0;
			Int_t shift = prim->GetEntriesFast();
			Int_t track_pos=0;
			for(int i=0;i<map_size;i++){
				Int_t good_track = map[i];
				for(int j=track_pos;j<good_track;j++){
					if(map[i]<shift)
						prim->RemoveAt(j);
					else
						glob->RemoveAt(j-shift);
				}
				track_pos = good_track+1;
			}
			prim->Compress();
			glob->Compress();
		}break;
		case kPrimaryTracks:{
			Int_t track_pos=0;
			for(int i=0;i<map_size;i++){
				Int_t good_track = map[i];
				for(int j=track_pos;j<good_track;j++){
					prim->RemoveAt(j);
				}
				track_pos = good_track+1;
			}
			prim->Compress();
		}break;
		case kGlobalTracks:{
			Int_t track_pos=0;
			for(int i=0;i<map_size;i++){
				Int_t good_track = map[i];
				for(int j=track_pos;j<good_track;j++){
					glob->RemoveAt(j);
				}
				track_pos = good_track+1;
			}
			glob->Compress();
		}break;
		}
}

void NicaMpdEventInterface::CopyAndCompress(NicaEventInterface* s, Int_t* map,
		Int_t map_size) {
	TClonesArray *prim_to = fEvent->GetPrimaryTracks();
	TClonesArray *glob_to = fEvent->GetGlobalTracks();
	TClonesArray *prim_from = ((NicaMpdEventInterface*)s)->fEvent->GetPrimaryTracks();
	TClonesArray *glob_from = ((NicaMpdEventInterface*)s)->fEvent->GetGlobalTracks();
	prim_to->Clear();
	glob_to->Clear();
	int shift = prim_from->GetEntries();
	switch(fMode){
	case kAllTracks:{
		Int_t glob_count = 0;
		Int_t prim_count = 0;
		for(int i=0;i<map_size;i++){
			Int_t pos = map[i];
			if(pos<shift){
				MpdTrack *track_from = (MpdTrack*)prim_from->UncheckedAt(pos);
				MpdTrack *track_to = (MpdTrack*)prim_to->ConstructedAt(prim_count++);
				*track_to = *track_from;
			}else{
				pos = pos-shift;
				MpdTrack *track_from = (MpdTrack*)glob_from->UncheckedAt(pos);
				MpdTrack *track_to = (MpdTrack*)glob_to->ConstructedAt(glob_count++);
				*track_to = *track_from;
			}
		}
	}break;
	case kPrimaryTracks:{
		for(int i=0;i<map_size;i++){
			MpdTrack *track_from = (MpdTrack*)prim_from->UncheckedAt(map[i]);
			MpdTrack *track_to = (MpdTrack*)prim_to->ConstructedAt(i);
			*track_to = *track_from;
		}
	}break;
	case kGlobalTracks:{
		for(int i=0;i<map_size;i++){
			MpdTrack *track_from = (MpdTrack*)glob_from->UncheckedAt(map[i]);
			MpdTrack *track_to = (MpdTrack*)glob_to->ConstructedAt(i);
			*track_to = *track_from;
		}
	}break;
	}
}

void NicaMpdEventInterface::ConnectToTree() {
	FairRootManager *manager = FairRootManager::Instance();
	fEvent = (MpdEvent*)manager->GetObject("MPDEvent.");
}

void NicaMpdEventInterface::FillTrackInterface(NicaTrackInterface* track,
		Int_t index) {
	NicaMpdTrackInterface *tr = (NicaMpdTrackInterface*)track;
	switch(fMode){
	case kAllTracks:{
		if(index<fEvent->GetEventInfoNofPrimaryTracks()){
			tr->SetRawTrack( fEvent->GetPrimaryTracks()->UncheckedAt(index));
		}else{
			index = index -fEvent->GetEventInfoNofPrimaryTracks();
			tr->SetRawTrack(fEvent->GetGlobalTracks()->UncheckedAt(index));
		}
	}break;
	case kPrimaryTracks:
		tr->SetRawTrack( fEvent->GetPrimaryTracks()->UncheckedAt(index));
		break;
	case kGlobalTracks:
		tr->SetRawTrack(fEvent->GetGlobalTracks()->UncheckedAt(index));
		break;
	}
}

void NicaMpdEventInterface::Register(Bool_t write) {
	FairRootManager *manager = FairRootManager::Instance();
	TString branchname = "MPDEvent.";
	manager->Register(branchname,"MPDEvent.",fEvent,write);
}

void NicaMpdEventInterface::SetRunInfoId(Int_t i) {
	fEvent->SetRunInfoRunId(i);
}

void NicaMpdEventInterface::CopyData(NicaEventInterface* s) {
	*fEvent = *((NicaMpdEventInterface*)s)->fEvent;
}

void NicaMpdEventInterface::SetMagneticField(TVector3 mag) const {
	fEvent->SetRunInfoMagneticFieldZ(mag.Z());
}

Int_t NicaMpdEventInterface::GetTotalTrackNo() const {
	switch(fMode){
		case kAllTracks:{
			return fEvent->GetGlobalTracks()->GetEntriesFast()+fEvent->GetPrimaryTracks()->GetEntriesFast();
		}break;
		case kGlobalTracks:
			return fEvent->GetGlobalTracks()->GetEntriesFast();
		break;
		case kPrimaryTracks:
			return fEvent->GetPrimaryTracks()->GetEntriesFast();
		break;
	}
	return 0;
}

Int_t NicaMpdEventInterface::GetRunInfoId() const {
	return fEvent->GetRunInfoRunId();
}

TVector3 NicaMpdEventInterface::GetMagneticField() const {
	Double_t mz =  fEvent->GetRunInfoMagneticFieldZ();
	return TVector3(0,0,mz);
}

TLorentzVector NicaMpdEventInterface::GetVertexError() const {
	TLorentzVector vx;
	vx.SetXYZT(fEvent->GetPrimaryVerticesXerr(),fEvent->GetPrimaryVerticesYerr(),fEvent->GetPrimaryVerticesZerr(),0);
	return vx;
}

NicaTrackInterface* NicaMpdEventInterface::GetTrackInterface() const {
	return new NicaMpdTrackInterface();
}

TObject* NicaMpdEventInterface::GetRawTrackPointer(Int_t index) const {
	switch(fMode){
	case kAllTracks:{
		if(index<fEvent->GetEventInfoNofPrimaryTracks()){
			return fEvent->GetPrimaryTracks()->UncheckedAt(index);
		}else{
			index = index -fEvent->GetEventInfoNofPrimaryTracks();
			return fEvent->GetGlobalTracks()->UncheckedAt(index);
		}
	}break;
	case kPrimaryTracks:
		return fEvent->GetPrimaryTracks()->UncheckedAt(index);
		break;
	case kGlobalTracks:
		return fEvent->GetGlobalTracks()->UncheckedAt(index);
		break;
	default:
		return NULL;
		break;
	}
}

TLorentzVector NicaMpdEventInterface::GetVertex() const {
	TLorentzVector v;
	v.SetXYZT(fEvent->GetPrimaryVerticesX(), fEvent->GetPrimaryVerticesY(), fEvent->GetPrimaryVerticesZ(),0);
	return v;
}


