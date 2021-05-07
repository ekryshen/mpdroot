/*
 * NicaFairSource.cxx
 *
 *  Created on: 2 sie 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdMcEventInterface.h"

#include "FairRootManager.h"
#include "MpdMCTrack.h"

#include "NicaMpdMcTrackInterface.h"

NicaMpdMcEventInterface::NicaMpdMcEventInterface():fEvent(NULL),fMcracks(NULL) {
	fEvent = new FairMCEventHeader();
	fMcracks = new NicaTrackClones("MpdMCTrack","MCTrack","MC");
}

void NicaMpdMcEventInterface::Clear(Option_t* opt) {
	fMcracks->GetArray()->Clear(opt);
}

void NicaMpdMcEventInterface::ConnectToTree() {
	FairRootManager *manager  = FairRootManager::Instance();
	if(CanDeleteEvent()){
		if(fEvent) delete fEvent;
		fMcracks->DeleteClones();
	}
	fEvent = (FairMCEventHeader*)manager->GetObject("MCEventHeader.");
	if(fEvent==nullptr)
		fEvent = (FairMCEventHeader*)manager->GetObject("EventHeader.");
	fMcracks->GetFromTree();
}

void NicaMpdMcEventInterface::Compress(Int_t* map, Int_t map_size) {
	fMcracks->Compress(map, map_size);
}

void NicaMpdMcEventInterface::CopyData(NicaEventInterface* s) {
	FairMCEventHeader *header = (FairMCEventHeader*)((NicaMpdMcEventInterface*)s)->fEvent;
	NicaTrackClones *tracks = (NicaTrackClones*)((NicaMpdMcEventInterface*)s)->fMcracks;
	fEvent->SetRunID(header->GetRunID());
	fEvent->SetEventID(header->GetEventID());
	fEvent->SetVertex(header->GetX(),header->GetY(),header->GetZ());
	fEvent->SetRotX(header->GetRotX());
	fEvent->SetRotY(header->GetRotY());
	fEvent->SetRotZ(header->GetRotZ());
	fEvent->SetB(header->GetB());
	fEvent->SetNPrim(header->GetNPrim());
	fMcracks->CopyFrom<MpdMCTrack>(tracks->GetArray());
}

void NicaMpdMcEventInterface::CopyAndCompress(NicaEventInterface* s, Int_t* map,
		Int_t map_size) {
	FairMCEventHeader *header = (FairMCEventHeader*)((NicaMpdMcEventInterface*)s)->fEvent;

	fEvent->SetRunID(header->GetRunID());
	fEvent->SetEventID(header->GetEventID());
	fEvent->SetVertex(header->GetX(),header->GetY(),header->GetZ());
	fEvent->SetRotX(header->GetRotX());
	fEvent->SetRotY(header->GetRotY());
	fEvent->SetRotZ(header->GetRotZ());
	fEvent->SetB(header->GetB());
	fEvent->SetNPrim(header->GetNPrim());
	fMcracks->CopyCompress<MpdMCTrack>(((NicaMpdMcEventInterface*)s)->fMcracks->GetArray(),
			map, map_size);
}

void NicaMpdMcEventInterface::Register(Bool_t write) {
	if(fEvent==NULL){
		fEvent = new FairMCEventHeader();
		fMcracks = new NicaTrackClones("MpdMCTrack","MCTrack","MC");
	}
	FairRootManager *manager = FairRootManager::Instance();
	TString event_name = "MCEventHeader.";
	TString track_name = "MCTrack.";
	manager->Register(event_name,"MCEventHeader",fEvent,write);
	fMcracks->Register(write);
}

NicaTrackInterface* NicaMpdMcEventInterface::GetTrackInterface() const {
	return new NicaMpdMcTrackInterface();
}

void NicaMpdMcEventInterface::FillTrackInterface(NicaTrackInterface* track,
		Int_t index) {
	track->SetRawTrack(fMcracks->UncheckedAt(index));
}

NicaMpdMcEventInterface::~NicaMpdMcEventInterface() {
	if(CanDeleteEvent()){
		if(fEvent) delete fEvent;
		if(fMcracks) delete fMcracks;
	}
}

TLorentzVector NicaMpdMcEventInterface::GetVertex() const {
	TLorentzVector vec;
	TVector3 vec3d;
	fEvent->GetVertex(vec3d);
	vec.SetXYZT(vec3d.X(),vec3d.Y(),vec3d.Z(),0);
	return vec;
}
