/*
 * NicaMpdMiniDstMcEventInterface.cxx
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "NicaMpdMiniDstMcEventInterface.h"

#include <TClonesArray.h>

#include "MpdMiniMcTrack.h"
#include "NicaEventInterface.h"
#include "NicaTrackInterface.h"
#include "NicaMpdMiniDstMcTrack.h"
#include "NicaTrackClones.h"

NicaMpdMiniDstMcEventInterface::NicaMpdMiniDstMcEventInterface() {
    fEvent = new NicaTrackClones("MpdMiniMcEvent","McEvent","mc");
    fTracks = new NicaTrackClones("MpdMiniMcTrack","McTrack","mc");
}

void NicaMpdMiniDstMcEventInterface::ConnectToTree() {
    if(CanDeleteEvent()){
        fTracks->DeleteClones();
        fEvent->DeleteClones();
    }
    fTracks->GetFromTree();
    fEvent->GetFromTree();
}

void NicaMpdMiniDstMcEventInterface::Register(Bool_t write) {
    fEvent->Register(write);
    fTracks->Register(write);
}

void NicaMpdMiniDstMcEventInterface::Clear(Option_t *opt) {
    fTracks->GetArray()->Clear();
}

void NicaMpdMiniDstMcEventInterface::Compress(Int_t *map, Int_t map_size) {
    fTracks->Compress(map, map_size);
}

void NicaMpdMiniDstMcEventInterface::CopyData(NicaEventInterface *s) {
    NicaMpdMiniDstMcEventInterface *interface = (NicaMpdMiniDstMcEventInterface*)s;
    fTracks->CopyFrom<MpdMiniMcTrack>(interface->fTracks->GetArray());
    fEvent->CopyFrom<MpdMiniMcEvent>(interface->fTracks->GetArray());
}

void NicaMpdMiniDstMcEventInterface::CopyAndCompress(NicaEventInterface *s,
        Int_t *map, Int_t map_size) {
    NicaMpdMiniDstMcEventInterface *interface = (NicaMpdMiniDstMcEventInterface*)s;
    fTracks->CopyCompress<MpdMiniMcTrack>(interface->fTracks->GetArray(), map, map_size);
    fEvent->CopyFrom<MpdMiniMcEvent>(interface->fTracks->GetArray());
}

void NicaMpdMiniDstMcEventInterface::FillTrackInterface(
        NicaTrackInterface *track, Int_t index) {
    track->SetRawTrack(fTracks->UncheckedAt(index));
}

NicaTrackInterface* NicaMpdMiniDstMcEventInterface::GetTrackInterface() const {
    return new NicaMpdMiniDstMcTrackInterface();
}

TLorentzVector NicaMpdMiniDstMcEventInterface::GetVertex() const {
    return TLorentzVector(GetEvent()->primaryVertexX(),GetEvent()->primaryVertexY(),
            GetEvent()->primaryVertexZ(),0);
}

NicaMpdMiniDstMcEventInterface::~NicaMpdMiniDstMcEventInterface() {
    // TODO Auto-generated destructor stub
}

