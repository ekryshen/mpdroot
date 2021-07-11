/*
 * MpdHitsOnTheFly.cxx
 *
 *  Created on: 25 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdHitsOnTheFly.h"

MpdHitsOnTheFly::MpdHitsOnTheFly() :
fEvent(NULL),fKalmans(NULL),fMpdTpcHits(NULL),
fHitMapSize(0),fHitMap(NULL)
{
}

InitStatus MpdHitsOnTheFly::Init() {
	FairRootManager *manager = FairRootManager::Instance();
	fEvent = (MpdEvent*)manager->GetObject("MPDEvent.");
	fKalmans = (TClonesArray*)manager->GetObject("TpcKalmanTrack");
	fMpdTpcHits = (TClonesArray*)manager->GetObject("TpcRecPoint");
	fHitMapSize = 1;
	fHitMap = new Int_t[fHitMapSize];
	return kSUCCESS;
}

void MpdHitsOnTheFly::Exec(Option_t* opt) {
	if(fMpdTpcHits->GetEntriesFast()>fHitMapSize){
		delete []fHitMap;
		std::cout<<fHitMapSize<<" -> "<<fMpdTpcHits->GetEntriesFast()<<std::endl;
		fHitMapSize = fMpdTpcHits->GetEntriesFast();
		fHitMap = new Int_t[fHitMapSize];
	}
	for(int i=0;i<fMpdTpcHits->GetEntriesFast();i++){
		fHitMap[i]=0;
	}
	for(int i=0;i<fKalmans->GetEntriesFast();i++){
		MpdTpcKalmanTrack *kalman = (MpdTpcKalmanTrack*) fKalmans->UncheckedAt(i);
		TObjArray *khits = kalman->GetTrHits();
		if(khits)
		for(int j=0;j<khits->GetEntriesFast();j++){
			MpdKalmanHit *hit = (MpdKalmanHit*)khits->UncheckedAt(j);
			Int_t id = hit->GetIndex();// index of MpdTpcHit
			fHitMap[id]++;
		}
	}
	TClonesArray *tracks = fEvent->GetGlobalTracks();
	for(int i=0;i<tracks->GetEntriesFast();i++){
		MpdTrack *track  = (MpdTrack*)tracks->UncheckedAt(i);
		ULong64_t layerMap, sharedMap;
		GetHitMaps(i,layerMap,sharedMap);
		track->SetLayerHitMap(layerMap);
		track->SetSharedHitMap(sharedMap);
	}

}

void MpdHitsOnTheFly::GetHitMaps(Int_t particle_index, ULong64_t& layerHit,
		ULong64_t& sharedHit) {
    MpdTpcKalmanTrack *kalman = (MpdTpcKalmanTrack*) fKalmans->UncheckedAt(particle_index);
    layerHit = 0;
    sharedHit = 0;
	TObjArray *khits = kalman->GetTrHits();
	if(khits)
	for(int j=0;j<khits->GetEntriesFast();j++){
		MpdKalmanHit *hit = (MpdKalmanHit*)khits->UncheckedAt(j);
		Int_t id = hit->GetIndex();// index of MpdTpcHit
		Int_t layer = hit->GetLayer();
		if(layer>=0){
			if(fHitMap[id]>1)
				SETBIT(sharedHit,layer);
			SETBIT(layerHit,layer);
    	}
	}
}

MpdHitsOnTheFly::~MpdHitsOnTheFly() {
	if(fHitMapSize>0)
		delete []fHitMap;
}

