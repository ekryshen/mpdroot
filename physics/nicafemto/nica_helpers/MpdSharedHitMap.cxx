/*
 * MpdSharedHitMap.cxx
 *
 *  Created on: 23 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdSharedHitMap.h"
#include "TObjArray.h"
#include "MpdTpcHit.h"
#include "TVector3.h"
#include "TArray.h"

MpdSharedHitMap::MpdSharedHitMap() {
	fEvent = NULL;
	fKalmanTracks = NULL;
	fShared = 0;
	fSharedPads =0;
	fRealShared =0;
}

InitStatus MpdSharedHitMap::Init() {
	FairRootManager *manager = FairRootManager::Instance();
	fEvent = (MpdEvent*)manager->GetObject("MPDEvent.");
	fKalmanTracks = (TClonesArray*)manager->GetObject("TpcKalmanTrack");
	return kSUCCESS;
}

void MpdSharedHitMap::Exec(Option_t* opt) {
	std::cout<<" CALC KALMAN MAPS"<<std::endl;
	fShared = 0;
	fSharedPads= 0;
	fRealShared =0;
	TClonesArray *tracks = fKalmanTracks;//->GetGlobalTracks();
	for(int i=0;i<tracks->GetEntriesFast();i++){
		fSharedPads = TMath::Max(fSharedPads,((MpdTpcKalmanTrack*)fKalmanTracks->UncheckedAt(i))->GetNofHits());
		for(int j=i+1;j<tracks->GetEntriesFast();j++){
			CheckTracks((MpdTpcKalmanTrack*)fKalmanTracks->UncheckedAt(i),
					(MpdTpcKalmanTrack*)fKalmanTracks->UncheckedAt(j)
					);
		}
	}
	std::cout<<"KALMAN MAPS CALCULATED\t"<<fShared<<" "<<fSharedPads<<" "<<fRealShared<<std::endl;
}

ULong64_t MpdSharedHitMap::CheckTracks(MpdTpcKalmanTrack *tr1, MpdTpcKalmanTrack *tr2) {
	TObjArray *hits1 = tr1->GetTrHits();
	TObjArray *hits2 = tr2->GetTrHits();
//	std::cout<<hits1->UncheckedAt(0)->ClassName()<<std::endl;
//	exit(0);
	for(int i=0;i<hits1->GetEntriesFast();i++){
		MpdKalmanHit *hit1 = (MpdKalmanHit*)hits1->UncheckedAt(i);
		for(int j=0;j<hits2->GetEntriesFast();j++){
			MpdKalmanHit *hit2 = (MpdKalmanHit*)hits2->UncheckedAt(j);
			if(hit1->GetDetectorID()==hit2->GetDetectorID()){
				if(hit1->GetSignal()==hit2->GetSignal()){
					hit1->SetLength(-1E+3);
					hit2->SetLength(-1E+3);
					fShared++;
				}else{
			//		std::cout<<"*******"<<std::endl;
					TArrayI *in1 = hit1->Index();
					TArrayI *in2 = hit2->Index();
					if(in1->GetSize()>0&&in2->GetSize()>0){
						if(in1->At(0)==in2->At(0)){
							fRealShared++;
						}
					//	std::cout<<in1->At(0)<<" "<<in2->At(0)<<std::endl;
					}
				}

			}
		}
	}
	return 0;
}

MpdSharedHitMap::~MpdSharedHitMap() {
	// TODO Auto-generated destructor stub
}

