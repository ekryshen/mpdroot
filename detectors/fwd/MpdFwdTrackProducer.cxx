//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdHitProducer
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#include "TClonesArray.h"
#include "FairRootManager.h"
#include "MpdFwdHit.h"
#include "MpdFwdTrack.h"
#include "MpdFwdPoint.h"
#include "MpdMCTrack.h"

#include "MpdFwdTrackProducer.h"
#include "map"
#include "unordered_set"

ClassImp(MpdFwdTrackProducer)

MpdFwdTrackProducer::MpdFwdTrackProducer() : FairTask("FwdTrackProducer"){}
MpdFwdTrackProducer::~MpdFwdTrackProducer(){}

InitStatus MpdFwdTrackProducer::Init(){
  fMcTracks = (TClonesArray*) FairRootManager::Instance()->GetObject("MCTrack");
  fFwdPoints = (TClonesArray*) FairRootManager::Instance()->GetObject("FwdPoint");
  fFwdHits = (TClonesArray*) FairRootManager::Instance()->GetObject("FwdHit");
  fFwdTracks = new TClonesArray("MpdFwdTrack");
  FairRootManager::Instance()->Register("FwdTrack", "Fwd", fFwdTracks, kTRUE);
  return kSUCCESS;
}

void MpdFwdTrackProducer::Exec(Option_t *opt){
  fFwdTracks->Clear();
  
//  std::map<int,MpdFwdTrack*> mapFwdTracks;
//  for (Int_t ihit=0; ihit<fFwdHits->GetEntriesFast(); ihit++){
//    MpdFwdHit* hit = (MpdFwdHit*) fFwdHits->UncheckedAt(ihit);
//    MpdFwdPoint* point = (MpdFwdPoint*) fFwdPoints->UncheckedAt(hit->GetRefIndex());
//    Int_t trackID = point->GetTrackID();
//    if (mapFwdTracks.find(trackID)==mapFwdTracks.end()){
//      MpdFwdTrack* tr = new ((*fFwdTracks)[fFwdTracks->GetEntriesFast()]) MpdFwdTrack();
//      mapFwdTracks[trackID] = tr;
//    }
//    MpdFwdTrack* track = mapFwdTracks[trackID];
//    track->AddHitIndex(ihit);
//  }

  std::unordered_map<Int_t, Int_t> mapFwdTracks;
  Int_t itrk = 0;
  for (Int_t ihit=0; ihit<fFwdHits->GetEntriesFast(); ihit++){
    MpdFwdHit* hit = (MpdFwdHit*) fFwdHits->UncheckedAt(ihit);
    MpdFwdPoint* point = (MpdFwdPoint*) fFwdPoints->UncheckedAt(hit->GetRefIndex());
    Int_t trackID = point->GetTrackID();
    MpdMCTrack* mcTrack = (MpdMCTrack*) fMcTracks->UncheckedAt(trackID);
    if (mcTrack->GetMotherId()>=0) continue;
    Double_t ptMC = mcTrack->GetPt();
    Double_t pdgCode = mcTrack->GetPdgCode();
    Double_t pMC = mcTrack->GetP();
    MpdFwdTrack* tr = nullptr;
    if (mapFwdTracks.find(trackID)==mapFwdTracks.end()){
      tr = new ((*fFwdTracks)[fFwdTracks->GetEntriesFast()]) MpdFwdTrack();
      tr->SetPtMC(ptMC);
      tr->SetPMC(pMC);
      tr->SetPdgCode(pdgCode);
      mapFwdTracks[trackID] = itrk;
      itrk++;
    } else {
      tr = (MpdFwdTrack*) fFwdTracks->UncheckedAt(mapFwdTracks.at(trackID));
    }
    tr->AddHitIndex(ihit);
  }

  printf("Event\n");
  for (Int_t itr=0; itr<fFwdTracks->GetEntriesFast(); itr++){
    MpdFwdTrack* tr = (MpdFwdTrack*) fFwdTracks->UncheckedAt(itr);
    //tr->Print();
    for (Int_t ihit=0;ihit<tr->GetNIndices();ihit++){
      Int_t hitIndex = tr->GetHitIndex(ihit);
      MpdFwdHit* hit = (MpdFwdHit*) fFwdHits->UncheckedAt(hitIndex);
      printf("%.2f ",hit->GetZ());
    }
    printf("\n");
  }
  
}

void MpdFwdTrackProducer::Finish(){
}
