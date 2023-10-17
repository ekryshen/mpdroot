//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdHitProducer
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#include "TClonesArray.h"
#include "FairRootManager.h"
#include "MpdFwdHit.h"
#include "MpdFwdPoint.h"
#include "TRandom.h"
#include "MpdFwdHitProducer.h"

ClassImp(MpdFwdHitProducer)

MpdFwdHitProducer::MpdFwdHitProducer() : FairTask("FwdHitProducer"){}
MpdFwdHitProducer::~MpdFwdHitProducer(){}

InitStatus MpdFwdHitProducer::Init(){
  fMcPoints = (TClonesArray*) FairRootManager::Instance()->GetObject("FwdPoint");
  fFwdHits = new TClonesArray("MpdFwdHit");
  FairRootManager::Instance()->Register("FwdHit", "Fwd", fFwdHits, kTRUE);
  return kSUCCESS;
}

void MpdFwdHitProducer::Exec(Option_t *opt){
  fFwdHits->Clear();
  for (Int_t ip=0; ip<fMcPoints->GetEntriesFast(); ip++){
    MpdFwdPoint* mcPoint = (MpdFwdPoint*) fMcPoints->UncheckedAt(ip);
    if (mcPoint->GetTrackID()<0)continue;
    MpdFwdHit* hit = new ((*fFwdHits)[fFwdHits->GetEntriesFast()]) MpdFwdHit();
    hit->SetRefIndex(ip);
//    double sigma = 0.000000005; //  cm
    double sigma = 0.01; //  cm
    hit->SetX(gRandom->Gaus(mcPoint->GetX(),sigma));
    hit->SetY(gRandom->Gaus(mcPoint->GetY(),sigma));
    hit->SetZ(mcPoint->GetZ());
    hit->SetTime(mcPoint->GetTime());
    hit->SetLength(mcPoint->GetLength());
  }
}

void MpdFwdHitProducer::Finish(){
}
