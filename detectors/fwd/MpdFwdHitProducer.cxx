//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdHitProducer
/// \brief Ideal hit producer
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#include "TClonesArray.h"
#include "FairRootManager.h"
#include "MpdFwdHit.h"
#include "MpdFwdPoint.h"
#include "TRandom.h"
#include "MpdFwdHitProducer.h"

ClassImp(MpdFwdHitProducer)

MpdFwdHitProducer::MpdFwdHitProducer(Double_t sigma, Double_t timeRes) : FairTask("FwdHitProducer"),
fMcPoints(0),fFwdHits(0),fSigma(sigma),fTimeRes(timeRes)
{}

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
    hit->SetX(gRandom->Gaus(mcPoint->GetX(),fSigma));
    hit->SetY(gRandom->Gaus(mcPoint->GetY(),fSigma));
    hit->SetZ(mcPoint->GetZ());
    hit->SetDx(fSigma);
    hit->SetDy(fSigma);
    hit->SetDz(0);
    hit->SetTimeStamp(gRandom->Gaus(mcPoint->GetTime()*1e9,fTimeRes)); // from s to ns
    hit->SetTimeStampError(fTimeRes);
  }
}

void MpdFwdHitProducer::Finish(){
}
