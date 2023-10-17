//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdPoint
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------

#include "MpdFwdTrack.h"
#include "vector"

ClassImp(MpdFwdTrack)
MpdFwdTrack::MpdFwdTrack():TObject(),fHitIndices(),fPtMC(), fPMC(), fPdgCode() {
}

MpdFwdTrack::~MpdFwdTrack() {}

void MpdFwdTrack::AddHitIndex(Int_t hitIndex){
  fHitIndices.push_back(hitIndex);
}

void MpdFwdTrack::Print(){
  printf("Track test: ");
  for (auto &hit : fHitIndices){
    printf("%d ",hit);
  }

  printf("\n");
}
