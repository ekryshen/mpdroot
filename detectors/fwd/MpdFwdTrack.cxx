//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdTrack
/// \brief Reconstructed track
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------

#include "MpdFwdTrack.h"
#include "vector"

ClassImp(MpdFwdTrack)
MpdFwdTrack::MpdFwdTrack():TObject(),fHitIndices() {
}

MpdFwdTrack::~MpdFwdTrack() {}

void MpdFwdTrack::AddHitIndex(Int_t hitIndex){
  fHitIndices.push_back(hitIndex);
}

void MpdFwdTrack::Print(){
  printf("Track: (z,x,y,ty,yz,q/pt) = (%.2f,%.2f,%.2f,%.2f,%2.2f,%.2f); ",fV[0],fV[1],fV[2],fV[3],fV[4],fV[5]);
  printf("hit indices: ");
  for (auto &hit : fHitIndices){
    printf("%d ",hit);
  }
  printf("\n");
}
