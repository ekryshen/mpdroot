#ifndef MPD_V0_FINDER_H
#define MPD_V0_FINDER_H

#include "MpdEvent.h"
#include "MpdV0.h"
#include "MpdV0Cuts.h"
#include "MpdHelix.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdKalmanTrack.h"
#include "MpdVertex.h"
#include "MpdTrack.h"
#include "MpdParticle.h"
#include "MpdTofMatching.h"
#include "FairMCTrack.h"
//#include "MpdEtofMatching.h"
//#include "MpdEctKalmanTrack.h"
#include "MpdParticleIdentification.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include <Riostream.h>

#include <set>
#include <map>
#include <vector>
#include <iomanip>
#include <ctime>

using std::cout;
using std::endl;
using namespace MpdV0Cuts;

class MpdV0Finder : public FairTask {
public:
   MpdV0Finder(const char *name = "MpdV0Finder", Int_t iVerbose = 1);
   virtual ~MpdV0Finder();
   InitStatus Init();
   InitStatus ReInit();
   void       Exec(Option_t *option);
   void       Finish();
   void       Reset();
   Double_t   PCA(MpdKalmanTrack *tr, TVector3 point);
   // TLorentzVector CalcMass(MpdKalmanTrack *tr1, MpdKalmanTrack *tr2, Int_t k, Int_t l);  // evaluate mass of mother
   // particle
   Int_t GetHighestProb(MpdTrack *track);
   void  FillMpdTracks();
   void  SelectTracks();
   void  FindMCVertices();
   Int_t CheckMCTracks(MpdKalmanTrack *tr1, MpdKalmanTrack *tr2);
   //  Double_t DCAV0(MpdKalmanTrack *tr1, MpdKalmanTrack *tr2)  ;
   // Double_t DCAtracks(MpdKalmanTrack *tr1, MpdKalmanTrack *tr2);

private:
   // Int_t protonPDG;
   // Int_t pionPDG;  //pi
   // Int_t lambdaPDG;
   Int_t npr;
   Int_t npi;
   Int_t nvert;
   Int_t Nreco;
   Int_t Nok;
   Int_t Nmc;
   Int_t Nlam;
   Int_t Nok1;
   Int_t Nmc1;
   Int_t Nlam1;

   TVector3 PVert;
   TVector3 vtx;
   // vector<MpdParticle*> pi;
   // vector<MpdParticle*> pr;
   TClonesArray *pi;
   TClonesArray *pr;

   Int_t         fNMPDtracks;
   Int_t         fNtracks;
   TClonesArray *fV0Cont;
   TClonesArray *fPrimVertex;
   TClonesArray *fKFTracks;
   TClonesArray *fTofMatching;
   TClonesArray *fMpdTracks;
   TClonesArray *fMCTracks;
   MpdV0        *V0;
   MpdEvent     *fEvent;

   TClonesArray *fPrimVertex1;
   TClonesArray *fKFTracks1;
   TClonesArray *fMpdTracks1;
   // TClonesArray* fMpdEvent;
   // TClonesArray* fKFEctTracks;
   // TClonesArray* fEtofMatching;
   // TObjArray* fKalmanTrs;

   // Mother cuts
   ClassDef(MpdV0Finder, 1);
};
#endif
