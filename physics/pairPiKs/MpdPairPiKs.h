#ifndef MPDPAIRPIKS_H
#define MPDPAIRPIKS_H

#include <deque>

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdHelix.h"
#include "FairMCEventHeader.h"
#include "MpdKalmanFilter.h"
#include "TpcPoint.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdPairPiKsTrack.h"
#include "MpdPairPiKsParams.h"
#include "MpdTofMatching.h"
#include "MpdTofMatchingData.h"

class MpdTpcKalmanTrack;

class MpdPairPiKs : public MpdAnalysisTask {

public:
   MpdPairPiKs() {}
   MpdPairPiKs(const char *name, const char *outputName = "taskName");
   ~MpdPairPiKs() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

   void setOutFile(std::string filename = "histos.root") { mOutFile = filename; }

protected:
   bool selectEvent(MpdAnalysisEvent &event);

   void selectPionTrack(MpdAnalysisEvent &event);

   bool selectTrack(MpdTrack *tr);

   bool selectTrackKs(MpdTrack *tr);

   void selectKsTrack(MpdAnalysisEvent &event);

   void processHistograms(MpdAnalysisEvent &event);

   long int IsSameParent(long int prim1, long int prim2) const;

   MpdHelix MakeHelix(const MpdTpcKalmanTrack *tr);

private:
   float cen;

   // Event properties
   bool     isInitialized = false;
   bool     isMC          = true;
   int      mCenBin       = 0;
   int      mZvtxBin      = 0;
   int      mRPBin        = 0;
   TVector3 mPrimaryVertex;

   std::string       mParamConfig;
   MpdPairPiKsParams mParams;

   std::string mOutFile = "histos.root";

   BaseTpcSectorGeo   *secGeo  = nullptr;
   MpdTpcKalmanFilter *recoTpc = nullptr;
   //   MpdKalmanFilter      *mKF              = nullptr;
   TClonesArray *mMCTracks        = nullptr;
   TClonesArray *mKalmanTracks    = nullptr;
   TClonesArray *mMpdGlobalTracks = nullptr;
   TClonesArray *eventM           = nullptr; // (V)
                                             //   MpdKalmanHit          mKHit;
   vector<MpdParticle *> vPartK;

   std::vector<MpdPairPiKsTrack> mP2; // (V) Negative tracks
   std::vector<MpdPairPiKsTrack> mP1; // (V) Positive tracks

   static constexpr short nCenBinsAna = 7; // (V) number of bins in centrality for analysis

   static constexpr short nMixEventZ    = 10; // (V) number of bins in z direction for mixing
   static constexpr short nMixEventCent = 10; // (V) number of bins in centrality for mixing
   static constexpr short nMixEventRP   = 5;  // (V) number of bins in event plane for mixing
   static constexpr short nMixTot       = nMixEventZ * nMixEventCent * nMixEventRP; // (V)

   int nMixed = 10; // (V) Depth of mixing
   int mixBin;
   int anaBin;

   TList *mixedEvents[nMixTot];

   // Histograms
   TList mHistoList;

   // General QA
   TH1F *mhEvents     = nullptr;
   TH1F *mhVertex     = nullptr;
   TH1F *mhCentrality = nullptr;
   TH1F *mhEvPlane    = nullptr;

   TH1F *mInvGenKstp;
   TH1F *mInvGenKstm;

   TH2F *mInvKs;
   TH2F *mInvKs1;
   TH2F *mInvKs2;

   TH2F *mDcaPi;
   TH2F *mDcaPi1;

   TH2F *mInvTwoPIDPip;
   TH2F *mInvTwoPIDPim;

   TH2F *mInvMixTwoPIDPip;
   TH2F *mInvMixTwoPIDPim;

   TH2F *mInvTrueTwoPIDPip;
   TH2F *mInvTrueTwoPIDPim;
   TH2F *mInvTrueTwoPIDKstp;
   TH2F *mMassResTwoPIDKstp;
   TH2F *mInvTrueTwoPIDKstm;
   TH2F *mMassResTwoPIDKstm;

   TH1F *mInvGenKstpBin[nCenBinsAna];
   TH1F *mInvGenKstmBin[nCenBinsAna];

   TH2F *mInvTrueTwoPIDKstpBin[nCenBinsAna];
   TH2F *mInvTrueTwoPIDKstmBin[nCenBinsAna];

   TH2F *mInvTwoPIDPipBin[nCenBinsAna];
   TH2F *mInvTwoPIDPimBin[nCenBinsAna];

   TH2F *mInvMixTwoPIDPipBin[nCenBinsAna];
   TH2F *mInvMixTwoPIDPimBin[nCenBinsAna];

   TH2F *mAccEffGen;
   TH2F *mAccEffRecTwoPID;

   TFile        *inFileSim;
   TTree        *inTreeSim;
   TClonesArray *tpcPoints;

   ClassDef(MpdPairPiKs, 1);
};
#endif
