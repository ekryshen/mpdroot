#ifndef MPDPAIRPILAMBDA_H
#define MPDPAIRPILAMBDA_H

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
#include "MpdPairPiLambdaTrack.h"
#include "MpdPairPiLambdaParams.h"
#include "MpdTofMatching.h"
#include "MpdTofMatchingData.h"

class MpdTpcKalmanTrack;

class MpdPairPiLambda : public MpdAnalysisTask {

public:
   MpdPairPiLambda() {}
   MpdPairPiLambda(const char *name, const char *outputName = "taskName");
   ~MpdPairPiLambda() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

   void setOutFile(std::string filename = "histos.root") { mOutFile = filename; }

protected:
   bool selectEvent(MpdAnalysisEvent &event);

   void selectPionTrack(MpdAnalysisEvent &event);

   bool selectTrack(MpdTrack *tr);

   bool selectTrackLambda(MpdTrack *tr);

   void selectLambdaTrack(MpdAnalysisEvent &event);

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

   std::string           mParamConfig;
   MpdPairPiLambdaParams mParams;

   std::string mOutFile = "histos.root";

   BaseTpcSectorGeo     *secGeo           = nullptr;
   MpdTpcKalmanFilter   *recoTpc          = nullptr;
   MpdKalmanFilter      *mKF              = nullptr;
   MpdPid               *mPID             = nullptr;
   TClonesArray         *mMCTracks        = nullptr;
   TObjArray            *mEMCClusters     = nullptr;
   TClonesArray         *mKalmanTracks    = nullptr;
   TClonesArray         *mMpdGlobalTracks = nullptr;
   TClonesArray         *mpdTofMatching   = nullptr;
   TClonesArray         *eventM           = nullptr; // (V)
   vector<MpdParticle *> mPartK;
   MpdKalmanHit          mKHit;

   std::vector<MpdPairPiLambdaTrack> mP2; // (V) Negative tracks
   std::vector<MpdPairPiLambdaTrack> mP1; // (V) Positive tracks

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

   TH1F *mInvGenSigmap;
   TH1F *mInvGenSigmam;

   TH2F *mInvLambda;
   TH2F *mInvLambda1;
   TH2F *mInvLambda2;

   TH2F *mDcaPi;
   TH2F *mDcaPi1;
   TH2F *mDcaP;
   TH2F *mDcaP1;

   TH2F *mInvTwoPIDPip;
   TH2F *mInvTwoPIDPim;

   TH2F *mInvMixTwoPIDPip;
   TH2F *mInvMixTwoPIDPim;

   TH2F *mInvTrueTwoPIDPip;
   TH2F *mInvTrueTwoPIDPim;
   TH2F *mInvTrueTwoPIDSigmap;
   TH2F *mMassResTwoPIDSigmap;
   TH2F *mInvTrueTwoPIDSigmam;
   TH2F *mMassResTwoPIDSigmam;

   TH1F *mInvGenSigmapBin[nCenBinsAna];
   TH1F *mInvGenSigmamBin[nCenBinsAna];

   TH2F *mInvTrueTwoPIDSigmapBin[nCenBinsAna];
   TH2F *mInvTrueTwoPIDSigmamBin[nCenBinsAna];

   TH2F *mInvTwoPIDPipBin[nCenBinsAna];
   TH2F *mInvTwoPIDPimBin[nCenBinsAna];

   TH2F *mInvMixTwoPIDPipBin[nCenBinsAna];
   TH2F *mInvMixTwoPIDPimBin[nCenBinsAna];

   TH2F *mAccEffGen;
   TH2F *mAccEffRecTwoPID;

   TFile        *inFileSim;
   TTree        *inTreeSim;
   TClonesArray *tpcPoints;

   ClassDef(MpdPairPiLambda, 1);
};
#endif
