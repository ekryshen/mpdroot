#ifndef MPDPAIRPK_H
#define MPDPAIRPK_H

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
#include "MpdPairPKTrack.h"
#include "MpdPairPKParams.h"
#include "MpdTofMatching.h"
#include "MpdTofMatchingData.h"

class MpdTpcKalmanTrack;

class MpdPairPK : public MpdAnalysisTask {

public:
   MpdPairPK() {}
   MpdPairPK(const char *name, const char *outputName = "taskName");
   ~MpdPairPK() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

   void setOutFile(std::string filename = "histos.root") { mOutFile = filename; }

protected:
   bool selectEvent(MpdAnalysisEvent &event);

   void selectPosTrack(MpdAnalysisEvent &event);

   bool selectTrack(MpdTrack *tr);

   void selectNegTrack(MpdAnalysisEvent &event);

   void processHistograms(MpdAnalysisEvent &event);

   long int IsSameParent(long int prim1, long int prim2) const;

private:
   float cen;

   // Event properties
   bool     isInitialized = false;
   bool     isMC          = true;
   int      mCenBin       = 0;
   int      mZvtxBin      = 0;
   int      mRPBin        = 0;
   TVector3 mPrimaryVertex;

   std::string     mParamConfig;
   MpdPairPKParams mParams;

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
   vector<MpdParticle *> mPartK;
   MpdKalmanHit          mKHit;
   TClonesArray         *eventM = nullptr; // (V)

   std::vector<MpdPairPKTrack> mP2; // (V) Negative tracks
   std::vector<MpdPairPKTrack> mP1; // (V) Positive tracks

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

   TH1F *mInvGen;

   TH2F *mInvNoPID;
   TH2F *mInvOnePID;
   TH2F *mInvTwoPID;

   TH2F *mInvMixNoPID;
   TH2F *mInvMixOnePID;
   TH2F *mInvMixTwoPID;

   TH2F *mInvTrueNoPID;
   TH2F *mInvTrueNoPIDLam1520;
   TH2F *mMassResNoPIDLam1520;

   TH2F *mInvTrueOnePID;
   TH2F *mInvTrueOnePIDLam1520;
   TH2F *mMassResOnePIDLam1520;

   TH2F *mInvTrueTwoPID;
   TH2F *mInvTrueTwoPIDLam1520;
   TH2F *mMassResTwoPIDLam1520;

   TH1F *mInvGenBin[nCenBinsAna];

   TH2F *mInvTrueNoPIDLam1520Bin[nCenBinsAna];
   TH2F *mInvTrueOnePIDLam1520Bin[nCenBinsAna];
   TH2F *mInvTrueTwoPIDLam1520Bin[nCenBinsAna];

   TH2F *mInvNoPIDBin[nCenBinsAna];
   TH2F *mInvOnePIDBin[nCenBinsAna];
   TH2F *mInvTwoPIDBin[nCenBinsAna];

   TH2F *mInvMixNoPIDBin[nCenBinsAna];
   TH2F *mInvMixOnePIDBin[nCenBinsAna];
   TH2F *mInvMixTwoPIDBin[nCenBinsAna];

   TH2F *mAccEffGen;
   TH2F *mAccEffRecNoPID;
   TH2F *mAccEffRecOnePID;
   TH2F *mAccEffRecTwoPID;

   TFile        *inFileSim;
   TTree        *inTreeSim;
   TClonesArray *tpcPoints;

   ClassDef(MpdPairPK, 1);
};
#endif
