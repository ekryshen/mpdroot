#ifndef MPDPAIRKK_H
#define MPDPAIRKK_H

#include <deque>

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdHelix.h"
#include "FairMCEventHeader.h"
#include "MpdKalmanFilter.h"
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdPairKKTrack.h"
#include "MpdPairKKParams.h"

class MpdTpcKalmanTrack;

class MpdPairKK : public MpdAnalysisTask {

public:
   MpdPairKK() {}
   MpdPairKK(const char *name, const char *outputName = "taskName");
   ~MpdPairKK() {} // TODO: normal descructor with cleaning off histos

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

   int centrality(int tpcMultiplicity) const;

   long int IsSameParent(long int prim1, long int prim2) const;

   float dEdx_sigma_K(float dEdx, float mom) const;

   float Beta_sigma_K(float beta, float mom) const;

private:
   // Event properties
   bool     isInitialized = false;
   bool     isMC          = true;
   int      mCenBin       = 0;
   int      mZvtxBin      = 0;
   TVector3 mPrimaryVertex;

   std::string     mParamConfig;
   MpdPairKKParams mParams;

   std::string mOutFile = "histos.root";

   MpdKalmanFilter      *mKF              = nullptr;
   MpdPid               *mPID             = nullptr;
   TClonesArray         *mMCTracks        = nullptr;
   TObjArray            *mEMCClusters     = nullptr;
   TClonesArray         *mKalmanTracks    = nullptr;
   TClonesArray         *mMpdGlobalTracks = nullptr;
   vector<MpdParticle *> mPartK;
   MpdKalmanHit          mKHit;
   TClonesArray         *eventM = nullptr; //(V)

   std::vector<MpdPairKKTrack> mP2; //(V) Negative tracks
   std::vector<MpdPairKKTrack> mP1; //(V) Positive tracks

   static constexpr short nCenBinsAna   = 7;                          //(V) number of bins in centrality for analysis
   static constexpr short nMixEventZ    = 10;                         //(V) number of bins in z direction for mixing
   static constexpr short nMixEventCent = 10;                         //(V) number of bins in centrality for mixing
   static constexpr short nMixTot       = nMixEventZ * nMixEventCent; //(V)

   const int nMixed = 10; //(V) Depth of mixing
   int       mixBin;
   int       anaBin;

   TList *mixedEvents[nMixTot];

   // Histograms
   TList                  mHistoList;
   static constexpr short mHistoCentBins = 5;
   // General QA
   TH1F *mhEvents     = nullptr;
   TH1F *mhVertex     = nullptr;
   TH1F *mhCentrality = nullptr;

   TH1F *mInvGen;

   TH2F *mInvNoPID;
   TH2F *mInvOnePID;
   TH2F *mInvTwoPID;

   TH2F *mInvMixNoPID;
   TH2F *mInvMixOnePID;
   TH2F *mInvMixTwoPID;

   TH2F *mInvTrueNoPID;
   TH2F *mInvTrueNoPIDPhi;

   TH2F *mInvTrueOnePID;
   TH2F *mInvTrueOnePIDPhi;

   TH2F *mInvTrueTwoPID;
   TH2F *mInvTrueTwoPIDPhi;

   TH1F *mInvGenBin[nCenBinsAna];

   TH2F *mInvTrueNoPIDPhiBin[nCenBinsAna];
   TH2F *mInvTrueOnePIDPhiBin[nCenBinsAna];
   TH2F *mInvTrueTwoPIDPhiBin[nCenBinsAna];

   TH2F *mInvNoPIDBin[nCenBinsAna];
   TH2F *mInvOnePIDBin[nCenBinsAna];
   TH2F *mInvTwoPIDBin[nCenBinsAna];

   TH2F *mInvMixNoPIDBin[nCenBinsAna];
   TH2F *mInvMixOnePIDBin[nCenBinsAna];
   TH2F *mInvMixTwoPIDBin[nCenBinsAna];

   ClassDef(MpdPairKK, 1);
};
#endif
