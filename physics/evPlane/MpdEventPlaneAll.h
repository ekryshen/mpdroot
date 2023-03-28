#ifndef MPDEVENTPLANEALL_H
#define MPDEVENTPLANEALL_H

#include <deque>

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdHelix.h"
#include "FairMCEventHeader.h"
#include "MpdKalmanFilter.h"
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdEventPlaneAllParams.h"
#include "TRandom.h"

class MpdTpcKalmanTrack;

class MpdEventPlaneAll : public MpdAnalysisTask {

public:
   MpdEventPlaneAll() {}
   MpdEventPlaneAll(const char *name, const char *outputName = "taskName");
   ~MpdEventPlaneAll();

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

   void setOutFile(std::string filename = "histos.root") { mOutFile = filename; }

protected:
   bool selectEvent(MpdAnalysisEvent &event);

   bool selectTrack(MpdTrack *tr);

   void processHistograms(MpdAnalysisEvent &event);

   float GetFHCalPhi(int iModule);

   int GetCentBin(float cent);

   std::map<int,float> SetZeroCorr();
   std::map<int,float> ReadEpCorrProfile(TProfile *const& prof);

private:
   // Event properties
   bool     isInitialized = false;
   TVector3 mPrimaryVertex;

   static constexpr short mFHCalModuleNum = 90;
   static constexpr short mFHCalMod1Side  = 45;
   static constexpr short nCentBins = 8; //from 0-10% up to 70-80%
   int mCorrStep; // 0 - raw, 1 - rec, 2 - shift

   std::string            mParamConfig;
   MpdEventPlaneAllParams mParams;

   std::string mOutFile = "histos.root";

   MpdKalmanFilter      *mKF              = nullptr;
   TClonesArray         *mMCTracks        = nullptr;
   TObjArray            *mEMCClusters     = nullptr;
   TClonesArray         *mKalmanTracks    = nullptr;
   TClonesArray         *mMpdGlobalTracks = nullptr;
   vector<MpdParticle *> mPartK;
   MpdKalmanHit          mKHit;

   //Maps to store info for the EP corrections
   // This is more optimal way to store information
   // in terms of RAM usage

   // For the recentering EP correction: <Qx>, <Qy> as a function of centrality bin
   std::map<int,float> mCorrQxFHCalFAll;
   std::map<int,float> mCorrQyFHCalFAll;
   std::map<int,float> mCorrQxFHCalNAll;
   std::map<int,float> mCorrQyFHCalNAll;
   std::map<int,float> mCorrQxFHCalSAll;
   std::map<int,float> mCorrQyFHCalSAll;
   std::map<int,float> mCorrQxTPCNAll;
   std::map<int,float> mCorrQyTPCNAll;
   std::map<int,float> mCorrQxTPCSAll;
   std::map<int,float> mCorrQyTPCSAll;

   // For the shift EP correction: <cosNphiEP>, <sinNphiEP> as a function of centrality bin
   std::map<int,float> mCorrCos1FHCalFAll;
   std::map<int,float> mCorrCos2FHCalFAll;
   std::map<int,float> mCorrCos3FHCalFAll;
   std::map<int,float> mCorrCos4FHCalFAll;
   std::map<int,float> mCorrCos5FHCalFAll;
   std::map<int,float> mCorrCos6FHCalFAll;
   std::map<int,float> mCorrCos7FHCalFAll;
   std::map<int,float> mCorrCos8FHCalFAll;
   std::map<int,float> mCorrSin1FHCalFAll;
   std::map<int,float> mCorrSin2FHCalFAll;
   std::map<int,float> mCorrSin3FHCalFAll;
   std::map<int,float> mCorrSin4FHCalFAll;
   std::map<int,float> mCorrSin5FHCalFAll;
   std::map<int,float> mCorrSin6FHCalFAll;
   std::map<int,float> mCorrSin7FHCalFAll;
   std::map<int,float> mCorrSin8FHCalFAll;

   std::map<int,float> mCorrCos1FHCalNAll;
   std::map<int,float> mCorrCos2FHCalNAll;
   std::map<int,float> mCorrCos3FHCalNAll;
   std::map<int,float> mCorrCos4FHCalNAll;
   std::map<int,float> mCorrCos5FHCalNAll;
   std::map<int,float> mCorrCos6FHCalNAll;
   std::map<int,float> mCorrCos7FHCalNAll;
   std::map<int,float> mCorrCos8FHCalNAll;
   std::map<int,float> mCorrSin1FHCalNAll;
   std::map<int,float> mCorrSin2FHCalNAll;
   std::map<int,float> mCorrSin3FHCalNAll;
   std::map<int,float> mCorrSin4FHCalNAll;
   std::map<int,float> mCorrSin5FHCalNAll;
   std::map<int,float> mCorrSin6FHCalNAll;
   std::map<int,float> mCorrSin7FHCalNAll;
   std::map<int,float> mCorrSin8FHCalNAll;

   std::map<int,float> mCorrCos1FHCalSAll;
   std::map<int,float> mCorrCos2FHCalSAll;
   std::map<int,float> mCorrCos3FHCalSAll;
   std::map<int,float> mCorrCos4FHCalSAll;
   std::map<int,float> mCorrCos5FHCalSAll;
   std::map<int,float> mCorrCos6FHCalSAll;
   std::map<int,float> mCorrCos7FHCalSAll;
   std::map<int,float> mCorrCos8FHCalSAll;
   std::map<int,float> mCorrSin1FHCalSAll;
   std::map<int,float> mCorrSin2FHCalSAll;
   std::map<int,float> mCorrSin3FHCalSAll;
   std::map<int,float> mCorrSin4FHCalSAll;
   std::map<int,float> mCorrSin5FHCalSAll;
   std::map<int,float> mCorrSin6FHCalSAll;
   std::map<int,float> mCorrSin7FHCalSAll;
   std::map<int,float> mCorrSin8FHCalSAll;

   std::map<int,float> mCorrCos1TPCNAll;
   std::map<int,float> mCorrCos2TPCNAll;
   std::map<int,float> mCorrCos3TPCNAll;
   std::map<int,float> mCorrCos4TPCNAll;
   std::map<int,float> mCorrSin1TPCNAll;
   std::map<int,float> mCorrSin2TPCNAll;
   std::map<int,float> mCorrSin3TPCNAll;
   std::map<int,float> mCorrSin4TPCNAll;

   std::map<int,float> mCorrCos1TPCSAll;
   std::map<int,float> mCorrCos2TPCSAll;
   std::map<int,float> mCorrCos3TPCSAll;
   std::map<int,float> mCorrCos4TPCSAll;
   std::map<int,float> mCorrSin1TPCSAll;
   std::map<int,float> mCorrSin2TPCSAll;
   std::map<int,float> mCorrSin3TPCSAll;
   std::map<int,float> mCorrSin4TPCSAll;

   // Histograms
   TList mHistoList;

   // General QA
   TH1F *mhEvents            = nullptr;
   TH1F *mhVertex            = nullptr;
   TH1F *mhHits              = nullptr;
   TH1F *mhEta               = nullptr;
   TH1F *mhDca               = nullptr;
   TH1F *mhPt                = nullptr;
   TH1F *mhQxRawFHCalFAll    = nullptr;
   TH1F *mhQyRawFHCalFAll    = nullptr;
   TH1F *mhPhiEPRawFHCalFAll = nullptr;
   TH1F *mhQxRawFHCalNAll    = nullptr;
   TH1F *mhQyRawFHCalNAll    = nullptr;
   TH1F *mhPhiEPRawFHCalNAll = nullptr;
   TH1F *mhQxRawFHCalSAll    = nullptr;
   TH1F *mhQyRawFHCalSAll    = nullptr;
   TH1F *mhPhiEPRawFHCalSAll = nullptr;
   TH1F *mhQxRawTPCNAll      = nullptr;
   TH1F *mhQyRawTPCNAll      = nullptr;
   TH1F *mhPhiEPRawTPCNAll   = nullptr;
   TH1F *mhQxRawTPCSAll      = nullptr;
   TH1F *mhQyRawTPCSAll      = nullptr;
   TH1F *mhPhiEPRawTPCSAll   = nullptr;
   TH1F *mhQxRecFHCalFAll    = nullptr;
   TH1F *mhQyRecFHCalFAll    = nullptr;
   TH1F *mhPhiEPRecFHCalFAll = nullptr;
   TH1F *mhQxRecFHCalNAll    = nullptr;
   TH1F *mhQyRecFHCalNAll    = nullptr;
   TH1F *mhPhiEPRecFHCalNAll = nullptr;
   TH1F *mhQxRecFHCalSAll    = nullptr;
   TH1F *mhQyRecFHCalSAll    = nullptr;
   TH1F *mhPhiEPRecFHCalSAll = nullptr;
   TH1F *mhQxRecTPCNAll      = nullptr;
   TH1F *mhQyRecTPCNAll      = nullptr;
   TH1F *mhPhiEPRecTPCNAll   = nullptr;
   TH1F *mhQxRecTPCSAll      = nullptr;
   TH1F *mhQyRecTPCSAll      = nullptr;
   TH1F *mhPhiEPRecTPCSAll   = nullptr;
   TH1F *mhPhiEPShfFHCalFAll = nullptr;
   TH1F *mhPhiEPShfFHCalNAll = nullptr;
   TH1F *mhPhiEPShfFHCalSAll = nullptr;
   TH1F *mhPhiEPShfTPCNAll   = nullptr;
   TH1F *mhPhiEPShfTPCSAll   = nullptr;
   TH1F *mhCorrStep          = nullptr;

   TProfile *mhCosFHCalNFHCalSAll = nullptr; 
   TProfile *mhCosTPCNTPCSAll     = nullptr;

   // For the recentering EP correction: <Qx>, <Qy>
   TProfile *mhCorrQxFHCalFAll   = nullptr;
   TProfile *mhCorrQyFHCalFAll   = nullptr;
   TProfile *mhCorrQxFHCalNAll   = nullptr;
   TProfile *mhCorrQyFHCalNAll   = nullptr;
   TProfile *mhCorrQxFHCalSAll   = nullptr;
   TProfile *mhCorrQyFHCalSAll   = nullptr;
   TProfile *mhCorrQxTPCNAll     = nullptr;
   TProfile *mhCorrQyTPCNAll     = nullptr;
   TProfile *mhCorrQxTPCSAll     = nullptr;
   TProfile *mhCorrQyTPCSAll     = nullptr;

   // For the shift EP correction: <cosNphiEP>, <sinNphiEP>
   TProfile *mhCorrCos1FHCalFAll = nullptr;
   TProfile *mhCorrCos2FHCalFAll = nullptr;
   TProfile *mhCorrCos3FHCalFAll = nullptr;
   TProfile *mhCorrCos4FHCalFAll = nullptr;
   TProfile *mhCorrCos5FHCalFAll = nullptr;
   TProfile *mhCorrCos6FHCalFAll = nullptr;
   TProfile *mhCorrCos7FHCalFAll = nullptr;
   TProfile *mhCorrCos8FHCalFAll = nullptr;
   TProfile *mhCorrSin1FHCalFAll = nullptr;
   TProfile *mhCorrSin2FHCalFAll = nullptr;
   TProfile *mhCorrSin3FHCalFAll = nullptr;
   TProfile *mhCorrSin4FHCalFAll = nullptr;
   TProfile *mhCorrSin5FHCalFAll = nullptr;
   TProfile *mhCorrSin6FHCalFAll = nullptr;
   TProfile *mhCorrSin7FHCalFAll = nullptr;
   TProfile *mhCorrSin8FHCalFAll = nullptr;

   TProfile *mhCorrCos1FHCalNAll = nullptr;
   TProfile *mhCorrCos2FHCalNAll = nullptr;
   TProfile *mhCorrCos3FHCalNAll = nullptr;
   TProfile *mhCorrCos4FHCalNAll = nullptr;
   TProfile *mhCorrCos5FHCalNAll = nullptr;
   TProfile *mhCorrCos6FHCalNAll = nullptr;
   TProfile *mhCorrCos7FHCalNAll = nullptr;
   TProfile *mhCorrCos8FHCalNAll = nullptr;
   TProfile *mhCorrSin1FHCalNAll = nullptr;
   TProfile *mhCorrSin2FHCalNAll = nullptr;
   TProfile *mhCorrSin3FHCalNAll = nullptr;
   TProfile *mhCorrSin4FHCalNAll = nullptr;
   TProfile *mhCorrSin5FHCalNAll = nullptr;
   TProfile *mhCorrSin6FHCalNAll = nullptr;
   TProfile *mhCorrSin7FHCalNAll = nullptr;
   TProfile *mhCorrSin8FHCalNAll = nullptr;

   TProfile *mhCorrCos1FHCalSAll = nullptr;
   TProfile *mhCorrCos2FHCalSAll = nullptr;
   TProfile *mhCorrCos3FHCalSAll = nullptr;
   TProfile *mhCorrCos4FHCalSAll = nullptr;
   TProfile *mhCorrCos5FHCalSAll = nullptr;
   TProfile *mhCorrCos6FHCalSAll = nullptr;
   TProfile *mhCorrCos7FHCalSAll = nullptr;
   TProfile *mhCorrCos8FHCalSAll = nullptr;
   TProfile *mhCorrSin1FHCalSAll = nullptr;
   TProfile *mhCorrSin2FHCalSAll = nullptr;
   TProfile *mhCorrSin3FHCalSAll = nullptr;
   TProfile *mhCorrSin4FHCalSAll = nullptr;
   TProfile *mhCorrSin5FHCalSAll = nullptr;
   TProfile *mhCorrSin6FHCalSAll = nullptr;
   TProfile *mhCorrSin7FHCalSAll = nullptr;
   TProfile *mhCorrSin8FHCalSAll = nullptr;

   TProfile *mhCorrCos1TPCNAll   = nullptr;
   TProfile *mhCorrCos2TPCNAll   = nullptr;
   TProfile *mhCorrCos3TPCNAll   = nullptr;
   TProfile *mhCorrCos4TPCNAll   = nullptr;
   TProfile *mhCorrSin1TPCNAll   = nullptr;
   TProfile *mhCorrSin2TPCNAll   = nullptr;
   TProfile *mhCorrSin3TPCNAll   = nullptr;
   TProfile *mhCorrSin4TPCNAll   = nullptr;

   TProfile *mhCorrCos1TPCSAll   = nullptr;
   TProfile *mhCorrCos2TPCSAll   = nullptr;
   TProfile *mhCorrCos3TPCSAll   = nullptr;
   TProfile *mhCorrCos4TPCSAll   = nullptr;
   TProfile *mhCorrSin1TPCSAll   = nullptr;
   TProfile *mhCorrSin2TPCSAll   = nullptr;
   TProfile *mhCorrSin3TPCSAll   = nullptr;
   TProfile *mhCorrSin4TPCSAll   = nullptr;

   ClassDef(MpdEventPlaneAll, 1);
};
#endif