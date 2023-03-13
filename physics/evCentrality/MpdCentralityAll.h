#ifndef MPDCENTRALITYALL_H
#define MPDCENTRALITYALL_H

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
#include "MpdCentralityAllParams.h"
#include "TRandom.h"

class MpdTpcKalmanTrack;

class MpdCentralityAll : public MpdAnalysisTask {

public:
   MpdCentralityAll() {}
   MpdCentralityAll(const char *name, const char *outputName = "taskName");
   ~MpdCentralityAll() {} // TODO: normal descructor with cleaning off histos

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

private:
   // Event properties
   bool     isInitialized = false;
   bool     isMC          = true;
   TVector3 mPrimaryVertex;
   TRandom  RND;
   float    TrigEffMult[15] = {0.,       0.745693, 0.806211, 0.860199, 0.912567, 0.951407, 0.96941, 0.983502,
                            0.992123, 0.995581, 0.998052, 0.99874,  0.999466, 0.999515, 0.999899};

   std::string            mParamConfig;
   MpdCentralityAllParams mParams;

   std::string mOutFile = "histos.root";

   MpdKalmanFilter      *mKF              = nullptr;
   TClonesArray         *mMCTracks        = nullptr;
   TObjArray            *mEMCClusters     = nullptr;
   TClonesArray         *mKalmanTracks    = nullptr;
   TClonesArray         *mMpdGlobalTracks = nullptr;
   vector<MpdParticle *> mPartK;
   MpdKalmanHit          mKHit;

   // Histograms
   TList mHistoList;

   // General QA
   TH1F *mhEvents          = nullptr;
   TH1F *mhVertex          = nullptr;
   TH1F *mhVertexAcc       = nullptr;
   TH1F *mhHits            = nullptr;
   TH1F *mhEta             = nullptr;
   TH1F *mhDca             = nullptr;
   TH1F *mhPt              = nullptr;
   TH1F *mhMultiplicity    = nullptr;
   TH1F *mhMultiplicityEff = nullptr;
   TH1F *mhCentrality      = nullptr;
   TH1F *mhCentConvert     = nullptr;
   TH2F *mhTrEff           = nullptr;

   ClassDef(MpdCentralityAll, 1);
};
#endif
