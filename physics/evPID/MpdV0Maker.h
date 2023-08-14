#ifndef MPDV0MAKER_H
#define MPDV0MAKER_H
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdHelix.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"

class MpdTrack;
class MpdTpcKalmanTrack;
class MpdV0;
class TH2F;

class TVector3;

class MpdV0Maker : public MpdAnalysisTask {

public:
   MpdV0Maker() = default;
   MpdV0Maker(const char *name, const char *outputName = "taskName");
   ~MpdV0Maker() = default;

   void            UserInit();
   void            ProcessEvent(MpdAnalysisEvent &event);
   void            Finish();
   static long int FindCommonParent(long int prim1, long int prim2, TClonesArray *MCTracks);

protected:
   bool selectTrack(MpdTrack *tr);

   bool createSelectV0(MpdTrack *tr1, MpdTpcKalmanTrack *ktr1, MpdTrack *tr2, MpdTpcKalmanTrack *ktr2, MpdV0 &v);

   MpdHelix MakeHelix(const MpdKalmanTrack &tr) const;

   void ArmenterosPodolanski(TVector3 &p1, TVector3 &p2, float &qt, float &alpha) const;

   bool ArmenterosQtCut(float qt, float alpha, MpdParticle &part) const;

   float CosPsiPair(TVector3 &p1, TVector3 &p2) const;

   bool AsymmetryCut(float asym, float pt) const;

   template <typename T>
   T addHist(T h)
   {
      fOutputList->Add(h);
      return h;
   }

private:
   bool                  mMC      = false;
   bool                  mFillEff = true; // Fill QA, efficiency and other histograms
   int                   mDefPDG  = 11;   // PDG code of particles contributing to V0 by default
   TVector3              mPrimaryVertex;  // vertex in current event
   vector<MpdParticle *> mPartK;          // transient list of candidates
   TClonesArray         *mMCTracks = nullptr;
   MpdKalmanHit          mKHit;
   MpdKalmanFilter      *mKF = nullptr;

   int   mNofHitsCut = 10; // Minimal number of TPC hits per track
   float mAlphaCut   = 3.1415;
   float mCosPsiCut  = -1.;
   float mDistCut    = 10.;
   float mEtaCut     = 1.;    // Maximal track rapidity
   float mMassCut    = 2.;    // 2 GeV?
   float mPtminCut   = 0.050; // track pT min
   float mChi2Cut    = 10.;   // Maximal Kalman chi2
   float mMinR2Cut   = 1.;    // Minimal conv radius
   float mMaxR2Cut   = 150.;  // Maximal conv. radius

   TH2F *mhCutEff = nullptr;

   TH2F *mhAlpha  = nullptr;
   TH2F *mhChi2   = nullptr;
   TH2F *mhDist   = nullptr;
   TH2F *mhMassEE = nullptr;
   TH2F *mhCosPsi = nullptr;
   TH2F *mhArmPo  = nullptr; // Armesteros-Podolanski plot
   TH2F *mhAsym   = nullptr; // electron asymmetry
   TH2F *mhConvSp = nullptr; // spectrum of converted photons

   // MCtrue
   TH2F *mhProbElTrue = nullptr; // Electron PID probability estimate for true electrons
   TH2F *mhdEdxTrue   = nullptr; // TPC dEdx for true electrons
   TH2F *mhAlphaTrue  = nullptr;
   TH2F *mhChi2True   = nullptr;
   TH2F *mhDistTrue   = nullptr;
   TH2F *mhMassEETrue = nullptr;
   TH2F *mhCosPsiTrue = nullptr;
   TH2F *mhArmPoTrue  = nullptr; // Armesteros-Podolanski plot
   TH2F *mhAsymTrue   = nullptr; // electron asymmetry
   TH2F *mhConvSpTrue = nullptr; // spectrum of converted photons

   ClassDef(MpdV0Maker, 1);
};
#endif
