#ifndef MPDV0MAKER_H
#define MPDV0MAKER_H
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdHelix.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "TMVA/Reader.h"

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
   void            setBDTCut(float cut = 0.9) { mBDTCut = cut; }
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

   template <typename T>
   T addHistClone(T h, const char *postfix = "")
   {
      auto hClone = (T)h->Clone(Form("%s%s", h->GetName(), postfix));
      fOutputList->Add(hClone);
      return hClone;
   }

private:
   bool                  mMC              = false;
   bool                  mFillEff         = true; // Fill QA, efficiency and other histograms
   bool                  mUseBDTEstimator = true; // Use TMVA BDT classifier to select conversion pairs
   int                   mDefPDG          = 11;   // PDG code of particles contributing to V0 by default
   TVector3              mPrimaryVertex;          // vertex in current event
   vector<MpdParticle *> mPartK;                  // transient list of candidates
   TClonesArray         *mMCTracks = nullptr;
   MpdKalmanHit          mKHit;
   MpdKalmanFilter      *mKF         = nullptr;
   TMVA::Reader         *mTMVAReader = nullptr; // TMVA classificator

   int   mNofHitsCut    = 10; // Minimal number of TPC hits per track
   float mCPACut        = 0.8;
   float mCosPsiCut     = -1.;
   float mDistCut       = 10.;
   float mDedxSigmaCut  = 3.;
   float mTofSigmaCut   = 3.;
   bool  mRequireTOFpid = false;
   float mEtaCut        = 1.;    // Maximal track rapidity
   float mMassCut       = 0.2;   // 2 GeV?
   float mPtminCut      = 0.050; // track pT min
   float mChi2Cut       = 10.;   // Maximal Kalman chi2
   float mMinR2Cut      = 1.;    // Minimal conv radius
   float mMaxR2Cut      = 150.;  // Maximal conv. radius
   float mBDTCut        = -0.26; // BDT maximal significance

   // transient values
   float mPt          = 0.; //! V0 pt
   float mEta         = 0.; //! V0 eta
   float mNcl1        = 0;  //! first track number of TPC clusters, BDT requires float
   float mdEdx1       = 0.; //! first track sigmafied dEdx vrt electron line
   float mNcl2        = 0;  //! second track number of TPC clusters, BDT requires float
   float mdEdx2       = 0.; //! second track sigmafied dEdx vrt electron line
   float mChi2        = 0.; //! Kalman fit chi2
   float mCPA         = 0.; //! Cos of Pointing angle
   float mCPsi        = 0.; //! Cos of pair orientation angle
   float mArmentAlpha = 0.; //! Armenteros alpha
   float mArmentQt    = 0.; //! Armenteros Qt
   float mMee         = 0.; //! pair mass
   float mV0r         = 0.; //! V0 vertex
   float mV0z         = 0.; //! V0 vertex z

   TH2F *mhCutEff = nullptr;

   TH2F *mhCPA    = nullptr;
   TH2F *mhChi2   = nullptr;
   TH2F *mhDist   = nullptr;
   TH2F *mhMassEE = nullptr;
   TH2F *mhCosPsi = nullptr;
   TH2F *mhArmPo  = nullptr; // Armesteros-Podolanski plot
   TH2F *mhAsym   = nullptr; // electron asymmetry
   TH2F *mhConvSp = nullptr; // spectrum of converted photons
   TH2F *mhBDT    = nullptr; // BDT decision

   // MCtrue
   TH2F *mhProbElTrue = nullptr; // Electron PID probability estimate for true electrons
   TH2F *mhdEdxTrue   = nullptr; // TPC dEdx for true electrons
   TH2F *mhCPATrue    = nullptr;
   TH2F *mhChi2True   = nullptr;
   TH2F *mhDistTrue   = nullptr;
   TH2F *mhMassEETrue = nullptr;
   TH2F *mhCosPsiTrue = nullptr;
   TH2F *mhArmPoTrue  = nullptr; // Armesteros-Podolanski plot
   TH2F *mhAsymTrue   = nullptr; // electron asymmetry
   TH2F *mhConvSpTrue = nullptr; // spectrum of converted photons
   TH2F *mhBDTTrue    = nullptr; // BDT decision

   ClassDef(MpdV0Maker, 1);
};
#endif
