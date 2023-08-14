#ifndef MPDCONVPI0_H
#define MPDCONVPI0_H

#include <deque>

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "MpdEvent.h"
#include "MpdPid.h"
#include "FairMCEventHeader.h"
#include "MpdKalmanFilter.h"
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdPhoton.h"
#include "MpdPhotonAnalysisParams.h"

class TProfile2D;
class TProfile3D;
class MpdV0;

namespace eventCutId {
enum eventCutId { kNone = 0, kVertexPresent, kVertexZ, kGoodCentrality, kNcuts };
}

namespace trackCutId {
enum trackCutId { kNone = 0, kNhits, kEtaPt, kProbElectron, kTpcTofPid, kNcuts };
}

namespace V0CutId {
enum V0CutId {
   kNone = 0,
   kGoodTracks,
   kNonZeroPt,
   kChi2,
   kRconv,
   kAlpha,
   kDist,
   kMass,
   kQt,
   kAsymmetry,
   kCosPsi,
   kPhotonAOD,
   kNcuts
};
}

namespace clustCutId {
enum clustCutId { kNone = 0, kE, kMult, kTof, kDisp, kCPV, kNcuts };
}

class V0 {
public:
   V0() = default;
   V0(bool it, float pt, float chi2, float mee, float rConv, float alpha, float dca, float qt, float qta, float asym,
      float psi)
      : isTrue(it), mpt(pt), mchi2(chi2), mmee(mee), mrConv(rConv), malpha(alpha), mdca(dca), mqt(qt), mqta(qta),
        masym(asym), mpsi(psi)
   {
   }

   ~V0() = default;
   bool  isTrue;
   float mpt;
   float mchi2;
   float mmee;
   float mrConv;
   float malpha;
   float mdca;
   float mqt;
   float mqta;
   float masym;
   float mpsi;
};

class MpdTpcKalmanTrack;

class MpdConvPi0 : public MpdAnalysisTask {

public:
   MpdConvPi0() {}
   MpdConvPi0(const char *name, const char *outputName = "taskName");
   ~MpdConvPi0() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

   void setOutFile(std::string filename = "histos.root") { mOutFile = filename; }

protected:
   bool selectEvent(MpdAnalysisEvent &event);

   void selectConversion(MpdAnalysisEvent &event);

   bool selectTrack(MpdTrack *tr);

   bool selectV0(MpdV0 *v0);

   void selectClusters(MpdAnalysisEvent &event);

   void processHistograms(MpdAnalysisEvent &event);

   bool ArmenterosQtCut(MpdV0 *v0) const;

   bool AsymmetryCut(float asym, float pt) const;

   bool TestHybrid(MpdPhoton &clu, MpdPhoton &v0) const;

   float distCPV(float dphi, float dz, float E) const;

   float lambdaCut(float l1, float l2, float E) const;

   float tofCut(float time, float E) const;

   float Nonlinearity(float oldE) const;

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
   std::vector<V0> mStorage;

   // Event properties
   bool     isMC           = true;
   bool     applySelection = true;
   int      mCenBin        = 0;
   int      mZvtxBin       = 0;
   int      mEPBin         = 0;
   TVector3 mPrimaryVertex;

   std::string             mParamConfig;
   MpdPhotonAnalysisParams mParams;

   std::string mOutFile = "histos.root";

   TClonesArray *mMCTracks        = nullptr;
   TObjArray    *mEMCClusters     = nullptr;
   TClonesArray *mKalmanTracks    = nullptr;
   TClonesArray *mMpdGlobalTracks = nullptr;

   std::vector<MpdPhoton> mClusters; /// EMC clusters in this event
   std::vector<MpdPhoton> mV0;       /// V0 in this event

   static constexpr short nMixEventZ    = 5; /// number of bins in z direction
   static constexpr short nMixEventCent = 5; /// number of bins of centrality
   static constexpr short nMixEventEP   = 5; /// number of bins of Reaction Plane orientation
   static constexpr short nMixTot       = nMixEventZ * nMixEventCent * nMixEventEP;

   int                                        mMaxMixSize = 10;
   std::array<std::deque<MpdPhoton>, nMixTot> mMixClu;
   std::array<std::deque<MpdPhoton>, nMixTot> mMixV0;
   vector<bool>                               isGoodTrack;

   // Histograms
   TList          mHistoList;
   int            mHistoCentBins;
   float          psiEP, psiEPs, psiEPn, psiRP;
   vector<double> mCentBins = {0, 10, 20, 30, 40, 60}; // TODO: make configurable
   TAxis         *mAxCent;
   // General QA
   TH1F *mhEventCutEff = nullptr;
   TH1F *mhVertex      = nullptr;
   TH1F *mhCentrality  = nullptr;
   // TH2F * mhEPvsCen = nullptr ;

   // Single Track selection
   TH2F *mhTrackCutEff               = nullptr;
   TH1F *mhTrackNhits                = nullptr; // Number of hits per track
   TH2F *mhTrackEtaPt                = nullptr; // track occupancy QA
   TH2F *mhTrackProbEl               = nullptr; // Electron PID probability estimate
   TH2F *mhTrackNsigDEdx             = nullptr; // TPC dEdx Nsigma
   TH2F *mhTrackNsigBeta             = nullptr; // TOF beta Nsigma
   TH2F *mhTrackNsigDEdxNsigBeta     = nullptr; // TOF beta Nsigma
   TH1F *mhTrackNhitsTrue            = nullptr; // Number of hits per track
   TH2F *mhTrackEtaPtTrue            = nullptr; // track occupancy QA
   TH2F *mhTrackProbElTrue           = nullptr; // Electron PID probability estimate for true electrons
   TH2F *mhTrackNsigDEdxTrue         = nullptr; // TPC dEdx Nsigma for true electrons
   TH2F *mhTrackNsigBetaTrue         = nullptr; // TOF beta Nsigma for true electrons
   TH2F *mhTrackNsigDEdxNsigBetaTrue = nullptr; // TOF beta Nsigma

   // V0selection
   TH2F *mhV0CutEff = nullptr;
   TH3F *mhConvMap  = nullptr; // Conversion map for all V0
   TH2F *mhAlpha    = nullptr;
   TH2F *mhV0rConv  = nullptr;
   TH2F *mhChi2     = nullptr;
   TH2F *mhDist     = nullptr;
   TH2F *mhMassEE   = nullptr;
   TH2F *mhCosPsi   = nullptr;
   TH2F *mhArmPo    = nullptr; // Armesteros-Podolanski plot
   TH2F *mhAsym     = nullptr; // electron asymmetry
   TH2F *mhConvSp   = nullptr; // spectrum of converted photons

   TH3F *mhConvMapTrue = nullptr; // Conversion map for true conversions
   TH2F *mhAlphaTrue   = nullptr;
   TH2F *mhChi2True    = nullptr;
   TH2F *mhV0rConvTrue = nullptr;
   TH2F *mhDistTrue    = nullptr;
   TH2F *mhMassEETrue  = nullptr;
   TH2F *mhCosPsiTrue  = nullptr;
   TH2F *mhArmPoTrue   = nullptr; // Armesteros-Podolanski plot
   TH2F *mhAsymTrue    = nullptr; // electron asymmetry
   TH2F *mhConvSpTrue  = nullptr; // spectrum of converted photons

   // Cluster selection
   TH2F *mhCluCutEff   = nullptr;
   TH2F *mhCluMult     = nullptr;
   TH2F *mhCluTofSigma = nullptr;
   TH2F *mhCluDistCPV  = nullptr;

   // Inv mass histos
   std::vector<TH2F *> mhRealCalo;
   std::vector<TH2F *> mhRealHybrid;
   std::vector<TH2F *> mhRealConv;

   std::vector<TH2F *> mhMixedCalo;
   std::vector<TH2F *> mhMixedHybrid;
   std::vector<TH2F *> mhMixedConv;

   std::vector<TH2F *> mhRealCaloTruePi;
   std::vector<TH2F *> mhRealHybridTruePi;
   std::vector<TH2F *> mhRealConvTruePi;

   std::vector<TH2F *> mhRealCaloTruePh;
   std::vector<TH2F *> mhRealHybridTruePh;
   std::vector<TH2F *> mhRealConvTruePh;

   std::vector<TH2F *> mhRealCaloTrueAll;
   std::vector<TH2F *> mhRealHybridTrueAll;
   std::vector<TH2F *> mhRealConvTrueAll;

   std::vector<TH2F *> hPrimPi;
   std::vector<TH2F *> hPrimPh;

   // flow
   std::vector<TProfile3D *> mp3V1etaPtMinvCalo;
   std::vector<TProfile3D *> mp3V1etaPtMinvCaloTruePi;
   std::vector<TProfile3D *> mp3V1etaPtMinvCaloTruePh;
   std::vector<TProfile3D *> mp3V1etaPtMinvHybrid;
   std::vector<TProfile3D *> mp3V1etaPtMinvHybridTruePi;
   std::vector<TProfile3D *> mp3V1etaPtMinvHybridTruePh;
   std::vector<TProfile3D *> mp3V1etaPtMinvConv;
   std::vector<TProfile3D *> mp3V1etaPtMinvConvTruePi;
   std::vector<TProfile3D *> mp3V1etaPtMinvConvTruePh;
   std::vector<TProfile2D *> mp2V1etaPtPrimPiEP;
   std::vector<TProfile2D *> mp2V1etaPtPrimPhEP;
   std::vector<TProfile2D *> mp2V1etaPtPrimPiRP;
   std::vector<TProfile2D *> mp2V1etaPtPrimPhRP;

   ClassDef(MpdConvPi0, 1);
};
#endif
