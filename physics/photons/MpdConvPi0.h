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
   kCPA,
   kDist,
   kMass,
   kQt,
   kAsymmetry,
   kCosPsi,
   kPhotonAOD,
   kBDT,
   kNcuts
};
}

namespace clustCutId {
enum clustCutId { kNone = 0, kE, kMult, kTof, kDisp, kCPV, kNcuts };
}

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

   void fillMC(MpdAnalysisEvent &event);

   void processHistograms(MpdAnalysisEvent &event);

   bool ArmenterosQtCut(MpdV0 *v0) const;

   bool AsymmetryCut(float asym, float pt) const;

   bool TestHybrid(MpdPhoton &clu, MpdPhoton &v0) const;

   float distCPV(float dphi, float dz, float E) const;

   float lambdaCut(float l1, float l2, float E) const;

   float tofCut(float time, float E) const;

   float smearTime(float time, float E) const;

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
   // Event properties
   bool     isMC           = true;
   bool     applySelection = true;
   int      mCenBin        = 0;
   float    mCentrality    = 0.;
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

   static constexpr short nMixEventZ    = 10; /// number of bins in z direction
   static constexpr short nMixEventCent = 6;  /// number of bins of centrality
   static constexpr short nMixEventEP   = 5;  /// number of bins of Reaction Plane orientation
   static constexpr short nMixBins      = nMixEventZ * nMixEventCent * nMixEventEP;

   int                                                      mMaxMixSize = 2;
   std::array<std::deque<std::vector<MpdPhoton>>, nMixBins> mMixEventsV0;
   std::array<std::deque<std::vector<MpdPhoton>>, nMixBins> mMixEventsClusters;

   vector<bool> isGoodTrack;

   // Histograms
   TList          mHistoList;
   int            mHistoCentBins;
   float          psiEP, psiEPs, psiEPn, psiRP;
   vector<double> mCentBins = {0, 10, 20, 30, 40, 60, 80}; // TODO: make configurable
   TAxis         *mAxCent;
   // General QA
   TH1F *hEventCutEff = nullptr;
   TH1F *hVertex      = nullptr;
   TH1F *hCentrality  = nullptr;
   // TH2F * hEPvsCen = nullptr ;

   // Single Track selection
   TH3F *h3trackEtaPtCutEff = nullptr;

   TH3F *h3trackEtaPtNhits       = nullptr;     // Number of hits per track
   TH2F *h2trackEtaPt            = nullptr;     // track occupancy QA
   TH3F *h3trackEtaPtProbEl      = nullptr;     // Electron PID probability estimate
   TH3F *h3trackEtaPtNsigDEdx    = nullptr;     // TPC dEdx Nsigma
   TH3F *h3trackEtaPtNsigBeta    = nullptr;     // TOF beta Nsigma
   TH2F *h2trackNsigDEdxNsigBeta = nullptr;     // TOF beta Nsigma
                                                //
   TH3F *h3trackEtaPtNhitsTrue       = nullptr; // Number of hits per track
   TH2F *h2trackEtaPtTrue            = nullptr; // track occupancy QA
   TH3F *h3trackEtaPtProbElTrue      = nullptr; // Electron PID probability estimate for true electrons
   TH3F *h3trackEtaPtNsigDEdxTrue    = nullptr; // TPC dEdx Nsigma for true electrons
   TH3F *h3trackEtaPtNsigBetaTrue    = nullptr; // TOF beta Nsigma for true electrons
   TH2F *h2trackNsigDEdxNsigBetaTrue = nullptr; // TOF beta Nsigma

   // V0selection
   TH3F *h3v0EtaPtCutEff = nullptr;

   TH3F *h3v0ConvMap     = nullptr; // Conversion map for all V0
   TH3F *h3v0EtaPtCPA    = nullptr;
   TH3F *h3v0EtaPtRconv  = nullptr;
   TH3F *h3v0EtaPtChi2   = nullptr;
   TH3F *h3v0EtaPtDist   = nullptr;
   TH3F *h3v0EtaPtMassEE = nullptr;
   TH3F *h3v0EtaPtCosPsi = nullptr;
   TH2F *h2v0ArmPo       = nullptr; // Armesteros-Podolanski plot
   TH3F *h3v0EtaPtAsym   = nullptr; // electron asymmetry
   TH2F *h2v0EtaPt       = nullptr; // spectrum of converted photons

   TH3F *h3v0ConvMapTrue     = nullptr; // Conversion map for true conversions
   TH3F *h3v0EtaPtCPATrue    = nullptr;
   TH3F *h3v0EtaPtChi2True   = nullptr;
   TH3F *h3v0EtaPtRconvTrue  = nullptr;
   TH3F *h3v0EtaPtDistTrue   = nullptr;
   TH3F *h3v0EtaPtMassEETrue = nullptr;
   TH3F *h3v0EtaPtCosPsiTrue = nullptr;
   TH2F *h2v0ArmPoTrue       = nullptr; // Armesteros-Podolanski plot
   TH3F *h3v0EtaPtAsymTrue   = nullptr; // electron asymmetry
   TH2F *h2v0EtaPtTrue       = nullptr; // spectrum of converted photons

   // Cluster selection
   TH3F *h3cluEtaECutEff   = nullptr;
   TH3F *h3cluEtaEMult     = nullptr;
   TH3F *h3cluEtaETofSigma = nullptr;
   TH3F *h3cluEtaEDistCPV  = nullptr;

   // Inv mass histos
   std::vector<TH3F *> h3yPtMinvRealCalo;
   std::vector<TH3F *> h3yPtMinvRealHybrid;
   std::vector<TH3F *> h3yPtMinvRealConv;

   std::vector<TH3F *> h3yPtMinvMixedCalo;
   std::vector<TH3F *> h3yPtMinvMixedHybrid;
   std::vector<TH3F *> h3yPtMinvMixedConv;
   TH1F               *hMixBinOccupancy;
   TH3F               *h3MixBinOccupancy;
   TH3F               *h3BinOccupancyRealCalo;
   TH3F               *h3BinOccupancyMixedCalo;
   TH3F               *h3BinOccupancyRealHybrid;
   TH3F               *h3BinOccupancyMixedHybrid;
   TH3F               *h3BinOccupancyRealConv;
   TH3F               *h3BinOccupancyMixedConv;

   std::vector<TH3F *> h3yPtMinvRealCaloTruePi;
   std::vector<TH3F *> h3yPtMinvRealHybridTruePi;
   std::vector<TH3F *> h3yPtMinvRealConvTruePi;

   std::vector<TH3F *> h3yPtMinvRealCaloTruePh;
   std::vector<TH3F *> h3yPtMinvRealHybridTruePh;
   std::vector<TH3F *> h3yPtMinvRealConvTruePh;

   std::vector<TH3F *> h3yPtMinvRealCaloTrueV0;
   std::vector<TH3F *> h3yPtMinvRealHybridTrueV0;
   std::vector<TH3F *> h3yPtMinvRealConvTrueV0;

   std::vector<TH2F *> h2yPtPrimPi;
   std::vector<TH2F *> h2yPtPrimPh;

   // flow
   TProfile                 *pR1cent;
   TProfile                 *pR1centTrue;
   std::vector<TProfile3D *> p3v1yPtMinvCalo;
   std::vector<TProfile3D *> p3v1yPtMinvCaloTruePi;
   std::vector<TProfile3D *> p3v1yPtMinvCaloTruePh;
   std::vector<TProfile3D *> p3v1yPtMinvHybrid;
   std::vector<TProfile3D *> p3v1yPtMinvHybridTruePi;
   std::vector<TProfile3D *> p3v1yPtMinvHybridTruePh;
   std::vector<TProfile3D *> p3v1yPtMinvConv;
   std::vector<TProfile3D *> p3v1yPtMinvConvTruePi;
   std::vector<TProfile3D *> p3v1yPtMinvConvTruePh;
   std::vector<TProfile2D *> p2v1yPtPrimPi;
   std::vector<TProfile2D *> p2v1yPtPrimPh;

   // flow RP
   std::vector<TProfile3D *> p3v1yPtMinvCaloRP;
   std::vector<TProfile3D *> p3v1yPtMinvCaloTruePiRP;
   std::vector<TProfile3D *> p3v1yPtMinvCaloTruePhRP;
   std::vector<TProfile3D *> p3v1yPtMinvHybridRP;
   std::vector<TProfile3D *> p3v1yPtMinvHybridTruePiRP;
   std::vector<TProfile3D *> p3v1yPtMinvHybridTruePhRP;
   std::vector<TProfile3D *> p3v1yPtMinvConvRP;
   std::vector<TProfile3D *> p3v1yPtMinvConvTruePiRP;
   std::vector<TProfile3D *> p3v1yPtMinvConvTruePhRP;
   std::vector<TProfile2D *> p2v1yPtPrimPiRP;
   std::vector<TProfile2D *> p2v1yPtPrimPhRP;

   ClassDef(MpdConvPi0, 1);
};
#endif
