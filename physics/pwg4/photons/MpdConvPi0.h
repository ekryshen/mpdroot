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
#include "MpdHelix.h"
#include "FairMCEventHeader.h"
#include "MpdKalmanFilter.h"
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdPhoton.h"
#include "MpdPhotonAnalysisParams.h"

  class V0{
  public:
    V0()=default ;
    V0(bool it,float pt,float chi2,float mee,float rConv,float alpha,float dca, float qt, float qta,float asym,float psi):
    isTrue(it),mpt(pt),mchi2(chi2),mmee(mee),mrConv(rConv),malpha(alpha),mdca(dca),mqt(qt),mqta(qta),masym(asym),mpsi(psi){}

    ~V0()=default ;
    bool isTrue;
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


class MpdConvPi0 : public MpdAnalysisTask {

  public : 
  MpdConvPi0() {}
  MpdConvPi0(const char *name, const char *outputName="taskName") ;
  ~MpdConvPi0() {}  //TODO: normal descructor with cleaning off histos

  void UserInit();
  void ProcessEvent(MpdAnalysisEvent &event) ;
  void Finish() ;

  void setOutFile(std::string filename = "histos.root"){mOutFile=filename ;}


 protected:

  bool selectEvent(MpdAnalysisEvent &event) ; 	

  void selectConversion(MpdAnalysisEvent &event) ;

  bool selectTrack(MpdTrack *tr) ;

  bool createSelectV0(MpdTrack *tr1, MpdTpcKalmanTrack *ktr1,  MpdTrack *tr2, MpdTpcKalmanTrack *ktr2, MpdPhoton &v );

  void selectClusters(MpdAnalysisEvent &event); 

  void processHistograms(MpdAnalysisEvent &event) ;

  int centrality(int tpcMultiplicity) const ;
 
  MpdHelix MakeHelix(const MpdKalmanTrack &tr) const ;

  float Beta_sigma(float beta, float mom) const ;
 
  float dEdx_sigma(float dEdx, float mom) const ;

  void ArmenterosPodolanski(TVector3& p1, TVector3& p2, float &qt, float &alpha ) const ;

  bool ArmenterosQtCut(float qt, float alpha, MpdParticle & part) const ;

  float CosPsiPair(TVector3& p1, TVector3& p2) const ;

  bool AsymmetryCut(float asym, float pt) const ;
  
  long int IsSameParent(long int prim1, long int prim2) const ;
  
  bool TestHybrid(MpdPhoton &clu, MpdPhoton &v0) const ;

  float distCPV(float dphi, float dz, float E) const ;

  float lambdaCut(float l1,float l2,float E) const; 

  float tofCut(float time,float E) const;

  float Nonlinearity(float oldE) const;
private: 

   std::vector<V0> mStorage ;

   //Event properties
   bool isInitialized = false ;
   bool isMC  =  true ;
   int mCenBin =0 ;
   int mZvtxBin=0 ;
   int mRPBin =0 ;
   TVector3 mPrimaryVertex;

   std::string mParamConfig ;
   MpdPhotonAnalysisParams mParams ;

   std::string mOutFile = "histos.root" ;

   MpdKalmanFilter* mKF          = nullptr ;
   MpdPid *mPID                  = nullptr ;
   TClonesArray * mMCTracks      = nullptr ;
   TObjArray * mEMCClusters      = nullptr ;
   TClonesArray * mKalmanTracks  = nullptr ;
   TClonesArray * mMpdGlobalTracks=nullptr ;
   vector<MpdParticle*> mPartK ;
   MpdKalmanHit   mKHit  ;

   std::vector<MpdPhoton> mClusters ; /// EMC clusters in this event
   std::vector<MpdPhoton> mV0 ;       /// EMC clusters in this event

   static constexpr short nMixEventZ = 20 ;      /// number of bins in z direction
   static constexpr short nMixEventCent = 10 ;   /// number of bins of centrality
   static constexpr short nMixEventRP = 5 ;      /// number of bins of Reaction Plane orientation
   static constexpr short nMixTot = nMixEventZ* nMixEventCent * nMixEventRP ;	

   int mMaxMixSize = 10 ; 
   std::array<std::deque<MpdPhoton>, nMixTot> mMixClu ;
   std::array<std::deque<MpdPhoton>, nMixTot> mMixV0 ;

   //Histograms
   TList mHistoList; 
   static constexpr short mHistoCentBins=5;
   //General QA
   TH1F * mhEvents = nullptr ;
   TH1F * mhVertex = nullptr ;
   TH1F * mhCentrality=nullptr ;
   //TH2F * mhEPvsCen = nullptr ;  

   //Single Track selection
   TH1F * mhNhits   = nullptr ;    // Number of hits per track
   TH2F * mhTracks  = nullptr ;    // track occupancy QA
   TH2F * mhProbEl  = nullptr ;    // Electron PID probability estimate
   TH2F * mhdEdx    = nullptr ;    // TPC dEdx
   //V0selection
   TH2F * mhCutEff  = nullptr ;
   TH3F * mhConvMap = nullptr ;    //Conversion map for all V0
   TH2F * mhAlpha = nullptr ;
   TH2F * mhChi2  = nullptr ;
   TH2F * mhDist  = nullptr ;
   TH2F * hmMassEE= nullptr ;
   TH2F * mhCosPsi= nullptr ;
   TH2F * mhArmPo = nullptr ;    //Armesteros-Podolanski plot
   TH2F * mhAsym  = nullptr ;    // electron asymmetry
   TH2F * mhConvSp= nullptr ;    // spectrum of converted photons


   TH2F * mhProbElTrue  = nullptr ;   // Electron PID probability estimate for true electrons
   TH2F * mhdEdxTrue    = nullptr ;   // TPC dEdx for true electrons
   TH3F * mhConvMapTrue = nullptr ;   //Conversion map for true conversions
   TH2F * mhAlphaTrue   = nullptr ;
   TH2F * mhChi2True    = nullptr ;
   TH2F * mhDistTrue    = nullptr ;
   TH2F * hmMassEETrue  = nullptr ;
   TH2F * mhCosPsiTrue  = nullptr ;
   TH2F * mhArmPoTrue   = nullptr ;    //Armesteros-Podolanski plot
   TH2F * mhAsymTrue    = nullptr ;    // electron asymmetry
   TH2F * mhConvSpTrue  = nullptr ;    // spectrum of converted photons

   //Cluster selection
   TH2F * mhCluCutEff  = nullptr ;

   // Inv mass histos
   TH2F * mhRealCalo[mHistoCentBins] ;
   TH2F * mhMixedCalo[mHistoCentBins] ;
   TH2F * mhRealHybrid[mHistoCentBins] ;
   TH2F * mhMixedHybrid[mHistoCentBins] ;
   TH2F * mhRealConv[mHistoCentBins] ;
   TH2F * mhMixedConv[mHistoCentBins] ;

   TH2F * mhRealCaloTrue[mHistoCentBins] ;
   TH2F * mhRealHybridTrue[mHistoCentBins] ;
   TH2F * mhRealConvTrue[mHistoCentBins] ;

   TH2F * mhRealCaloTrueEl[mHistoCentBins] ;
   TH2F * mhRealHybridTrueEl[mHistoCentBins] ;
   TH2F * mhRealConvTrueEl[mHistoCentBins] ;

   TH2F * mhRealCaloTrueAll[mHistoCentBins] ;
   TH2F * mhRealHybridTrueAll[mHistoCentBins] ;
   TH2F * mhRealConvTrueAll[mHistoCentBins] ;

   TH2F * hPrimPi0[mHistoCentBins] ;

  ClassDef(MpdConvPi0, 1) ;

} ;
#endif
