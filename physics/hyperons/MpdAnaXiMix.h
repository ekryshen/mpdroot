#ifndef MPDANAXIMIX_H
#define MPDANAXIMIX_H

/// \ingroup physics
/// \class MpdAnaXiMix
/// \brief Lambda and Xi reconstruction with event mixing
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - 21-03-2023

#include "MpdLambda.h"
#include "MpdXi.h"
#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdHelix.h"
#include "FairMCEventHeader.h"
#include "MpdKalmanFilter.h"
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdTpcKalmanTrack.h"

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

//#include <deque>
#include <map>

class MpdTrackFinderIts5spd;
class MpdTpcKalmanFilter;
class TMinuit;

//#define ITS
#ifdef ITS
typedef MpdItsKalmanTrack AzTrack;
#else
typedef MpdTpcKalmanTrack AzTrack;
#endif

class MpdAnaXiMix : public MpdAnalysisTask {

public:
   struct Params {
   public:
      Params() {}
      const Int_t    pdgCodePr  = 2212; // proton
      const Int_t    pdgCodeNeg = -211; // pi-
      const Int_t    pdgCodePos = 211;  // pi+
      const Int_t    pdgCodeL0  = 3122; // lambda (1.11568 GeV)
      const Int_t    pdgCodeXi  = 3312; // Xi- (1.3213 GeV)
      const Int_t    pdgCodeOm  = 3334; // Omega- (1.67245 GeV)
      const Int_t    pdgCodeKm  = -321; // K- (0.4937 GeV)
      const Int_t    pdgCodeK0  = 310;  // K0 (0.4976 GeV)
      const int      nMix       = 5;    // number of events to mix
      const Double_t gC2p       = 3.;   // 4.; //9.;           //chi2 of p to PV
      const Double_t gC2pi      = 5.;   // 5.; //11.;          //chi2 of pion to PV
      const Double_t gPathL     = 0.0;  // cm - path to Lambda decay
      const Double_t gPathXi    = 0.;   // cm - path to Xi decay
      const Double_t gC2L       = 25.;  // 9999.;  //chi2 between pion & p in V0
      const Double_t gC2Xi      = 25.;  // 16.;//9999.;  //chi2 between pion & L in Xi
      const Double_t gDecayOm   = 0.;   // 1.0;
   };

   MpdAnaXiMix();
   MpdAnaXiMix(const char *name, const char *outputName = "taskName");
   ~MpdAnaXiMix() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

   void SetOutFile(std::string filename = "histos.root"){}; // mOutFile = filename; }

protected:
   bool        SelectEvent(MpdAnalysisEvent &event);
   void        SelectTracks(MpdAnalysisEvent &event);
   void        RecoEff(vector<Int_t> &vecP, vector<Int_t> &vecPi, Int_t pid = 0);
   void        ApplyPid(vector<Int_t> &vecP, vector<Int_t> &vecPi);
   void        BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle *> &vecL);
   MpdHelix    MakeHelix(const MpdKalmanTrack *tr);
   MpdHelix    MakeHelix(const MpdParticle *part);
   void        BuildCascade(vector<Int_t> &vecK, vector<Int_t> &vecPi, vector<MpdParticle *> &vecL);
   Double_t    DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu);
   static void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
   /*
   bool selectEvent(MpdAnalysisEvent &event);

   void selectPosTrack(MpdAnalysisEvent &event);

   bool selectTrack(MpdTrack *tr);

   void selectNegTrack(MpdAnalysisEvent &event);

   void processHistograms(MpdAnalysisEvent &event);

   int centrality(int tpcMultiplicity) const;

   long int IsSameParent(long int prim1, long int prim2) const;

   float dEdx_sigma_K(float dEdx, float mom) const;

   float Beta_sigma_K(float beta, float mom) const;
   */

private:
   bool                   fisInitialized = false;
   TClonesArray          *fMcTracks, *fItsTracks, *fTofMatches, *fMpdTracks;
   MpdTrackFinderIts5spd *fRecoIts;
   MpdTpcKalmanFilter    *fRecoTpc;
   MpdPid                *fPid;

   std::vector<MpdLambda>                        fvLambdas;
   std::vector<MpdLambda>                       *fvvvL;
   std::vector<MpdXi>                            fvXis;
   std::vector<MpdXi>                           *fvvvXi;
   std::vector<tuple<int, float, float, float>>  fvLambMpdgPtEtaY;
   std::vector<tuple<int, float, float, float>> *fvvvLpt;
   std::vector<tuple<int, float, float, float>>  fvXiMpdgPtEtaY;
   std::vector<tuple<int, float, float, float>> *fvvvXipt;
   map<int, MpdVertex>                           fMapVertexEvent; // for event mixing
   multimap<Int_t, AzTrack>                      fMapPiEvent;     // for event mixing
   vector<pair<Double_t, Double_t>>              fVecL1, fVecL2;
   map<int, int>                                 fLays;

   TTree     *fTree;
   Float_t    fB0, fCentr, fZv, fZvgen;
   int        fNtr13, fNv;
   int        fEvNo;
   MpdVertex *fMpdVert;
   Params     fParams;
   TMinuit   *fMinuit;
   // static TVector3 fVtxN, fMomN;
   // static MpdHelix *fTrC;

   std::unordered_map<string, TH1F *> fHistosF;

   ClassDef(MpdAnaXiMix, 1);
};
#endif
