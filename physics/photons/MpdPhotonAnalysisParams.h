#ifndef MPDPHOTONANALYSISPARAMS_H
#define MPDPHOTONANALYSISPARAMS_H

#include <map>
#include <string>

#include "Rtypes.h"
#include "TObject.h"

class MpdPhotonAnalysisParams : public TObject {

public:
   bool mApplySelection = true;
   // Event selection cuts
   float mZvtxCut  = 40.; /// event selection cut (cm)
   int   mNhitsCut = 10;  /// number of hits in TPC tracks used for centrality

   // V0 cuts
   float       mMinR2Cut      = 5.;   // (cm) Minimal conversion radius (to exclude Dalitz)
   float       mMaxR2Cut      = 120.; // (cm) Maximal conversion radius (to exclude poorly reconstructed tracks)
   float       mPIDsigM       = 4.0;
   float       mPIDsigE       = 4.0;
   float       mPIDenergy     = 11.;
   float       mPIDkoeff      = 1.;
   std::string mPIDgenerator  = "NSIG";
   std::string mPIDtracking   = "CF";
   std::string mPIDparticles  = "elpikapr";
   int         mNofHitsCut    = 10;   // minimal number of hits to accept track
   float       mEtaCut        = 1.0;  // maximal pseudorapidity accepted
   float       mPtminCut      = 0.05; // minimal pt used in analysis
   float       mProbElCut     = -1;   // minimal dEdx probability for electrons (does not seem to be filled)
   float       mdEdxSigmaCut  = 3.0;  // dEdx cut in sigmas
   float       mBetaSigmaCut  = 3.0;  // beta cut
   bool        mRequireTOFpid = true; // RequireTOFpid
   float       mMassCut       = 0.051;
   float       mDistCut       = 2.8;        // maximal closest distance between daughters
   float       mCosPsiCut     = 0.96242520; // e+e- pair orientation wrt B-filed
   float       mCPACut        = 0.95;       // Cos Pointing Angle
   float       mChi2Cut       = 10;
   bool        mUseBDT        = true; // Use Boosted Decision Tree
   bool        mUseBDTRegP    = true; // Use BDT estimate (regression) of P
   float       mBDTCut        = 0.17; // BDT classificator cut [-1.,1.], recomended values [-0.2,0.2]

   // Cluster cuts
   float mCluEmin     = 0.05;
   int   mCluMult     = 2;
   float mCluTofMin   = -4.;  // cluster time cut in sigmas
   float mCluTofMax   = 0.;   // cluster time cut in sigmas
   float mCluDisp     = 6.25; // cluster disp cut in sigma squared
   float mCluDispEmin = 0.5;  // Min energy to apply Disp cut (GeV)
   float mCluCPV      = 6.25; // neutrality cut in sigma squared

   void ReadFromFile(std::string fname = "pi0.par");
   void Print() const;

protected:
   void read(std::string name, bool &b);
   void read(std::string name, int &b);
   void read(std::string name, float &b);
   void read(std::string name, std::string &b);

   std::map<std::string, std::string> mMap;

   ClassDef(MpdPhotonAnalysisParams, 1);
};
#endif // MPDPHOTONANALYSISPARAMS_H
