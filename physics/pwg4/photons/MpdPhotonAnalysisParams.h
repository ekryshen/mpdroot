#ifndef MPDPHOTONANALYSISPARAMS_H
#define MPDPHOTONANALYSISPARAMS_H

#include <map> 
#include <string>

#include "Rtypes.h"
#include "TObject.h"

class MpdPhotonAnalysisParams: public TObject{

  public:
  	//
   //Event selection cuts
   float mZvtxCut = 40.;        /// event selection cut (cm) 
   int mNhitsCut = 10;          /// number of hits in TPC tracks used for centrality

   //V0 cuts
   float  mMinR2Cut     = 5. ;   // (cm) Minimal conversion radius (to exclude Dalitz)
   float  mMaxR2Cut     = 120. ; // (cm) Maximal conversion radius (to exclude poorly reconstructed tracks)
   float  mPIDsigM 	 	= 4.0;
   float  mPIDsigE 	 	= 4.0;
   float  mPIDenergy 	= 11.;
   float  mPIDkoeff	 	= 1.;
   std::string mPIDgenerator = "NSIG";
   std::string mPIDtracking  = "CF";
   std::string mPIDparticles = "elpikapr";
   int    mNofHitsCut   = 10 ;  // minimal number of hits to accept track
   float  mEtaCut       = 1.0 ; // maximal pseudorapidity accepted
   float  mPtminCut     = 0.05 ;// minimal pt used in analysis
   float  mProbElCut    = 0.75; 
   float  mdEdxSigmaCut = 3.0 ; //dEdx cut in sigmas
   float  mBetaSigmaCut = 3.0 ; //beta cut
   float  mMassCut      = 0.051; 
   float  mDistCut      = 2.8;   // maximal closest distance between daughters
   float  mCosPsiCut    = 0.96242520 ;  // e+e- pair orientation wrt B-filed
   float  mAlphaCut     = 0.102 ; //r vs p angle
   float mChi2Cut=10 ;

   //Cluster cuts
   float mCluEmin = 0.05 ;
   int   mCluMult = 2 ;
   float mCluTof  = 5.;     //cluster time cut in sigmas 
   float mCluDisp = 6.25;   //cluster disp cut in sigma squared
   float mCluDispEmin=0.5;  //Min energy to apply Disp cut (GeV)
   float mCluCPV  = 6.25;   //neutrality cut in sigma squared

   void ReadFromFile(std::string fname="ConvDef.txt") ;
   void Print() const ;

  protected:
   void read(std::string name, bool &b) ;
   void read(std::string name, int  &b) ;
   void read(std::string name, float &b) ; 
   void read(std::string name, std::string &b);

   std::map<std::string,std::string> mMap ;
  
  ClassDef(MpdPhotonAnalysisParams, 1) ;

};
#endif //MPDPHOTONANALYSISPARAMS_H
