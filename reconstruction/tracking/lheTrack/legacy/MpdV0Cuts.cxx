#include "MpdV0Cuts.h"
using namespace MpdV0Cuts;

/***Mother cuts***/                    // less or more for cuts
Double_t MpdV0Cuts::V0Chi2Cut = 1000.; // 10     //*less than chi2 are selected
Double_t MpdV0Cuts::cosaCut =
   0.996; //-0.8   //*more pointing angle (P,R)               MC vertices are near -1 but many true are at 1
Double_t MpdV0Cuts::RminCut    = 0.;   // 6   //distance between primary and secondary
Double_t MpdV0Cuts::RmaxCut    = 1000; // 50
Double_t MpdV0Cuts::minmassCut = 1.08; // 1.08
Double_t MpdV0Cuts::maxmassCut = 1.15; // 1.15
Double_t MpdV0Cuts::DCAV0Cut   = 50;   // 10
Double_t MpdV0Cuts::PminCut    = 0.;   // 0.3

/***Daughter cuts***/
Double_t MpdV0Cuts::NofHitsCut        = 10;    // 10
Double_t MpdV0Cuts::protonChi2vminCut = 0;     // 40
Double_t MpdV0Cuts::protonChi2vmaxCut = 1000;  // 500 proton pion or pos neg?
Double_t MpdV0Cuts::pionChi2vminCut   = 0.;    // 30
Double_t MpdV0Cuts::pionChi2vmaxCut   = 1000.; // 500
Double_t MpdV0Cuts::IPposCut = 0.; // 0.3*more dca of track to primary vertex.  before v0fit or after???     Should
                                   // probably be used IPpid instead of IPcharge
Double_t MpdV0Cuts::IPnegCut       = 0.; // 0*more
Double_t MpdV0Cuts::protonIPCut    = 0.;
Double_t MpdV0Cuts::pionIPCut      = 0.;
Double_t MpdV0Cuts::ProbCut        = 0.; // 0.25
Double_t MpdV0Cuts::protonPtminCut = 0.; // 0.1*less momentum cuts     min cut
Double_t MpdV0Cuts::protonPminCut  = 0.; // 0.1    //min cut
Double_t MpdV0Cuts::pionPtminCut   = 0.; // 0.03
Double_t MpdV0Cuts::pionPminCut    = 0.; // 0.03
Double_t MpdV0Cuts::protonPtmaxCut = 10.;
Double_t MpdV0Cuts::protonPmaxCut  = 10.;
Double_t MpdV0Cuts::pionPtmaxCut   = 10.;
Double_t MpdV0Cuts::pionPmaxCut    = 10.;

Int_t MpdV0Cuts::protonPDG = 2212; // proton+
Int_t MpdV0Cuts::pionPDG   = 211;  // pi+
Int_t MpdV0Cuts::lambdaPDG = 3122;

/*
 /***Mother cuts***/    // less or more for cuts
/*
 Double_t MpdV0Cuts::V0Chi2Cut           = 11.; //10     //*less than chi2 are selected
  Double_t MpdV0Cuts::cosaCut             = 0.996;    //-0.8   //*more pointing angle (P,R)               MC vertices
 are near -1 but many true are at 1 Double_t MpdV0Cuts::RminCut             = 2.6;    //6   //distance between primary
 and secondary Double_t MpdV0Cuts::RmaxCut             = 1000;     //50 Double_t MpdV0Cuts::minmassCut          = 1.08;
 //1.08 Double_t MpdV0Cuts::maxmassCut          = 1.15; //1.15 Double_t MpdV0Cuts::DCAV0Cut            = 3;     //10
  Double_t MpdV0Cuts::PminCut             = 0.; //0.3

*/

/***Daughter cuts***/
/*
Double_t MpdV0Cuts::NofHitsCut          = 10; //10
Double_t MpdV0Cuts::protonChi2vminCut   = 7.6;        //40
Double_t MpdV0Cuts::protonChi2vmaxCut   = 1000; //500 proton pion or pos neg?
Double_t MpdV0Cuts::pionChi2vminCut     = 6.;    //30
Double_t MpdV0Cuts::pionChi2vmaxCut     = 1000.;  //500



*/