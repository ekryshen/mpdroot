// =======================================================================================================================================
// Supplementary class for Forward-backward correlations analysis
// Author: Igor Altsybeev (SPbSU), 10.2020
// =======================================================================================================================================

#include "TH1.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TString.h"

#include "WinPairBase.h"

#include <iostream>
#include <map>
using namespace std;


const char *varNames[] = {
    "Nevents",

    "f_Nevents",
    "b_Nevents",

//    "x_Nevents",
//    "y_Nevents",

    "fb_Nevents", // WAS COMMENTED BY MISTAKE! uncommented 23.04.2019
//    "xy_Nevents",

//    "fx_Nevents",
//    "fy_Nevents",
//    "bx_Nevents",
//    "by_Nevents",

    // f, b
    "Nf",
    "Nb",
    "N2_f",
    "N2_b",
    "Nf_Nb",

    // x, y
//    "Nx",
//    "Ny",
//    "Nx2",
//    "Ny2",
//    "Nx_Ny",

    //
    "NfPb_Nf",
    "NfPb_Pb",
    "NfPb_Nf2",
    "NfPb_Pb2",
    "NfPb_Nf_Pb",

    "PfNb_Pf",
    "PfNb_Nb",
    "PfNb_Pf2",
    "PfNb_Nb2",
    "PfNb_Pf_Nb",

    "PfPb_Pf",
    "PfPb_Pb",
    "PfPb_Pf2",
    "PfPb_Pb2",
    "PfPb_Pf_Pb",



//    "NxPy_Nx",
//    "NxPy_Py",
//    "NxPy_Nx2",
//    "NxPy_Py2",
//    "NxPy_Nx_Py",

//    "PxNy_Px",
//    "PxNy_Ny",
//    "PxNy_Px2",
//    "PxNy_Ny2",
//    "PxNy_Px_Ny",

//    "PxPy_Px",
//    "PxPy_Py",
//    "PxPy_Px2",
//    "PxPy_Py2",
//    "PxPy_Px_Py",


    // ##### mixes:
//    "Nf_Nx",
//    "Nf_Ny",
//    "Nb_Nx",
//    "Nb_Ny",

//    "Nf_Px",
//    "Nf_Py",
//    "Nb_Px",
//    "Nb_Py",

//    "Pf_Nx",
//    "Pf_Ny",
//    "Pb_Nx",
//    "Pb_Ny",

//    "Pf_Px",
//    "Pf_Py",
//    "Pb_Px",
//    "Pb_Py",


    // RATIOS:
//    "Nf_OVER_Nb",
//    "Nf_OVER_Nx",
//    "Nf_OVER_Ny",

//    "Nb_OVER_Nf",
//    "Nb_OVER_Nx",
//    "Nb_OVER_Ny",

//    "Nx_OVER_Nf",
//    "Nx_OVER_Nb",
//    "Nx_OVER_Ny",

//    "Ny_OVER_Nf",
//    "Ny_OVER_Nb",
//    "Ny_OVER_Nx",

//    // RATIO - RATIO:
//    "Nf_OVER_Nb_vs_Nx_OVER_Ny",
//    "Nf_OVER_Nx_vs_Nb_OVER_Ny",

//    "Nb_OVER_Nf_vs_Ny_OVER_Nx",
//    "Nx_OVER_Nf_vs_Ny_OVER_Nb",

//    // RATIO - pT:
//    "Nf_OVER_Nx_vs_Pb",
//    "Nf_OVER_Nx_vs_Py",
//    "Nb_OVER_Ny_vs_Pf",
//    "Nb_OVER_Ny_vs_Px",

//    "Nx_OVER_Nf_vs_Pb",
//    "Nx_OVER_Nf_vs_Py",
//    "Ny_OVER_Nb_vs_Pf",
//    "Ny_OVER_Nb_vs_Px",



    // added on 04.04.2019: for dptdpt etc (from /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/AliForwardBackwardAnalysis.cxx:1581)
    "sumPtAllEvF"  ,
    "sumPtAllEvB"  ,
    "piFpjB"       ,
    "nF*PB"    ,
    "nB*PF"    ,
    "pipjF"        ,
    "(nF-1)*PF",
    "nF*(nF-1)"    ,
    "pipjB"        ,
    "(nB-1)*PB",
    "nB*(nB-1)"    ,


    // to check VV comparison ptpt vs dptdpt: special terms:
    "PF*PB",
//    "nB*PF", // already have it above!
//    "nF*PB", // already have it above!
    "PF",
    "PB",

    "sum_pT2_F",
    "sum_pT2_B",

};




/*
const char *varNames[] = {
    "Nevents",

    "f_Nevents",
    "b_Nevents",

//    "x_Nevents",
//    "y_Nevents",

    "fb_Nevents",
//    "xy_Nevents",

//    "fx_Nevents",
//    "fy_Nevents",
//    "bx_Nevents",
//    "by_Nevents",

    // f, b
    "Nf",
    "Nb",
    "N2_f",
    "N2_b",
    "Nf_Nb",

//    // x, y
//    "Nx",
//    "Ny",
//    "Nx2",
//    "Ny2",
//    "Nx_Ny",

    //
    "NfPb_Nf",
    "NfPb_Pb",
    "NfPb_Nf2",
    "NfPb_Pb2",
    "NfPb_Nf_Pb",

    "PfNb_Pf",
    "PfNb_Nb",
    "PfNb_Pf2",
    "PfNb_Nb2",
    "PfNb_Pf_Nb",

    "PfPb_Pf",
    "PfPb_Pb",
    "PfPb_Pf2",
    "PfPb_Pb2",
    "PfPb_Pf_Pb",



//    "NxPy_Nx",
//    "NxPy_Py",
//    "NxPy_Nx2",
//    "NxPy_Py2",
//    "NxPy_Nx_Py",

//    "PxNy_Px",
//    "PxNy_Ny",
//    "PxNy_Px2",
//    "PxNy_Ny2",
//    "PxNy_Px_Ny",

//    "PxPy_Px",
//    "PxPy_Py",
//    "PxPy_Px2",
//    "PxPy_Py2",
//    "PxPy_Px_Py",


    // ##### mixes:
//    "Nf_Nx",
//    "Nf_Ny",
//    "Nb_Nx",
//    "Nb_Ny",

//    "Nf_Px",
//    "Nf_Py",
//    "Nb_Px",
//    "Nb_Py",

//    "Pf_Nx",
//    "Pf_Ny",
//    "Pb_Nx",
//    "Pb_Ny",

//    "Pf_Px",
//    "Pf_Py",
//    "Pb_Px",
//    "Pb_Py",


    // RATIOS:
//    "Nf_OVER_Nb",
//    "Nf_OVER_Nx",
//    "Nf_OVER_Ny",

//    "Nb_OVER_Nf",
//    "Nb_OVER_Nx",
//    "Nb_OVER_Ny",

//    "Nx_OVER_Nf",
//    "Nx_OVER_Nb",
//    "Nx_OVER_Ny",

//    "Ny_OVER_Nf",
//    "Ny_OVER_Nb",
//    "Ny_OVER_Nx",

    // RATIO - RATIO:
//    "Nf_OVER_Nb_vs_Nx_OVER_Ny",
//    "Nf_OVER_Nx_vs_Nb_OVER_Ny",

//    "Nb_OVER_Nf_vs_Ny_OVER_Nx",
//    "Nx_OVER_Nf_vs_Ny_OVER_Nb",

    // RATIO - pT:
//    "Nf_OVER_Nx_vs_Pb",
//    "Nf_OVER_Nx_vs_Py",
//    "Nb_OVER_Ny_vs_Pf",
//    "Nb_OVER_Ny_vs_Px",

//    "Nx_OVER_Nf_vs_Pb",
//    "Nx_OVER_Nf_vs_Py",
//    "Ny_OVER_Nb_vs_Pf",
//    "Ny_OVER_Nb_vs_Px",


    // added on 04.04.2019: for dptdpt etc (from /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/AliForwardBackwardAnalysis.cxx:1581)
    "sumPtAllEvF"  ,
    "sumPtAllEvB"  ,
    "piFpjB"       ,
    "nF*PB"    ,
    "nB*PF"    ,
    "pipjF"        ,
    "(nF-1)*PF",
    "nF*(nF-1)"    ,
    "pipjB"        ,
    "(nB-1)*PB",
    "nB*(nB-1)"    ,


    // to check VV comparison ptpt vs dptdpt: special terms:
    "PF*PB",
//    "nB*PF", // already have it above!
//    "nF*PB", // already have it above!
    "PF",
    "PB",

};
*/

const int nVars = sizeof(varNames)/sizeof(*varNames);




// ###############





// ###############
class WinPairFilling : public WinPairBase //WinPair2019
{
public:
    // e-by-e data:
//    double _nF; double _nB;
    int _nF; int _nB;
    double _ptF; double _ptB;
    double _ptF2; double _ptB2;  // sums of pi^2 for each particle - to check expression for variance

//    double _nX; double _nY;
//    int _nX; int _nY;
//    double _ptX; double _ptY;


    // extra info:
    Double_t pipjF; //, sum_piF;
    Double_t pipjB; //, sum_piB;
    // to calc C:
      double *arr_pF; //[FB_max_tracks_in_win];
      double *arr_pB; //[FB_max_tracks_in_win];
//      float arr_pF[ 10 ];//1200 ];
//      float arr_pB[ 10 ];//1200 ];
    int maxNtracks;


    // for dptdpt:
    Double_t piFpjB;

    TH1D *histAccumulatedValues;  //!    accumulated values for observables



    WinPairFilling():
        _nF(0), _nB(0),
        _ptF(0), _ptB(0),
        _ptF2(0), _ptB2(0),
//        _nX(0), _nY(0),
//        _ptX(0), _ptY(0),
        histAccumulatedValues(0x0)

    {}

    void init(int _nMaxTracks) //setMaxNumberOfTracks(int _n)
    {
          arr_pF = new double[_nMaxTracks];
          arr_pB = new double[_nMaxTracks];

        if(0)cout << "initializing window: " << strAccumHistName.Data() << endl;
        histAccumulatedValues = new TH1D( strAccumHistName //"histAccumulatedValues"
                                          , strAccumHistName, nVars,-0.5,nVars-0.5);

//            TString gArrayMemberNames[nVars];
        for( int i=0; i < nVars; i++ )
            histAccumulatedValues->GetXaxis()->SetBinLabel( i+1, varNames[i] );
    }



    void addTrack( int pid, double eta, double phi, double pt, int charge )
    {
//        cout << "test1" << endl;
        int pidAbs = abs(pid);
//        int charge = (pid!=0 ? pid / pidAbs : 0);

//        cout << pidAbs << " " << charge << endl;

//        cout << "test11" << endl;

        if ( pt < ptWin[0] || pt > ptWin[1] )
            return;

        // backward window:
        if ( eta > eWin[0] && eta < eWin[1]
             && ( phiWin[0]<-999 ? true : ( phi > phiWin[0] && phi < phiWin[1] ) ) // true if NO PHI WINS, else -> check phi win
             )
        {
            // B
            if ( partTypes[1]==0 || pidAbs == partTypes[1] )
                if ( partCharges[1]==0 || charge == partCharges[1] )
                {
                      arr_pB[_nB] = pt;
                    _nB++;
                    _ptB += pt;
                    _ptB2 += pt*pt;
                }
            // Y
//            if ( partTypes[3]==0 || pidAbs == partTypes[3] )
//                if ( partCharges[3]==0 || charge == partCharges[3] )
//                {
//                    _nY++;
//                    _ptY += pt;
//                }
        }
//        cout << "test12" << endl;

        // forward window:
        if ( eta > eWin[2] && eta < eWin[3]
             && ( phiWin[2]<-999 ? true : ( phi > phiWin[2] && phi < phiWin[3] ) ) // true if NO PHI WINS, else -> check phi win
             )
        {
            // F
            if ( partTypes[0]==0 || pidAbs == partTypes[0] )
                if ( partCharges[0]==0 || charge == partCharges[0] )
                {
                      arr_pF[_nF] = pt;
                    _nF++;
                    _ptF += pt;
                    _ptF2 += pt*pt;
                }
            // X
//            if ( partTypes[2]==0 || pidAbs == partTypes[2] )
//                if ( partCharges[2]==0 || charge == partCharges[2] )
//                {
//                    _nX++;
//                    _ptX += pt;
//                }
        }

//        cout << "test2" << endl;

    }

    void finishEvent()
    {
        double meanPtF = -1;
        double meanPtB = -1;
//        double meanPtX = -1;
//        double meanPtY = -1;

//        double rBtoF = -1;
//        double rXtoF = -1;
//        double rYtoF = -1;

//        double rFtoB = -1;
//        double rXtoB = -1;
//        double rYtoB = -1;

//        double rFtoX = -1;
//        double rBtoX = -1;
//        double rYtoX = -1;

//        double rFtoY = -1;
//        double rBtoY = -1;
//        double rXtoY = -1;

//        if ( _nF > 0 ) { meanPtF = _ptF/_nF;    rBtoF=_nB/_nF;  rXtoF=_nX/_nF;  rYtoF=_nY/_nF;  }
//        if ( _nB > 0 ) { meanPtB = _ptB/_nB;    rFtoB=_nF/_nB;  rXtoB=_nX/_nB;  rYtoB=_nY/_nB;  }
//        if ( _nX > 0 ) { meanPtX = _ptX/_nX;    rFtoX=_nF/_nX;  rBtoX=_nB/_nX;  rYtoX=_nY/_nX;  }
//        if ( _nY > 0 ) { meanPtY = _ptY/_nY;    rFtoY=_nF/_nY;  rBtoY=_nB/_nY;  rXtoY=_nX/_nY;  }
        if ( _nF > 0 ) { meanPtF = _ptF/_nF; }
        if ( _nB > 0 ) { meanPtB = _ptB/_nB; }
//        if ( _nX > 0 ) { meanPtX = _ptX/_nX; }
//        if ( _nY > 0 ) { meanPtY = _ptY/_nY; }

        // !!! .. i.e. we keep ptF etc as SUM pT in each event, not mean pT!


        // to calc C:
        // F
          for( Int_t i = 0; i < _nF; i++ )
              for( Int_t j = i+1; j < _nF; j++ )
                  pipjF += arr_pF[i]*arr_pF[j];

          // B
          for( Int_t i = 0; i < _nB; i++ )
              for( Int_t j = i+1; j < _nB; j++ )
                  pipjB += arr_pB[i]*arr_pB[j];

          // for dpt-dpt
          for( Int_t i = 0; i < _nF; i++ )
              for( Int_t j = 0; j < _nB; j++ )
                  piFpjB += arr_pF[i]*arr_pB[j];




        histAccumulatedValues->Fill( "Nevents",   1                  );

        // NfNb:
        histAccumulatedValues->Fill( "Nf"     ,   _nF             );
        histAccumulatedValues->Fill( "Nb"     ,   _nB             );
        histAccumulatedValues->Fill( "N2_f"   ,   _nF*_nF      );
        histAccumulatedValues->Fill( "N2_b"   ,   _nB*_nB      );
        histAccumulatedValues->Fill( "Nf_Nb"  ,   _nF*_nB      );

        // NxNy:
//        histAccumulatedValues->Fill( "Nx"     ,   _nX             );
//        histAccumulatedValues->Fill( "Ny"     ,   _nY             );
//        histAccumulatedValues->Fill( "Nx2"   ,   _nX*_nX      );
//        histAccumulatedValues->Fill( "Ny2"   ,   _nY*_nY      );
//        histAccumulatedValues->Fill( "Nx_Ny"  ,   _nX*_nY      );


        // Nf_Nx Nf_Ny Nb_Nx Nb_Ny:
//        histAccumulatedValues->Fill( "Nf_Nx",   _nF*_nX       );
//        histAccumulatedValues->Fill( "Nf_Ny",   _nF*_nY       );
//        histAccumulatedValues->Fill( "Nb_Nx",   _nB*_nX       );
//        histAccumulatedValues->Fill( "Nb_Ny",   _nB*_nY       );


        // FROM /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/AliForwardBackwardAnalysis.cxx:1583
        histAccumulatedValues->Fill( "sumPtAllEvF", _ptF );
        histAccumulatedValues->Fill( "sumPtAllEvB", _ptB );

        // for dptdpt:
        histAccumulatedValues->Fill( "piFpjB",     piFpjB );
        histAccumulatedValues->Fill( "nF*PB",  _nF*_ptB );
        histAccumulatedValues->Fill( "nB*PF",  _nB*_ptF );

        // for C when av over pairs is OUTSIDE sum:
        histAccumulatedValues->Fill( "pipjF",              pipjF );
        histAccumulatedValues->Fill( "(nF-1)*PF",      (_nF-1)*_ptF ); // (n-1)*sumPt
        histAccumulatedValues->Fill( "nF*(nF-1)",          (_nF-1)*_nF ); // (n-1)*n

        histAccumulatedValues->Fill( "pipjB",         pipjB );
        histAccumulatedValues->Fill( "(nB-1)*PB", (_nB-1)*_ptB ); // (n-1)*sumPt
        histAccumulatedValues->Fill( "nB*(nB-1)",     (_nB-1)*_nB ); // (n-1)*n

        // to check varF vs C_F (IA, Nov 2019):
        histAccumulatedValues->Fill( "sum_pT2_F",     _ptF2 ); // (n-1)*n
        histAccumulatedValues->Fill( "sum_pT2_B",     _ptB2 ); // (n-1)*n


        // to check VV comparison ptpt vs dptdpt: special terms:
        histAccumulatedValues->Fill( "PF*PB",      _ptF*_ptB );
//        histAccumulatedValues->Fill( "nB*PF",            _nB*_ptF );
//        histAccumulatedValues->Fill( "nF*PB",            _nF*_ptB );  //_nF*_nB*_ptB );
        histAccumulatedValues->Fill( "PF",              _ptF );
        histAccumulatedValues->Fill( "PB",              _ptB );


        // NfPb:
        if ( _nB > 0 )
        {
            histAccumulatedValues->Fill( "b_Nevents",    1                  );

            histAccumulatedValues->Fill( "NfPb_Nf",        _nF             );
            histAccumulatedValues->Fill( "NfPb_Pb",       meanPtB             );
            histAccumulatedValues->Fill( "NfPb_Nf2",      _nF*_nF      );
            histAccumulatedValues->Fill( "NfPb_Pb2",      meanPtB*meanPtB      );
            histAccumulatedValues->Fill( "NfPb_Nf_Pb",    _nF*meanPtB      );

//            histAccumulatedValues->Fill( "Nf_OVER_Nb",    _nF/(double)_nB      );
//            histAccumulatedValues->Fill( "Nx_OVER_Nb",    _nX/(double)_nB      );
//            histAccumulatedValues->Fill( "Ny_OVER_Nb",    _nY/(double)_nB      );

//            histAccumulatedValues->Fill( "Pb_Nx",   meanPtB*_nX      );
//            histAccumulatedValues->Fill( "Pb_Ny",   meanPtB*_nY      );


//            // for C when av over pairs is OUTSIDE sum:
//            histAccumulatedValues->Fill( "pipjB",         pipjB );
//            histAccumulatedValues->Fill( "(nB-1)*sum_pB", (_nB-1)*(_nB*_ptB) ); // (n-1)*sumPt
//            histAccumulatedValues->Fill( "nB*(nB-1)",     (_nB-1)*_nB ); // (n-1)*n

            // for C when av over pairs is INSIDE sum: (like in GOOD_Ck_definition_STAR_2005_0504031.pdf)
            int nPairsB = (_nB-1)*_nB;
            if ( nPairsB > 0 )
            {
//                histAccumulatedValues->Fill( "pipjB_avPerEv",         pipjB / nPairsB );
//                histAccumulatedValues->Fill( "(nB-1)*PB_avPerEv", (_nB-1)*_ptB / nPairsB );
            }

        }

        // PfNb:
        if ( _nF > 0 )
        {
            histAccumulatedValues->Fill( "f_Nevents",    1                  );

            histAccumulatedValues->Fill( "PfNb_Pf",         meanPtF           );
            histAccumulatedValues->Fill( "PfNb_Nb",        _nB            );
            histAccumulatedValues->Fill( "PfNb_Pf2",     meanPtF*meanPtF      );
            histAccumulatedValues->Fill( "PfNb_Nb2",      _nB*_nB      );
            histAccumulatedValues->Fill( "PfNb_Pf_Nb",    meanPtF*_nB      );

//            histAccumulatedValues->Fill( "Nb_OVER_Nf",    _nB/(double)_nF      );
//            histAccumulatedValues->Fill( "Nx_OVER_Nf",    _nX/(double)_nF      );
//            histAccumulatedValues->Fill( "Ny_OVER_Nf",    _nY/(double)_nF      );

//            histAccumulatedValues->Fill( "Pf_Nx",   meanPtF*_nX      );
//            histAccumulatedValues->Fill( "Pf_Ny",   meanPtF*_nY      );

//            // for C when av over pairs is OUTSIDE sum:
//            histAccumulatedValues->Fill( "pipjF",              pipjF );
//            histAccumulatedValues->Fill( "(nF-1)*sum_pF",      (_nF-1)*(_nF*_ptF) ); // (n-1)*sumPt
//            histAccumulatedValues->Fill( "nF*(nF-1)",          (_nF-1)*_nF ); // (n-1)*n

            // for C when av over pairs is INSIDE sum:
            int nPairsF = (_nF-1)*_nF;
            if ( nPairsF > 0 )
            {
//                histAccumulatedValues->Fill( "pipjF_avPerEv",          pipjF / nPairsF );
//                histAccumulatedValues->Fill( "(nF-1)*PF_avPerEv",  (_nF-1)*_ptF / nPairsF );
            }

        }

        // fb:
        if ( _nF > 0 && _nB > 0 )
        {
            histAccumulatedValues->Fill( "fb_Nevents",   1                  );

            histAccumulatedValues->Fill( "PfPb_Pf",       meanPtF             );
            histAccumulatedValues->Fill( "PfPb_Pb",       meanPtB             );
            histAccumulatedValues->Fill( "PfPb_Pf2",     meanPtF*meanPtF      );
            histAccumulatedValues->Fill( "PfPb_Pb2",     meanPtB*meanPtB      );
            histAccumulatedValues->Fill( "PfPb_Pf_Pb",   meanPtF*meanPtB      );

//            histAccumulatedValues->Fill( "Nx_OVER_Nf_vs_Pb",   _nX/(double)_nF *meanPtB      );
//            histAccumulatedValues->Fill( "Ny_OVER_Nb_vs_Pf",   _nY/(double)_nB *meanPtF      );

//            histAccumulatedValues->Fill( "Nx_OVER_Nf_vs_Ny_OVER_Nb",   _nX/(double)_nF * _nY/(double)_nB      );

//            // for dptdpt:
//            histAccumulatedValues->Fill( "piFpjB",     piFpjB );
//            histAccumulatedValues->Fill( "nF*sum_pB",  _nF*(_nB*_ptB) );
//            histAccumulatedValues->Fill( "nB*sum_pF",  _nB*(_nF*_ptF) );

        }


        // y:
//        if ( _nY > 0 )
//        {
//            histAccumulatedValues->Fill( "y_Nevents",    1                  );

//            histAccumulatedValues->Fill( "NxPy_Nx",        _nX             );
//            histAccumulatedValues->Fill( "NxPy_Py",       meanPtY             );
//            histAccumulatedValues->Fill( "NxPy_Nx2",      _nX*_nX      );
//            histAccumulatedValues->Fill( "NxPy_Py2",      meanPtY*meanPtY      );
//            histAccumulatedValues->Fill( "NxPy_Nx_Py",    _nX*meanPtY      );


//            histAccumulatedValues->Fill( "Nf_OVER_Ny",    _nF/(double)_nY      );
//            histAccumulatedValues->Fill( "Nb_OVER_Ny",    _nB/(double)_nY      );
//            histAccumulatedValues->Fill( "Nx_OVER_Ny",    _nX/(double)_nY      );

//            histAccumulatedValues->Fill( "Nf_Py",   _nF*meanPtY       );
//            histAccumulatedValues->Fill( "Nb_Py",   _nB*meanPtY       );
//        }

//        // x:
//        if ( _nX > 0 )
//        {
//            histAccumulatedValues->Fill( "x_Nevents",    1                  );

//            histAccumulatedValues->Fill( "PxNy_Px",         meanPtX           );
//            histAccumulatedValues->Fill( "PxNy_Ny",        _nY            );
//            histAccumulatedValues->Fill( "PxNy_Px2",     meanPtX*meanPtX      );
//            histAccumulatedValues->Fill( "PxNy_Ny2",      _nY*_nY      );
//            histAccumulatedValues->Fill( "PxNy_Px_Ny",    meanPtX*_nY      );

//            histAccumulatedValues->Fill( "Nf_OVER_Nx",    _nF/(double)_nX      );
//            histAccumulatedValues->Fill( "Nb_OVER_Nx",    _nB/(double)_nX      );
//            histAccumulatedValues->Fill( "Ny_OVER_Nx",    _nY/(double)_nX      );

//            histAccumulatedValues->Fill( "Nf_Px",   _nF*meanPtX       );
//            histAccumulatedValues->Fill( "Nb_Px",   _nB*meanPtX       );
//        }


//        // xy:
//        if ( _nX > 0 && _nY > 0 )
//        {
//            histAccumulatedValues->Fill( "xy_Nevents",   1                  );
//            histAccumulatedValues->Fill( "PxPy_Px",       meanPtX             );
//            histAccumulatedValues->Fill( "PxPy_Py",       meanPtY             );
//            histAccumulatedValues->Fill( "PxPy_Px2",     meanPtX*meanPtX      );
//            histAccumulatedValues->Fill( "PxPy_Py2",     meanPtY*meanPtY      );
//            histAccumulatedValues->Fill( "PxPy_Px_Py",   meanPtX*meanPtY      );

//            histAccumulatedValues->Fill( "Nf_OVER_Nx_vs_Py",   _nF/(double)_nX*meanPtY      );
//            histAccumulatedValues->Fill( "Nb_OVER_Ny_vs_Px",   _nB/(double)_nY*meanPtX      );

//            histAccumulatedValues->Fill( "Nf_OVER_Nx_vs_Nb_OVER_Ny",   _nF/(double)_nX * _nB/(double)_nY      );
//        }

//        // fx:
//        if ( _nF > 0 && _nX > 0 )
//        {
//            histAccumulatedValues->Fill( "fx_Nevents",   1                  );

//            histAccumulatedValues->Fill( "Pf_Px",   meanPtF*meanPtX      );

//            histAccumulatedValues->Fill( "Nb_OVER_Nf_vs_Ny_OVER_Nx",   _nB/(double)_nF * _nY/(double)_nX      );
//        }
//        // fy:
//        if ( _nF > 0 && _nY > 0 )
//        {
//            histAccumulatedValues->Fill( "fy_Nevents",   1                  );

//            histAccumulatedValues->Fill( "Pf_Py",   meanPtF*meanPtY      );
//            histAccumulatedValues->Fill( "Nx_OVER_Nf_vs_Py",   _nX/(double)_nF *meanPtY      );
//            histAccumulatedValues->Fill( "Nb_OVER_Ny_vs_Pf",   _nB/(double)_nY *meanPtF      );
//        }
//        // bx:
//        if ( _nB > 0 && _nX > 0 )
//        {
//            histAccumulatedValues->Fill( "bx_Nevents",   1                  );

//            histAccumulatedValues->Fill( "Pb_Px",   meanPtB*meanPtX      );
//            histAccumulatedValues->Fill( "Nf_OVER_Nx_vs_Pb",   _nF/(double)_nX *meanPtB      );
//            histAccumulatedValues->Fill( "Ny_OVER_Nb_vs_Px",   _nY/(double)_nB *meanPtX      );
//        }
//        // by:
//        if ( _nB > 0 && _nY > 0 )
//        {
//            histAccumulatedValues->Fill( "by_Nevents",   1                  );
//            histAccumulatedValues->Fill( "Pb_Py",   meanPtB*meanPtY      );

//            histAccumulatedValues->Fill( "Nf_OVER_Nb_vs_Nx_OVER_Ny",   _nF/(double)_nB * _nX/(double)_nY      );
//        }


        // reset values:
        _nF = 0; _nB = 0;
        _ptF = 0; _ptB = 0;
        _ptF2 = 0; _ptB2 = 0;

//        _nX = 0; _nY = 0;
//        _ptX = 0; _ptY = 0;


        // extra info:
        pipjF = 0; //sum_piF = 0;
        pipjB = 0; //sum_piB = 0;

        piFpjB = 0;

    }

//    ClassDef(WinPairFilling, 1);


};


