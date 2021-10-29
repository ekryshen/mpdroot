// =======================================================================================================================================
// Supplementary class for Forward-backward correlations analysis
// Author: Igor Altsybeev (SPbSU), 10.2020
// =======================================================================================================================================

#include "TH1.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TString.h"

#include <iostream>
#include <map>
using namespace std;


// #############
struct WinPairBase
{
public:
    // acceptance of the windows:
    double eWin[4]; // eta for B (min, max), F (min, max),
    double phiWin[4]; // phi for B (min, max), F (min, max),
    double ptWin[2];
    int partTypes[4];   // F,B, X,Y  // 0 - NO PID SELECTION
    int partCharges[4]; // F,B, X,Y  // 0 - NO CHARGE SELECTION

    int subsampleId;
    // main hist for data accumulation:
//    TH1D *histAccumulatedValues;   //! accumulated values for observables

    TString strAccumHistName;

    WinPairBase() // int _subsId = -1)
    {
        subsampleId = 0;
//        cout << "constructor WinPairBase" << endl;
        for( int i = 0; i < 4; i++ )
        {
            partTypes[i] = 0;
            partCharges[i] = 0;
        }


    }

    ~WinPairBase()
    {
    }

    void setSubsampleId( int id )
    {
        subsampleId = id;
    }

    void setParticleTypes( int *_pTypes, int *_pCharges )
    {
        for( int i = 0; i < 4; i++ )
        {
            partTypes[i] = _pTypes[i];
            partCharges[i] = _pCharges[i];
        }
    }

    void setWindows( const char* strPrefix,
                     int iCW, int cBin, double _eMinB, double _eMaxB, double _eMinF, double _eMaxF, double _ptMin, double _ptMax,
                     double _phiMinB=-1000, double _phiMaxB=-1000, double _phiMinF=-1000, double _phiMaxF=-1000 )
    {
        // fill map with name-number pairs
//        for( int i=0; i < nVars; i++ )
//            mapVar.insert( pair<const char*,int>(varNames[i], i) );

        eWin[0] = _eMinB;
        eWin[1] = _eMaxB;
        eWin[2] = _eMinF;
        eWin[3] = _eMaxF;

        phiWin[0] = _phiMinB;
        phiWin[1] = _phiMaxB;
        phiWin[2] = _phiMinF;
        phiWin[3] = _phiMaxF;



        ptWin[0] = _ptMin;
        ptWin[1] = _ptMax;

        TString strPID = Form( "pid_%d_%d_%d_%d_charge_%d_%d_%d_%d",
                               partTypes[0], partTypes[1], partTypes[2], partTypes[3],
                partCharges[0], partCharges[1], partCharges[2], partCharges[3]
                );

        //        cout << "dir = " << dir << endl;
        strAccumHistName = Form( "histAccumulatedValues_%s_%s_cW%d_cBin%d_eta_B_%.1f_%.1f_F_%.1f_%.1f_pt_%.1f_%.1f",
                                strPrefix, strPID.Data(), iCW, cBin, eWin[0], eWin[1], eWin[2], eWin[3], ptWin[0], ptWin[1]
//                strAccumHistName = Form( "histAccumulatedValues_%s_cW%d_cBin%d_eta_B_%.1f_%.1f_F_%.1f_%.1f_pt_%.1f_%.1f",
//                                        strPID.Data(), iCW, cBin, eWin[0], eWin[1], eWin[2], eWin[3], ptWin[0], ptWin[1]
                );
        if( phiWin[0] > -999 )
        {
            TString strPhiWin = Form( "phi_B_%.1f_%.1f_F_%.1f_%.1f", phiWin[0], phiWin[1], phiWin[2], phiWin[3] );
            strAccumHistName.Append( Form( "_%s", strPhiWin.Data() ) );

        }

        if ( subsampleId >=0 ) // i.e. work with subsamples
            strAccumHistName.Append( Form( "_subs%d", subsampleId ) );


//    ClassDef(WinPairBase, 1);


};


