// =======================================================================================================================================
// This analysis code does several things:
// 1) it performs the analysis of Forward-backward correlations at SIM, REC_PRIMARY and REC_ALL levels
// 2) it extracts a special tree with tracks for Identity Method (see fTreeIM below)
// 3) a lot of supplementary information is extracted, in particular, various DCA distributions.
// Author: Igor Altsybeev (SPbSU), 10.2020
// =======================================================================================================================================

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THn.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"

#include "TChain.h"

#include "TStyle.h"
#include "TProfile.h"

#include "TStopwatch.h"
#include "TLine.h"

#include "TRandom3.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TMath.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TClonesArray.h"

#include "TDatabasePDG.h"

TDatabasePDG *db_PDG = new TDatabasePDG;

#include "WinPairFilling.h"


// for FB analysis:
const int nWinPairWithSubsamples = 15;//3;//15;//8;//30;

double eRange = 1.0;//1.4;//1.6;//0.8;
double eSize = 0.2;//4;
double eStep = 0.2;

//const int nPhiWins = 16;
//const int nEtaBins = /*2**/nPhiWins;//23; //(eRange-eSize) / eStep + 2;//1;
const int nEtaBins = (eRange-eSize) / eStep + 1 + 0.0000001;//2;//1;
//double phiStep = TMath::TwoPi()/nPhiWins;

// FOR ITS+MFT:


//    WinPair2019 *winPairs_REC[10];
const int nPtBins = 1;
double ptmin[nPtBins] = { 0.2 };//0.1,  };//0.2 }; //0.3, };//0.2 };
double ptmax[nPtBins] = { 2.0,};// 2.0 }; //1.5,};// 5.0 };
//    double ptmin = 0.2;
//    double ptmax = 5.0;


//
const int nPartTypes = 3;//4;
int arrPartTypes[nPartTypes][4] =
{ // F,B, X,Y
  { /*421*/0, 0, 0, 0 },         // 0
  { 0, 0, 0, 0 },         // 1
  { 0, 0, 0, 0 },         // 2
//  { 321, 321, 211, 211 }, // 6

  //      { 321, 321, 211, 211 }, // 7
  //      { 321, 321, 211, 211 }, // 8
  //      { 321, 321, 211, 211 }, // 9
};

int arrCharges[nPartTypes][4] =
{ // F,B, X,Y
  {  0,  0,  0,  0 }, // 0
  { +1, +1, +1, +1 }, // 1
  { +1, -1, +1, -1 }, // 2
//  {  0,  0,  0,  0 }, // 0
};


const int nCW = 1;
const int nCentrBins[] = { 10 };//9 };
const double cBinWidths[] = { 100 }; // in %


WinPairFilling winPairs_SIM[nWinPairWithSubsamples][nPartTypes][nCW][ 10/*nCentrBins[nCW-1]*/ ][nPtBins][10/*nEtaBins*/];
WinPairFilling winPairs_REC[nWinPairWithSubsamples][nPartTypes][nCW][ 10/*nCentrBins[nCW-1]*/ ][nPtBins][10/*nEtaBins*/];
WinPairFilling winPairs_REC_PRIM[nWinPairWithSubsamples][nPartTypes][nCW][ 10/*nCentrBins[nCW-1]*/ ][nPtBins][10/*nEtaBins*/];





// ### DCA histos:
const int N_MAX_TYPES_FOR_DCA = 10;

//
TH2D *fHist2D_vZ_vs_DCAz;     //!
//
TH1D *fHist1D_DCA_XY[N_MAX_TYPES_FOR_DCA];        //!
TH1D *fHist1D_DCA_Z[N_MAX_TYPES_FOR_DCA];         //!
TH1D *fHist1D_DCA_Z_SPEC_QA[N_MAX_TYPES_FOR_DCA]; //!
TH2D *fHist2D_DCA_XYvsZ[N_MAX_TYPES_FOR_DCA];     //!

TH1D *fHist1D_DCA_XYZ[N_MAX_TYPES_FOR_DCA];          //!
TH1D *fHist1D_DCA_XYZ_Artif_FB_Sign[N_MAX_TYPES_FOR_DCA];          //!

// vs pt
TH2D *fHist2D_DCA_XY_vs_pt[N_MAX_TYPES_FOR_DCA];          //!
TH2D *fHist2D_DCA_Z_vs_pt[N_MAX_TYPES_FOR_DCA];           //!
TH2D *fHist2D_DCA_Z_SPEC_QA_vs_pt[N_MAX_TYPES_FOR_DCA];   //!
TH2D *fHist2D_DCA_XYZ_vs_pt[N_MAX_TYPES_FOR_DCA];          //!
TH2D *fHist2D_DCA_XYZ_vs_pt_Artif_FB_Sign[N_MAX_TYPES_FOR_DCA];          //!
TH2D *fHist2D_DCA_XYZ_vs_pt_FB_02_06_Artif_FB_Sign[N_MAX_TYPES_FOR_DCA];          //!

// vs pt, eta
TH3D *fHist3D_DCA_XY_vs_PtEta[N_MAX_TYPES_FOR_DCA];          //!
TH3D *fHist3D_DCA_Z_vs_PtEta[N_MAX_TYPES_FOR_DCA];           //!
TH3D *fHist3D_DCA_XYZ_vs_PtEta[N_MAX_TYPES_FOR_DCA];          //!
TH3D *fHist3D_DCA_XYZ_vs_PtEta_Artif_FB_Sign[N_MAX_TYPES_FOR_DCA];          //!

TH1D *fHist1D_IsDCAok;        //!





void FillHistosDCA( int id, double dcaXY, double dcaZ, double pt, double eta )
{
    double dca_r = sqrt(dcaXY*dcaXY + dcaZ*dcaZ);

    fHist1D_DCA_XY[id]          ->Fill( dcaXY );
    fHist1D_DCA_Z[id]           ->Fill( dcaZ );
    fHist1D_DCA_Z_SPEC_QA[id]   ->Fill( dcaZ );
    fHist2D_DCA_XYvsZ[id]       ->Fill( dcaXY, dcaZ );

    fHist1D_DCA_XYZ[id]           ->Fill( dca_r );
    fHist1D_DCA_XYZ_Artif_FB_Sign[id]           ->Fill( eta > 0 ? dca_r : -dca_r );


    fHist2D_DCA_XY_vs_pt[id]          ->Fill( dcaXY, pt );
    fHist2D_DCA_Z_vs_pt[id]           ->Fill( dcaZ, pt );
    fHist2D_DCA_Z_SPEC_QA_vs_pt[id]   ->Fill( dcaZ, pt );
    fHist2D_DCA_XYZ_vs_pt[id]           ->Fill( dca_r, pt );
    fHist2D_DCA_XYZ_vs_pt_Artif_FB_Sign[id]           ->Fill( eta > 0 ? dca_r : -dca_r, pt );
    if( fabs(eta) > 0.2 && fabs(eta) < 0.6 )
        fHist2D_DCA_XYZ_vs_pt_FB_02_06_Artif_FB_Sign[id]           ->Fill( eta > 0 ? dca_r : -dca_r, pt );


    fHist3D_DCA_XY_vs_PtEta[id]          ->Fill( dcaXY, pt, eta );
    fHist3D_DCA_Z_vs_PtEta[id]           ->Fill( dcaZ, pt, eta );
    fHist3D_DCA_XYZ_vs_PtEta[id]          ->Fill( dca_r, pt, eta );
    fHist3D_DCA_XYZ_vs_PtEta_Artif_FB_Sign[id] ->Fill( eta > 0 ? dca_r : -dca_r, pt, eta );
}





// ####################
void get_info_and_identity_tree()
{
    TChain *fEventTree = new TChain("tree_sim_rec");
    if(0)
    {
        fEventTree->Add( "output_tree_2020_09_02_SMASH_0_3/output_with_tree_fId_0_3.root" );
    }
    else
    {
        int step = 20;
        for ( int iFile = 0; iFile < 3/*20*//*20*/; iFile++)
        {
            int minId = iFile*step;
            int maxId = (iFile+1)*step-1;
            fEventTree->Add( Form("../output_tree_2020_09_02_SMASH_%d_%d/output_with_tree_fId_%d_%d.root", minId, maxId, minId, maxId ) );
        }
    }



//    TTree *fEventTree = (TTree*)file->Get("tree_sim_rec");
    if(!fEventTree)
    {
        cout << "NO SUCH TREE!" << endl;
        return;
    }

    //    fEventTree->Draw( "trackCharge" );
    //    fEventTree->Draw( "trackChargeUsingPt", "", "same" );



    // output file:
//    TString outname = Form( "/zfs/hybrilit.jinr.ru/user/a/altsybee/analysis_IA/output_tree2/analyzed_with_FHCAL_classes/an_results_%s", arr_str_filenames[file_id_to_analyze].Data() );
    TString outname = "an_results.root";
    cout << "output file name: " << outname << endl;




    cout << "nEtaBins = " << nEtaBins << endl;
    cout << "(eRange-eSize) / eStep + 1 = " << (eRange-eSize) / eStep + 1 << endl;
    //    return 0;

    // ### FB:
    for ( int iSub = 0; iSub < nWinPairWithSubsamples; iSub++)
        for ( int iType = 0; iType < nPartTypes; iType++)
            for ( int iCW = 0; iCW < nCW; iCW++)
                for ( int cBin = 0; cBin < nCentrBins[iCW]; ++cBin )
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                        for ( int iEta = 0; iEta < nEtaBins; ++iEta )
                        {
                            // SIM
                            WinPairFilling *wp = &winPairs_SIM[iSub][iType][iCW][cBin][iPt][iEta];
                            //wp = new WinPairFilling();
                            wp->setSubsampleId( iSub );
                            wp->setParticleTypes( arrPartTypes[iType], arrCharges[iType] );
                            wp->setWindows(  "SIM", iCW, cBin, -eRange+iEta*eStep, -eRange+eSize+iEta*eStep, eRange-eSize-iEta*eStep, eRange-iEta*eStep, ptmin[iPt], ptmax[iPt] );
                            wp->init( 1500 );

                            // REC
                            wp = &winPairs_REC[iSub][iType][iCW][cBin][iPt][iEta];
                            //wp = new WinPairFilling();
                            wp->setSubsampleId( iSub );
                            wp->setParticleTypes( arrPartTypes[iType], arrCharges[iType] );
                            wp->setWindows(  "REC", iCW, cBin, -eRange+iEta*eStep, -eRange+eSize+iEta*eStep, eRange-eSize-iEta*eStep, eRange-iEta*eStep, ptmin[iPt], ptmax[iPt] );
                            wp->init( 1500 );

                            // REC PRIM
                            wp = &winPairs_REC_PRIM[iSub][iType][iCW][cBin][iPt][iEta];
                            //wp = new WinPairFilling();
                            wp->setSubsampleId( iSub );
                            wp->setParticleTypes( arrPartTypes[iType], arrCharges[iType] );
                            wp->setWindows(  "REC_PRIM", iCW, cBin, -eRange+iEta*eStep, -eRange+eSize+iEta*eStep, eRange-eSize-iEta*eStep, eRange-iEta*eStep, ptmin[iPt], ptmax[iPt] );
                            wp->init( 1500 );
                        }


    //    double bImpRanges[] = { 4.69658, 6.61771, 8.08484, 9.32845, 10.4377, 11.427, 12.3386, 13.194, 19.912, };
    double bImpRanges[]   = { 4.54532, 6.42462, 7.8728, 9.08507, 10.161, 11.1334, 12.0291, 12.8691, 13.7706, };
    //    int nCbins = sizeof(bImpRanges)/sizeof(*bImpRanges);

    double energyRanges[] = { 14.8081, 19.1973, 22.0429, 23.7457, 24.7246, 25.4155, 26.0418, 26.7621, 99.6,  };
//    int nCbins = sizeof(energyRanges)/sizeof(*energyRanges);

    double multRanges[] = { 3199.5, 301.076, 216.702, 155.07, 108.548, 73.5102, 47.7274, 29.0818, 16.1795, 8.16926,  };
    //{ 3199.5, 309.307, 219.908, 155.215, 107.217, 71.3813, 45.524, 27.1389, 15.9874, 9.32887, -1 };
    int nCbins = sizeof(multRanges)/sizeof(*multRanges);



    // ##### histos IA:
    TH1D *h_ImpactPar = new TH1D("h_ImpactPar",";b (fm);Entries",2000,-2,20); //440,-2,20);

    TH1D *h_nTracksInEtaCuts_SIM = new TH1D("h_nTracksInEtaCuts_SIM",";n tracks;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);
    TH1D *h_nTracksInEtaCuts_REC = new TH1D("h_nTracksInEtaCuts_REC",";n tracks;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);
    TH1D *h_nTracksInEtaCuts_REC_ACCEPTED = new TH1D("h_nTracksInEtaCuts_REC_ACCEPTED",";n tracks;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);

    TH1D *h_nTracksForMultCentrality = new TH1D("h_nTracksForMultCentrality",";n tracks;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);

    TH1D *h_QA_centr_bins = new TH1D("h_QA_centr_bins",";centrality bin;Entries", 101, -0.5, 100.5 );




    TH1D *h_QA_meanMultInFiles = new TH1D("h_QA_meanMultInFiles",";file id;Entries", 3001, -0.5, 3000.5 );


    int LogBinsX = 581;
    double StartValueX = 0.1;//
    double StopValueX = 20;//
    double LogWidthX[LogBinsX+1];

    int LogBinsY = 581;
    double StartValueY = 1000;//
    double StopValueY = 20000;//
    double LogWidthY[LogBinsY+1];

    //calculate bins
    for(int i = 0; i <= LogBinsX; i++) LogWidthX[i] = pow(10,log10(StartValueX)+(log10(StopValueX)-log10(StartValueX))/double(LogBinsX)*double(i));
    for(int i = 0; i <= LogBinsY; i++) LogWidthY[i] = pow(10,log10(StartValueY)+(log10(StopValueY)-log10(StartValueY))/double(LogBinsY)*double(i));


    TH1D*  h_dEdx = new TH1D("h_dEdx","h_dEdx", LogBinsY, LogWidthY );


    // ### vertex
    TH2D* h_VertexXY = new TH2D("h_VertexXY", "h_VertexXY",200, -1, 1, 200, -1, 1 );
    TH1D* h_VertexZ = new TH1D("h_VertexZ", "h_VertexZ", 800, -80, 80);
    TH1D* h_VertexZ_accepted = new TH1D("h_VertexZ_accepted", "h_VertexZ_accepted", 800, -80, 80);
//    TH1D* h_VertexZ_accepted_excl_0 = new TH1D("h_VertexZ_accepted_excl_0", "h_VertexZ_accepted_excl_0", 800, -80, 80);


    // ### vertexZ spec
    //    TH1D* h1D_vertexZ = new TH1D( "h1D_vertexZ", ";vertex Z;n events", 400, -40, 40  );
    TH1D* h1D_vertexZspec = new TH1D( "h1D_vertexZspec", ";h1D_vertexZspec Z;n events", 800, 22, 25  );
    TH2D* h2D_eta_vs_vZ = new TH2D("h2D_eta_vs_vZ", ";v_{Z} (cm);#eta", 50, -50, 50, 40, -2, 2 );
    TH2D* h2D_eta_vs_vZ_if_nClustAbove30 = new TH2D("h2D_eta_vs_vZ_if_nClustAbove30", ";v_{Z} (cm);#eta", 50, -50, 50, 40, -2, 2 );
    TH2D* h2D_chi2_perCluster_vs_vZ = new TH2D("h2D_chi2_perCluster_vs_vZ", ";v_{Z} (cm);#chi^{2}/cluster", 50, -50, 50, 40, 0, 8 );
    TH2D* h2D_hits_vs_vZ = new TH2D("h2D_hits_vs_vZ", ";v_{Z} (cm);n TPC hits", 50, -50, 50, 56, -0.5, 55.5 );
    TH2D* h2D_hits_vs_vZ_if_etaWithin1 = new TH2D("h2D_hits_vs_vZ_if_etaWithin1", ";v_{Z} (cm);n TPC hits", 50, -50, 50, 56, -0.5, 55.5 );

    TH2D* h2D_hits_vs_eta_if_vZ_above_20 = new TH2D("h2D_hits_vs_eta_if_vZ_above_20", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );
    TH2D* h2D_hits_vs_eta_if_vZ_below_minus20 = new TH2D("h2D_hits_vs_eta_if_vZ_below_minus20", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );
    TH2D* h2D_hits_vs_eta_if_vZ_WITHIN_20 = new TH2D("h2D_hits_vs_eta_if_vZ_WITHIN_20", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );
    TH2D* h2D_hits_vs_eta_if_vZ_WITHIN_10 = new TH2D("h2D_hits_vs_eta_if_vZ_WITHIN_10", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );
    TH2D* h2D_hits_vs_eta_if_vZ_WITHIN_2_excl_0 = new TH2D("h2D_hits_vs_eta_if_vZ_WITHIN_2_excl_0", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );

    // if STRICT DCAz
    TH2D* h2D_STRICT_DCAz_hits_vs_eta_if_vZ_above_20 = new TH2D("h2D_STRICT_DCAz_hits_vs_eta_if_vZ_above_20", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );
    TH2D* h2D_STRICT_DCAz_hits_vs_eta_if_vZ_below_minus20 = new TH2D("h2D_STRICT_DCAz_hits_vs_eta_if_vZ_below_minus20", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );
    TH2D* h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_20 = new TH2D("h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_20", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );
    TH2D* h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_10 = new TH2D("h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_10", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );
    TH2D* h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_2_excl_0 = new TH2D("h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_2_excl_0", ";#eta;n TPC hits", 50, -2.5, 2.5, 56, -0.5, 55.5 );




    // ### track info
    int _nYBins = 800;
    double etaMin = -8;
    double etaMax = 8;
    int _nPtBins = 500;
    double ptMin = 0;
    double ptMax = 5;

    // ### SIM:
    TH1D *h_SIM_Pt = new TH1D("h_SIM_Pt", "h_SIM_Pt", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_SIM_Eta = new TH1D("h_SIM_Eta", "h_SIM_Eta", _nYBins, etaMin, etaMax );
    TH1D *h_SIM_Phi = new TH1D("h_SIM_Phi", "h_SIM_Phi", 400, -2*TMath::TwoPi(), 2*TMath::TwoPi());
    TH1D *h_SIM_Rapidity = new TH1D("h_SIM_Rapidity", "h_SIM_Rapidity", _nYBins, etaMin, etaMax );

    TH1D *h_SIM_Pt_in_eta1 = new TH1D("h_SIM_Pt_in_eta1", "h_SIM_Pt_in_eta1", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_SIM_Eta_in_pt015_2 = new TH1D("h_SIM_Eta_in_pt015_2", "h_SIM_Eta_in_pt015_2", _nYBins, etaMin, etaMax );

    // ### REC:
    TH1D* h_Pt = new TH1D("h_Pt", "h_Pt", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D* h_Eta = new TH1D("h_Eta", "h_Eta", _nYBins, etaMin, etaMax );
    TH1D* h_Phi = new TH1D("h_Phi", "h_Phi", 400, -2*TMath::TwoPi(), 2*TMath::TwoPi());
    TH1D* h_Charge = new TH1D("h_Charge", "h_Charge", 7, -3.5, 3.5);
    //    TH1D* h_Y = new TH1D("h_Y", "h_Y", _nYBins, etaMin, etaMax ); //500, -1, 1);

    TH1D *h_Pt_in_eta1 = new TH1D("h_Pt_in_eta1", "h_Pt_in_eta1", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_Eta_in_pt015_2 = new TH1D("h_Eta_in_pt015_2", "h_Eta_in_pt015_2", _nYBins, etaMin, etaMax );
    TH1D *h_RecPrim_Pt_in_eta1 = new TH1D("h_RecPrim_Pt_in_eta1", "h_RecPrim_Pt_in_eta1", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_RecPrim_Eta_in_pt015_2 = new TH1D("h_RecPrim_Eta_in_pt015_2", "h_RecPrim_Eta_in_pt015_2", _nYBins, etaMin, etaMax );


    TH1D* h1D_chi2_perCluster = new TH1D("h1D_chi2_perCluster", ";#chi^{2}/cluster;n tracks", 200, 0, 20 );
    TH1D* h1D_nClusters = new TH1D("h1D_nClusters", ";n clusters;n tracks", 56, -0.5, 55.5 );


    // ### REC ACCEPTED:
    TH1D* h_ACCEPTED_Pt = new TH1D("h_ACCEPTED_Pt", "h_ACCEPTED_Pt", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D* h_ACCEPTED_Eta = new TH1D("h_ACCEPTED_Eta", "h_ACCEPTED_Eta", _nYBins, etaMin, etaMax );
    TH1D* h_ACCEPTED_Phi = new TH1D("h_ACCEPTED_Phi", "h_ACCEPTED_Phi", 400, -2*TMath::TwoPi(), 2*TMath::TwoPi());
    TH1D* h_ACCEPTED_Charge = new TH1D("h_ACCEPTED_Charge", "h_ACCEPTED_Charge", 7, -3.5, 3.5);
    //    TH1D* h_Y = new TH1D("h_Y", "h_Y", _nYBins, etaMin, etaMax ); //500, -1, 1);

    TH1D *h_ACCEPTED_Pt_in_eta1 = new TH1D("h_ACCEPTED_Pt_in_eta1", "h_ACCEPTED_Pt_in_eta1", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_ACCEPTED_Eta_in_pt015_2 = new TH1D("h_ACCEPTED_Eta_in_pt015_2", "h_ACCEPTED_Eta_in_pt015_2", _nYBins, etaMin, etaMax );
    TH1D *h_ACCEPTED_RecPrim_Pt_in_eta1 = new TH1D("h_ACCEPTED_RecPrim_Pt_in_eta1", "h_ACCEPTED_RecPrim_Pt_in_eta1", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_ACCEPTED_RecPrim_Eta_in_pt015_2 = new TH1D("h_ACCEPTED_RecPrim_Eta_in_pt015_2", "h_ACCEPTED_RecPrim_Eta_in_pt015_2", _nYBins, etaMin, etaMax );

    TH1D* h1D_BEFORE_CUT_chi2_perCluster = new TH1D("h1D_BEFORE_CUT_chi2_perCluster", ";#chi^{2}/cluster;n tracks", 200, 0, 20 );
    TH1D* h1D_ACCEPTED_chi2_perCluster = new TH1D("h1D_ACCEPTED_chi2_perCluster", ";#chi^{2}/cluster;n tracks", 200, 0, 20 );
    TH1D* h1D_ACCEPTED_nClusters = new TH1D("h1D_ACCEPTED_nClusters", ";n clusters;n tracks", 56, -0.5, 55.5 );



    // N SIGMA plots:
    TH1D* h_nSigmaPion   = new TH1D("h_nSigmaPion", "h_nSigmaPion", 200, 0, 10 );
    TH1D* h_nSigmaKaon   = new TH1D("h_nSigmaKaon", "h_nSigmaKaon", 200, 0, 10 );
    TH1D* h_nSigmaProton   = new TH1D("h_nSigmaProton", "h_nSigmaProton", 200, 0, 10 );

    TH2D* h_2D_nSigma_Pion_Kaon   = new TH2D("h_2D_nSigma_Pion_Kaon", "h_2D_nSigma_Pion_Kaon", 20, 0, 10, 20, 0, 10 );
    TH2D* h_2D_nSigma_Pion_Proton   = new TH2D("h_2D_nSigma_Pion_Proton", "h_2D_nSigma_Pion_Proton", 20, 0, 10, 20, 0, 10 );
    TH2D* h_2D_nSigma_Kaon_Proton   = new TH2D("h_2D_nSigma_Kaon_Proton", "h_2D_nSigma_Kaon_Proton", 20, 0, 10, 20, 0, 10 );




    // 2D histos:
    TH2D* h2D_chargeSimVsRec = new TH2D("h2D_chargeSimVsRec", "h2D_chargeSimVsRec;charge rec;charge sim", 5, -2.5, 2.5, 5, -2.5, 2.5 );
    TH2D* h2D_etaSimVsRec = new TH2D("h2D_etaSimVsRec", "h2D_etaSimVsRec;#eta rec;#eta sim", 80, -2, 2, 80, -2, 2 );
    TH2D* h2D_ptSimVsRec = new TH2D("h2D_ptSimVsRec", "h2D_ptSimVsRec;pT rec;pT sim", 80, 0, 4, 80, 0, 4 );

    TH1D* h1D_RESOLUTION_etaSimVsRec = new TH1D("h1D_RESOLUTION_etaSimVsRec", "h1D_RESOLUTION_etaSimVsRec;#Delta#eta;entries", 800, -4, 4 );
    TH1D* h1D_RESOLUTION_phiSimVsRec = new TH1D("h1D_RESOLUTION_phiSimVsRec", "h1D_RESOLUTION_phiSimVsRec;#Delta#varphi;entries", 700, -7, 7 );
    TH1D* h1D_RESOLUTION_ptSimVsRec = new TH1D("h1D_RESOLUTION_ptSimVsRec", "h1D_RESOLUTION_ptSimVsRec;#Deltap_{T};entries", 800, -2, 2 );



    //  ### basic plots IA:
    const int nPIDs = 4;
    TString strPID[] = { "pions", "kaons", "protons", "electrons" };
    int arrPidIds[nPIDs] = { 211, 321, 2212, 11 };

    // energy ZDC
    TH1D *h_zdc_total_energy  = new TH1D("h_zdc_total_energy","h_zdc_total_energy", 2000, 0, 100);

    // rapidity:
    TH1D *h_y_rec_all  = new TH1D("h_y_rec_all","h_y_rec_all", _nYBins,etaMin,etaMax);
    TH1D *h_y_rec_PID[nPIDs];
    TH1D *h_y_rec_PID_plus[nPIDs];
    TH1D *h_y_rec_PID_minus[nPIDs];

    // eta:
    TH1D *h_eta_rec_all  = new TH1D("h_eta_rec_all","h_eta_rec_all", _nYBins,etaMin,etaMax);
    TH1D *h_eta_rec_PID[nPIDs];
    TH1D *h_eta_rec_PID_plus[nPIDs];
    TH1D *h_eta_rec_PID_minus[nPIDs];

    // pt:
    TH1D *h_pt_rec_all  = new TH1D("h_pt_rec_all","h_pt_rec_all", _nPtBins, ptMin, ptMax );
    TH1D *h_pt_rec_PID[nPIDs];
    TH1D *h_pt_rec_PID_plus[nPIDs];
    TH1D *h_pt_rec_PID_minus[nPIDs];

    // dEdx:
    TH2D*  h_dEdx_vs_mom_all = new TH2D("h_dEdx_vs_mom","h_dEdx_vs_mom", LogBinsX, LogWidthX,LogBinsY, LogWidthY);
    TH2D*  h_dEdx_vs_mom_PID[nPIDs];
    TH2D*  h_dEdx_vs_mom_PID_plus[nPIDs];
    TH2D*  h_dEdx_vs_mom_PID_minus[nPIDs];

    for( int iPid = 0; iPid < nPIDs; iPid++)
    {
        const char *_strPID = strPID[iPid].Data();
        h_y_rec_PID[iPid] = new TH1D( Form("h_y_%s", _strPID ), Form("h_y_%s", _strPID), _nYBins,etaMin,etaMax);
        h_y_rec_PID_plus[iPid] = new TH1D( Form("h_y_%s_plus", _strPID ), Form("h_y_%s_plus", _strPID), _nYBins,etaMin,etaMax);
        h_y_rec_PID_minus[iPid] = new TH1D( Form("h_y_%s_minus", _strPID ), Form("h_y_%s_minus", _strPID), _nYBins,etaMin,etaMax);

        h_eta_rec_PID[iPid] = new TH1D( Form("h_eta_%s", _strPID ), Form("h_eta_%s", _strPID), _nYBins,etaMin,etaMax);
        h_eta_rec_PID_plus[iPid] = new TH1D( Form("h_eta_%s_plus", _strPID ), Form("h_eta_%s_plus", _strPID), _nYBins,etaMin,etaMax);
        h_eta_rec_PID_minus[iPid] = new TH1D( Form("h_eta_%s_minus", _strPID ), Form("h_eta_%s_minus", _strPID), _nYBins,etaMin,etaMax);

        h_pt_rec_PID[iPid] = new TH1D( Form("h_pt_%s", _strPID ), Form("h_pt_%s", _strPID), _nPtBins,ptMin,ptMax);
        h_pt_rec_PID_plus[iPid] = new TH1D( Form("h_pt_%s_plus", _strPID ), Form("h_pt_%s_plus", _strPID), _nPtBins,ptMin,ptMax);
        h_pt_rec_PID_minus[iPid] = new TH1D( Form("h_pt_%s_minus", _strPID ), Form("h_pt_%s_minus", _strPID), _nPtBins,ptMin,ptMax);

        h_dEdx_vs_mom_PID[iPid] = new TH2D( Form("h_dEdx_vs_mom_PID_%s", _strPID ), Form("h_dEdx_vs_mom_PID_%s", _strPID), LogBinsX, LogWidthX,LogBinsY, LogWidthY);
        h_dEdx_vs_mom_PID_plus[iPid]   = new TH2D( Form("h_dEdx_vs_mom_PID_plus_%s", _strPID ), Form("h_dEdx_vs_mom_PID_plus_%s", _strPID), LogBinsX, LogWidthX,LogBinsY, LogWidthY);
        h_dEdx_vs_mom_PID_minus[iPid]  = new TH2D( Form("h_dEdx_vs_mom_PID_minus_%s", _strPID ), Form("h_dEdx_vs_mom_PID_minus_%s", _strPID), LogBinsX, LogWidthX,LogBinsY, LogWidthY);
    }

    // now init DCA histos:
    TString arr_DCA_hist_labels[] = { "allRec", "recPrim", "sec", };// "weak", "material",    "weakLambda", "weakK0S" };
    int nLabels = sizeof(arr_DCA_hist_labels)/sizeof(*arr_DCA_hist_labels);
    cout << "nLabels = " << nLabels << endl;
    const int nDCAzBins = 400;//1000;
    const int nPtBinsForDCA = 90;//19;//24;
    double fPtMin = 0.2;//1.2;
    double fPtMax = 2.0;//5.0;

    for ( int iType = 0; iType < nLabels; iType++ )
    {
        // all pt
        fHist1D_DCA_XY          [iType] =  new TH1D( Form( "hist1D_DCA_XY_%s", arr_DCA_hist_labels[iType].Data() ), "DCA xy;xy,cm;n tracks", 1000,0,4 );
        fHist1D_DCA_Z           [iType] =  new TH1D( Form( "hist1D_DCA_Z_%s", arr_DCA_hist_labels[iType].Data() ), "DCA z;z,cm;n tracks", nDCAzBins, -4,4 );
        fHist1D_DCA_Z_SPEC_QA   [iType] =  new TH1D( Form( "hist1D_DCA_Z_SPEC_QA_%s", arr_DCA_hist_labels[iType].Data() ), "DCA z;z,cm;n tracks", nDCAzBins, -1100, 1100 );
        fHist2D_DCA_XYvsZ       [iType] =  new TH2D( Form( "hist2D_DCA_XYvsZ_%s", arr_DCA_hist_labels[iType].Data() ), "DCA (xy vs. z);xy,cm;z,cm", 400,0,4, 400, -4,4 );

        fHist1D_DCA_XYZ                 [iType] =  new TH1D( Form( "hist1D_DCA_XYZ_%s", arr_DCA_hist_labels[iType].Data() ), "DCA r, cm;n tracks", nDCAzBins, -4,4 );
        fHist1D_DCA_XYZ_Artif_FB_Sign   [iType] =  new TH1D( Form( "hist1D_DCA_XYZ_Artif_FB_Sign_%s", arr_DCA_hist_labels[iType].Data() ), "DCA r, cm;n tracks", nDCAzBins, -4,4 );

    //        fOutputList->Add( fHist1D_DCA_XY[iType] );
    //        fOutputList->Add( fHist1D_DCA_Z[iType] );
    //        fOutputList->Add( fHist1D_DCA_Z_SPEC_QA[iType] );
    //        fOutputList->Add( fHist2D_DCA_XYvsZ[iType] );

        // vs pt:
        fHist2D_DCA_XY_vs_pt          [iType] =  new TH2D( Form( "hist2D_DCA_XY_vs_pt_%s", arr_DCA_hist_labels[iType].Data() ), "DCA xy;xy,cm;n tracks", 1000,0,4, nPtBinsForDCA, fPtMin, fPtMax );
        fHist2D_DCA_Z_vs_pt           [iType] =  new TH2D( Form( "hist2D_DCA_Z_vs_pt_%s", arr_DCA_hist_labels[iType].Data() ), "DCA z;z,cm;n tracks", nDCAzBins, -4,4, nPtBinsForDCA, fPtMin, fPtMax );
        fHist2D_DCA_Z_SPEC_QA_vs_pt   [iType] =  new TH2D( Form( "hist2D_DCA_Z_SPEC_QA_vs_pt_%s", arr_DCA_hist_labels[iType].Data() ), "DCA z;z,cm;n tracks", nDCAzBins, -1100, 1100, nPtBinsForDCA, fPtMin, fPtMax );

        fHist2D_DCA_XYZ_vs_pt           [iType] =  new TH2D( Form( "hist2D_DCA_XYZ_vs_pt_%s", arr_DCA_hist_labels[iType].Data() ), "DCA r, cm;p_{T};n tracks", nDCAzBins, -4,4, nPtBinsForDCA, fPtMin, fPtMax );
        fHist2D_DCA_XYZ_vs_pt_Artif_FB_Sign           [iType] =  new TH2D( Form( "hist2D_DCA_XYZ_vs_pt_Artif_FB_Sign_%s", arr_DCA_hist_labels[iType].Data() ), "DCA r, cm;p_{T};n tracks", nDCAzBins, -4,4, nPtBinsForDCA, fPtMin, fPtMax );
        fHist2D_DCA_XYZ_vs_pt_FB_02_06_Artif_FB_Sign  [iType] =  new TH2D( Form( "hist2D_DCA_XYZ_vs_pt_FB_02_06_Artif_FB_Sign_%s", arr_DCA_hist_labels[iType].Data() ), "DCA r, cm;p_{T};n tracks", nDCAzBins, -4,4, nPtBinsForDCA, fPtMin, fPtMax );



    //        fOutputList->Add( fHist2D_DCA_XY_vs_pt[iType] );
    //        fOutputList->Add( fHist2D_DCA_Z_vs_pt[iType] );
    // //        fOutputList->Add( fHist2D_DCA_Z_SPEC_QA_vs_pt[iType] );

        // vs pt, eta:
        fHist3D_DCA_XY_vs_PtEta          [iType] =  new TH3D( Form( "hist3D_DCA_XY_vs_PtEta_%s",  arr_DCA_hist_labels[iType].Data() ), "DCA xy;xy,cm;n tracks", 1000,0,4, nPtBinsForDCA, fPtMin, fPtMax, /*nEtaBins*/10, -1.0, 1.0 );
        fHist3D_DCA_Z_vs_PtEta           [iType] =  new TH3D( Form( "hist3D_DCA_Z_vs_PtEta_%s",   arr_DCA_hist_labels[iType].Data() ), "DCA z;z,cm;n tracks", nDCAzBins, -4,4, nPtBinsForDCA, fPtMin, fPtMax, /*nEtaBins*/10, -1.0, 1.0 );
        fHist3D_DCA_XYZ_vs_PtEta         [iType] =  new TH3D( Form( "hist3D_DCA_XYZ_vs_PtEta_%s", arr_DCA_hist_labels[iType].Data() ), "DCA r;r,cm;n tracks", 1000,0,4, nPtBinsForDCA, fPtMin, fPtMax, /*nEtaBins*/10, -1.0, 1.0 );
        fHist3D_DCA_XYZ_vs_PtEta_Artif_FB_Sign         [iType] =  new TH3D( Form( "hist3D_DCA_XYZ_vs_PtEta_Artif_FB_Sign_%s", arr_DCA_hist_labels[iType].Data() ), "DCA r;r,cm;n tracks", 2000,-4,4, nPtBinsForDCA, fPtMin, fPtMax, /*nEtaBins*/10, -1.0, 1.0 );

    //        fOutputList->Add( fHist3D_DCA_XY_vs_PtEta[iType] );
    //        fOutputList->Add( fHist3D_DCA_Z_vs_PtEta[iType] );

    } // end of particle type loop
    // ### end of DCA histos:


    //  ### end of IA histos



//    return;




    // ##### EVENT TREE:
    static const Int_t NMaxTrack = 12000;

    Int_t fEventId; //!
    Float_t fImpactPar;   //!
    Float_t fVertexX,fVertexY,fVertexZ; //!
    Float_t f_zdc_total_energy; //!

    fEventTree->SetBranchAddress("fEventId",&fEventId );
    fEventTree->SetBranchAddress("fImpactPar",&fImpactPar );

    fEventTree->SetBranchAddress("vertexX",&fVertexX );
    fEventTree->SetBranchAddress("vertexY",&fVertexY );
    fEventTree->SetBranchAddress("vertexZ",&fVertexZ );

    fEventTree->SetBranchAddress("zdc_total_energy", &f_zdc_total_energy );

    // SIM:
    Int_t f_SIM_NumberOfTracks; //!
    Int_t f_SIM_TrackId[NMaxTrack];
    Float_t f_SIM_TrackPt[NMaxTrack],f_SIM_TrackEta[NMaxTrack],f_SIM_TrackPhi[NMaxTrack], f_SIM_TrackY[NMaxTrack]; //
    Int_t f_SIM_TrackCharge[NMaxTrack];
    Int_t f_SIM_pid[NMaxTrack];
    Int_t f_SIM_MotherId[NMaxTrack];

    fEventTree->SetBranchAddress("SIM_NumberOfTracks",&f_SIM_NumberOfTracks );
    fEventTree->SetBranchAddress("SIM_TrackId",f_SIM_TrackId );

    fEventTree->SetBranchAddress("SIM_TrackPt",f_SIM_TrackPt   );
    fEventTree->SetBranchAddress("SIM_TrackPhi",f_SIM_TrackPhi );
    fEventTree->SetBranchAddress("SIM_TrackEta",f_SIM_TrackEta );
    fEventTree->SetBranchAddress("SIM_TrackY",f_SIM_TrackY );
    fEventTree->SetBranchAddress("SIM_TrackCharge",f_SIM_TrackCharge );
    fEventTree->SetBranchAddress("SIM_pid",f_SIM_pid );
    fEventTree->SetBranchAddress("SIM_MotherId",f_SIM_MotherId );

    // REC:
    Int_t fNumberOfTracks; //!
    Int_t fTrackSimID[NMaxTrack];
    Int_t fTrackSimMotherId[NMaxTrack];

    Float_t fTrackPt[NMaxTrack],fTrackEta[NMaxTrack],fTrackPhi[NMaxTrack]; //
    Float_t fTrackPtBeforeAbs[NMaxTrack];
    Int_t fTrackCharge[NMaxTrack];
    Int_t fTrackChargeUsingPt[NMaxTrack];
    Float_t fTrackMom[NMaxTrack];
    Float_t fdEdx[NMaxTrack];
    Float_t fDCA_x[NMaxTrack];
    Float_t fDCA_y[NMaxTrack];
    Float_t fDCA_z[NMaxTrack];
    Int_t fnTPChits[NMaxTrack];
//    Int_t fnKalmanHits[NMaxTrack];
    Float_t fTrackChi2[NMaxTrack];

    fEventTree->SetBranchAddress("nTracks",&fNumberOfTracks );
    fEventTree->SetBranchAddress("trackSimID",fTrackSimID );
    fEventTree->SetBranchAddress("trackSimMotherId",fTrackSimMotherId );

    fEventTree->SetBranchAddress("trackPt",fTrackPt );
    fEventTree->SetBranchAddress("trackPhi",fTrackPhi );
    fEventTree->SetBranchAddress("trackEta",fTrackEta );
    fEventTree->SetBranchAddress("trackPtBeforeAbs",fTrackPtBeforeAbs );

    fEventTree->SetBranchAddress("trackCharge",fTrackCharge );
//    fEventTree->SetBranchAddress("trackChargeUsingPt",fTrackChargeUsingPt );
    fEventTree->SetBranchAddress("trackMom",fTrackMom );
    fEventTree->SetBranchAddress("dEdx",fdEdx );

    fEventTree->SetBranchAddress("DCA_x",fDCA_x );
    fEventTree->SetBranchAddress("DCA_y",fDCA_y );
    fEventTree->SetBranchAddress("DCA_z",fDCA_z );

    fEventTree->SetBranchAddress("nTPChits",fnTPChits );
//    fEventTree->SetBranchAddress("fnKalmanHits",fnKalmanHits );
    fEventTree->SetBranchAddress( "trackChi2",fTrackChi2 );




//return;




    // ##### TREE WITH TRACKS FOR IDENTITY METHOD:
    ULong64_t IM_event; //!
    Float_t IM_vZ; //!
    Float_t IM_dca; //!
    Float_t IM_dcaR; //!
    Int_t IM_sign; //!
    Float_t IM_centrality; //!
    Float_t IM_pt; //!
    Float_t IM_mom; //!
    Float_t IM_eta; //!
    Float_t IM_dEdx; //!
    Int_t IM_prim_or_sec; //!
    Int_t IM_sim_PID; //!

    TTree       *fTreeIM = new TTree("treeIM","treeIM"); //! output Tree
    fTreeIM->Branch("event",        &IM_event,      "event/l");
//    fTreeIM->Branch("vZ",           &IM_vZ,        "vZ/F");
//    fTreeIM->Branch("dca",          &IM_dca,        "dca/F");
//    fTreeIM->Branch("dcaR",         &IM_dcaR,        "dcaR/F");
//    fTreeIM->Branch("sign",         &IM_sign,       "sign/I");
    fTreeIM->Branch("centrality",   &IM_centrality,    "centrality/F");
    fTreeIM->Branch("pt",   &IM_pt,    "pt/F");
    fTreeIM->Branch("mom",   &IM_mom,    "mom/F");
    fTreeIM->Branch("dEdx",   &IM_dEdx,    "dEdx/F");
    fTreeIM->Branch("eta",   &IM_eta,    "eta/F");

//    fTreeIM->Branch("prim_or_sec",   &IM_prim_or_sec,    "eta/I");
    fTreeIM->Branch("sim_PID",   &IM_sim_PID,    "sim_PID/I");






    // ### event loop
    int counterTracksAllRec = 0;
    int counterTracksRecSimChargeMatch = 0;
    int counterNoSimFound = 0;

    double meanMultInRun = 0;
    int nEvInRun = 0;
    int runCounter = 0;

    Int_t nEvents = fEventTree->GetEntries();
//    Int_t nEvents = 2000;
    cout << " Number of events in tree = " << nEvents << endl;
    cout << endl;

//    bool firstRunFirstEvent = true;
    int prevEventId = -1;


    int QA_HaveSimPID = 0;
    int QA_NoSimPID = 0;
    int QA_mom_outside= 0;

    // ### EVENT LOOP
    for ( int iEv = 0; iEv < nEvents; iEv++ )
        //        for ( int iEv = 0; iEv < 2000; iEv++ )
    {
        if( iEv%500 == 0 )
            cout << "analysing " << (int)iEv << "/" << nEvents << " event... \r"; cout.flush();
        fEventTree->GetEntry(iEv);


        // raw primary vertex:
        h_VertexXY->Fill( fVertexX, fVertexY );
        h_VertexZ->Fill( fVertexZ );


        h_ImpactPar->Fill( fImpactPar );


        // !!! trigger-like cut:
        if ( fNumberOfTracks < 2 )
            continue;

        // event-level vZ cut:
        if ( fabs(fVertexZ) > 50 )
            continue;


        // !!! ELIMINATE zero vZ for both SIM and REC!
        if ( fabs(fVertexZ) < 0.01 ) // = no REC vertex (?!)
            continue;








        // ### pre-loop REC tracks for quality plots:
        int nRecInCuts = 0;
        for ( Int_t iP = 0; iP < fNumberOfTracks; iP++ )
        {
            double eta = fTrackEta[iP];
            //            double AbsEta = TMath::Abs(eta);
            double phi = fTrackPhi[iP];
            if ( phi < 0 )
                phi += TMath::TwoPi();



//            double mom = fTrackMom[iP];
            double pt = fTrackPt[iP];
            int charge = fTrackCharge[iP];
            pt = TMath::Abs(pt);
//            double dedx = fdEdx[iP];


            // ### IA:
            h_Pt->Fill(pt);                     // plot the pt value of the track in a histogram
            h_Eta->Fill(eta);
            h_Phi->Fill(phi);
            h_Charge->Fill(charge);
            //            h_Y->Fill(trackMC->Y());

            // FOR EFFICIENCY PLOTS:
            int motherId = fTrackSimMotherId[iP];
            if( fabs(eta)<1.0 )//1.6 )
            {
                h_Pt_in_eta1->Fill( pt );
                if( motherId == -1 ) // primary
                    h_RecPrim_Pt_in_eta1->Fill( eta );
            }
            if( pt > 0.15 && pt < 2.0 )
            {
                h_Eta_in_pt015_2->Fill( eta );
                if( motherId == -1 ) // primary
                    h_RecPrim_Eta_in_pt015_2->Fill( eta );
            }




            //                h_y_rec_all->Fill( rapidity_MC );
            h_eta_rec_all->Fill( eta );
            h_pt_rec_all->Fill( pt );

            int nTPChits = fnTPChits[iP];

            if( fabs(fVertexZ)>0.001 )
            {
                h2D_eta_vs_vZ->Fill( fVertexZ, eta );
                if ( nTPChits > 30 )
                    h2D_eta_vs_vZ_if_nClustAbove30->Fill( fVertexZ, eta );
            }


            // number of hits VS vZ
            const double kCutDCAz = 1.0;
            double dcaZ = fDCA_z[iP];
            h2D_hits_vs_vZ->Fill( fVertexZ, nTPChits );
            if( fabs(eta)<1.0 )
            {
                h2D_hits_vs_vZ_if_etaWithin1->Fill( fVertexZ, nTPChits );
                if(nTPChits>0)
                {
                    h2D_chi2_perCluster_vs_vZ->Fill( fVertexZ, fTrackChi2[iP]/nTPChits );
                    h1D_chi2_perCluster->Fill( fTrackChi2[iP]/nTPChits );
                    h1D_nClusters->Fill( nTPChits );
                }
            }


            if( fVertexZ > 20 )
            {
                h2D_hits_vs_eta_if_vZ_above_20->Fill( eta, nTPChits );
                if( fabs(dcaZ)<kCutDCAz ) h2D_STRICT_DCAz_hits_vs_eta_if_vZ_above_20->Fill( eta, nTPChits );
            }
            else if( fVertexZ < -20 )
            {
                h2D_hits_vs_eta_if_vZ_below_minus20->Fill( eta, nTPChits );
                if( fabs(dcaZ)<kCutDCAz ) h2D_STRICT_DCAz_hits_vs_eta_if_vZ_below_minus20->Fill( eta, nTPChits );
            }
            else
            {
                h2D_hits_vs_eta_if_vZ_WITHIN_20->Fill( eta, nTPChits );
                if( fabs(dcaZ)<kCutDCAz ) h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_20->Fill( eta, nTPChits );
            }

            if( fabs(fVertexZ) < 10 )
            {
                h2D_hits_vs_eta_if_vZ_WITHIN_10->Fill( eta, nTPChits );
                if( fabs(dcaZ)<kCutDCAz ) h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_10->Fill( eta, nTPChits );
            }
            if( fabs(fVertexZ) < 2 && fabs(fVertexZ) > 0.01 )
            {
                h2D_hits_vs_eta_if_vZ_WITHIN_2_excl_0->Fill( eta, nTPChits );
                if( fabs(dcaZ)<kCutDCAz ) h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_2_excl_0->Fill( eta, nTPChits );
            }



            if( fabs(eta)<1.0/*1.6*/ && (pt > 0.15 && pt < 2.0) )
            {
                nRecInCuts++;
            }

        } // end of pre-loop REC tracks for quality plots:


        h_nTracksInEtaCuts_REC->Fill( nRecInCuts );




        // ### !!! vZ cut for analysis
        if ( fabs(fVertexZ) > 20 )
//            if ( fVertexZ < 1 || fVertexZ > 5)//20 )
            continue;
        h_VertexZ_accepted->Fill( fVertexZ );
//        h_VertexZ_accepted_excl_0->Fill( fVertexZ );



        // ### ZDC info
        h_zdc_total_energy->Fill( f_zdc_total_energy );



        // ### determine centrality by REC mult in TPC:
        int nRecForMultCentrality = 0;
        for ( Int_t iP = 0; iP < fNumberOfTracks; iP++ )
            if( fabs( fTrackEta[iP] )<1.0 && ( fTrackPt[iP] > 0.15) && ( fnTPChits[iP] > 30 ) )
            {
                nRecForMultCentrality++;
            }
        h_nTracksForMultCentrality->Fill( nRecForMultCentrality );

        // if zero number of tracks - continue!
//        if ( nRecForMultCentrality==0 )
//            continue;



        // determine centrality bin!
        int cBin = -1;
        for( int j = 0; j < nCbins; j++ )
            //            if ( bImp < bImpRanges[j] )
            //            if ( f_zdc_total_energy < energyRanges[j] )
            if ( nRecForMultCentrality < multRanges[j] && nRecForMultCentrality > multRanges[j+1] )
            {
                cBin = j;
                break;
            }
        if ( cBin >= nCbins ) // just outside of max b_imp thresh.
            continue;
        if ( cBin == -1 )
        {
            //            cout << "AHTUNG!!! cBin == -1" << endl;
            continue;
        }
        int cBinId[10];
        cBinId[0] = cBin; // tmp

        //        cout << "cBin : " << cBinId[0] << endl;





        h_QA_centr_bins->Fill( cBin );


        if ( nRecInCuts > 25 )
        {
            //            h1D_vertexZ->Fill(vZ);
            h1D_vertexZspec->Fill(fVertexZ);
        }










        // assign subsample id
        int subsampleId = gRandom->Integer(nWinPairWithSubsamples);


        // ### SIM particle loop
        int nSimInCuts = 0;
        for ( Int_t iP = 0; iP < f_SIM_NumberOfTracks; iP++ )
        {
            int pdg = f_SIM_pid[iP];

            int charge = f_SIM_TrackCharge[iP];

            // consider only primary for SIM level!
            if ( f_SIM_MotherId[iP] != -1 )
                continue;

            if ( charge==0 )
                continue;

            double eta = f_SIM_TrackEta[iP];
            double pt = f_SIM_TrackPt[iP];
            double phi = f_SIM_TrackPhi[iP];
            double rapidity = f_SIM_TrackY[iP];

            h_SIM_Pt->Fill( pt );
            h_SIM_Eta->Fill( eta );
            h_SIM_Rapidity->Fill( rapidity );
            h_SIM_Phi->Fill( phi );

            // FOR EFFICIENCY PLOTS:
            if( fabs(eta)<1.0 )//1.6 )
                h_SIM_Pt_in_eta1->Fill( pt );
            if( pt > 0.15 && pt < 2 )
                h_SIM_Eta_in_pt015_2->Fill( eta );

            if( fabs(eta)<1.0 && (pt > 0.15 && pt < 10) )
                nSimInCuts++;


            // add track to winPairs_SIM
            for ( int iType = 0; iType < nPartTypes; iType++)
                for ( int iCW = 0; iCW < nCW; iCW++)
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                        for ( int iEta = 0; iEta < nEtaBins; ++iEta )
                            winPairs_SIM[subsampleId][iType][iCW][ cBinId[iCW] ][iPt][iEta].addTrack( pdg, eta, phi, pt, charge );

        } // end of loop over SIM particles
        h_nTracksInEtaCuts_SIM->Fill( nSimInCuts );


//        cout << "iEv = " << iEv << ", fEventId = " << fEventId << endl;
        // QA: write av mult in this file (if 1000 events per file!)
//        if ( (iEv+1)%1000 == 0 ) // new file will start next event!
        if ( fEventId < prevEventId ) // && !firstRunFirstEvent )
        {
            cout << endl;
            cout << "QA: iEv = " << iEv << ", fEventId = " << fEventId << endl;
            if(nEvInRun>0)
            {
                h_QA_meanMultInFiles->SetBinContent( runCounter+1, meanMultInRun/nEvInRun );
                cout << "meanMultInRun/nEvInRun = " << meanMultInRun/nEvInRun << endl;
            }
            else
                cout << "!!! nEvInRun=0!" << endl;

            meanMultInRun = 0;
            nEvInRun = 0;

            runCounter++;
        }
//        firstRunFirstEvent = false;
        prevEventId = fEventId;


        meanMultInRun += nSimInCuts;
        nEvInRun++;




        // ### main REC track loop for analysis
        int nAccTracks = 0;
        for ( Int_t iP = 0; iP < fNumberOfTracks; iP++ )
        {
            //            mpdTrack = (MpdTrack*) mpdTracks->UncheckedAt(trackId);

            //            int ID = mpdTrack->GetID();
            //            if ( flag_array[ID] == 1 )
            //                continue;
            //            flag_array[ID] = 1;

            //            mcTrack = (FairMCTrack*) mcTracks->UncheckedAt(ID);
            //            int pdgc = mcTrack->GetPdgCode();
            //            int mother = mcTrack->GetMotherId();
            //            double rapidity_MC = mcTrack->GetRapidity();

            // cut on n TPC hits:
            int nTPChits = fnTPChits[iP];
            if ( nTPChits < 30 )
                continue;

            double eta = fTrackEta[iP];
            //            double AbsEta = TMath::Abs(eta);
            double phi = fTrackPhi[iP];
            if ( phi < 0 )
                phi += TMath::TwoPi();



            double mom = fTrackMom[iP];
            double pt = fTrackPt[iP];
            int charge = fTrackCharge[iP];
            pt = TMath::Abs(pt);
            double dedx = fdEdx[iP];
            int motherId = fTrackSimMotherId[iP];


            // ### IA:
            h_ACCEPTED_Pt->Fill(pt);                     // plot the pt value of the track in a histogram
            h_ACCEPTED_Eta->Fill(eta);
            h_ACCEPTED_Phi->Fill(phi);
            h_ACCEPTED_Charge->Fill(charge);
            //            h_Y->Fill(trackMC->Y());



            double dcaXY = sqrt( fDCA_x[iP]*fDCA_x[iP] + fDCA_y[iP]*fDCA_y[iP] );
            double dcaZ = fDCA_z[iP];


            //cout << "test1"<< endl;


            // chi2perCluster cut:
            double chi2perCluster = fTrackChi2[iP]/nTPChits;
            h1D_BEFORE_CUT_chi2_perCluster->Fill( chi2perCluster );
            if( chi2perCluster > 5)
                continue;


            // ### Fill DCA
            if( fabs(eta) < 1.0 && cBinId[0] == 0 ) // take only most central events!
            {
                // all rec
                FillHistosDCA( 0, dcaXY, dcaZ, pt, eta );
                if ( motherId == -1 ) // primary
                    FillHistosDCA( 1, dcaXY, dcaZ, pt, eta );
                else // secondary
                    FillHistosDCA( 2, dcaXY, dcaZ, pt, eta );
            }


            // #### TREE for Identity Method:
//            if( fabs(eta) < 1.0 && cBinId[0] == 0) //1.6 )
            if( fabs(eta) < 1.0 && ( cBin == 0 || cBin == 1 ) ) //1.6 )
            {
                // find id of corresponding sim track in sim array
                int simPID = -1000;
                if(1)
                {
                    int simId = fTrackSimID[iP];

                    int idInSimArr = -1;
                    for ( Int_t iSim = 0; iSim < f_SIM_NumberOfTracks; iSim++ )
                        if ( simId == f_SIM_TrackId[iSim] )
                        {
                            idInSimArr = iSim;
                            break;
                        }
                    if ( idInSimArr != -1 )
                        simPID = f_SIM_pid[idInSimArr];

                } // end of matching sim-rec


//                bool isB = ( eta > -1.0 && eta < -0.4 ) ? true : false;
//                bool isF = ( eta > 0.4  && eta < 1.0  ) ? true : false;
                bool isB = ( eta < 0 ) ? true : false;
//                bool isF = ( eta > 0 ) ? true : false;
//                if( isB || isF )  // fixed FB wins! //1.6 )
//                {
                if( simPID == -1000 )
                {
                    QA_NoSimPID++;//cout << "AHTUNG!" << endl;
                }
                else
                {
                    QA_HaveSimPID++;
                    if( mom>0 && mom<10 && fabs(dcaZ)<2 && fabs(dcaXY)<2 )
                    {
                        IM_event = fEventId;
                        IM_vZ = fVertexZ;
                        IM_dca = dcaZ;
                        IM_dcaR = sqrt( dcaZ*dcaZ + dcaXY*dcaXY );
                        IM_sign = isB ? -1 : +1; // artificial "sign" for IM!  //charge;
                        IM_centrality = cBin;
                        IM_pt = pt;
                        IM_mom = mom; //charge > 0 ? mom : -mom;
                        IM_dEdx = charge > 0 ? dedx : -dedx;
                        IM_eta = eta;

                        IM_sim_PID = simPID;
                        IM_prim_or_sec = ( motherId == -1 ) ? 0 : 1; // primary or secondary: 0-prim, 1-sec


                        fTreeIM->Fill();

                    }
                    else
                        QA_mom_outside++;

                }
            }



            // !!! DCAz cut!
            if( fabs(dcaZ) > 2 )
                continue;

            // !!! DCAxy cut!
            if( fabs(dcaXY) > 1 )
                continue;



            // ### FOR EFFICIENCY PLOTS:
            if( fabs(eta) < 1.0 )//1.6 )
            {
                h_ACCEPTED_Pt_in_eta1->Fill( pt );
                if( motherId == -1 ) // primary
                    h_ACCEPTED_RecPrim_Pt_in_eta1->Fill( pt );
            }
            if( pt > 0.15 && pt < 2 )
            {
                h_ACCEPTED_Eta_in_pt015_2->Fill( eta );
                if( motherId == -1 ) // primary
                    h_ACCEPTED_RecPrim_Eta_in_pt015_2->Fill( eta );
            }




            // ... continue kinematic cuts:
            if( fabs(eta) > 1.0 )
                continue;

            if( pt < 0.15 || pt > 2.0 )
                continue;

            nAccTracks++;
            counterTracksAllRec++;

            h1D_ACCEPTED_chi2_perCluster->Fill( chi2perCluster );
            h1D_ACCEPTED_nClusters->Fill( nTPChits );


            // add track to winPairs_REC
            for ( int iType = 0; iType < nPartTypes; iType++)
                for ( int iCW = 0; iCW < nCW; iCW++)
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                        for ( int iEta = 0; iEta < nEtaBins; ++iEta )
                        {
                            winPairs_REC[subsampleId][iType][iCW][ cBinId[iCW] ][iPt][iEta].addTrack( 0/*pid*/, eta, phi, pt, charge );
                            if ( motherId == -1 ) // prim
                                winPairs_REC_PRIM[subsampleId][iType][iCW][ cBinId[iCW] ][iPt][iEta].addTrack( 0/*pid*/, eta, phi, pt, charge );
                        }





//            cout << "test2"<< endl;

            //            bool ret = kFALSE;

            if (1) //(( pidFlag == 0) || (pidFlag == 2))
            {
                // FOR nSigma
                //                ret = pid->FillProbs(mom, dedx, charge); /// only dE/dx available
                //                if(ret==kTRUE)
                //                {
                //                    double nSigmaPion   = pid->GetNsigmaToBetheBloch("pi");
                //                    double nSigmaKaon   = pid->GetNsigmaToBetheBloch("ka");
                //                    double nSigmaProton = pid->GetNsigmaToBetheBloch("pr");
                //        //            cout << "n sigma: " << nSigmaPion << " " << nSigmaKaon << " " << nSigmaProton << endl;


                //                    h_nSigmaPion   ->Fill( nSigmaPion   );
                //                    h_nSigmaKaon   ->Fill( nSigmaKaon   );
                //                    h_nSigmaProton ->Fill( nSigmaProton );

                //                    h_2D_nSigma_Pion_Kaon   ->Fill( nSigmaPion,  nSigmaKaon  );
                //                    h_2D_nSigma_Pion_Proton   ->Fill( nSigmaPion,  nSigmaProton  );
                //                    h_2D_nSigma_Kaon_Proton   ->Fill( nSigmaKaon,  nSigmaProton  );
                //                }



                // other fillings:
                h_dEdx->Fill(dedx);
                //                h_dEdx_vs_mom->Fill( mom/charge, dedx );
                //                cout << "dedx = " << dedx << endl;




                double mom_to_charge = fabs(mom/charge);

                h_dEdx_vs_mom_all->Fill( mom_to_charge, dedx );

                // loop over PIDs:
                //                for( int iPid = 0; iPid < nPIDs; iPid++)
                //                    if ( fabs(pdgc) == arrPidIds[iPid] )
                //                    {
                //                        // rapidity
                //                        h_y_rec_PID[iPid]->Fill( rapidity_MC );
                //                        if (charge>0)       h_y_rec_PID_plus[iPid]->Fill( rapidity_MC );
                //                        else if (charge<0)  h_y_rec_PID_minus[iPid]->Fill( rapidity_MC );

                //                        // eta
                //                        h_eta_rec_PID[iPid]->Fill( eta );
                //                        if (charge>0)       h_eta_rec_PID_plus[iPid]->Fill( eta );
                //                        else if (charge<0)  h_eta_rec_PID_minus[iPid]->Fill( eta );

                //                        // pt
                //                        h_pt_rec_PID[iPid]->Fill( pt );
                //                        if (charge>0)       h_pt_rec_PID_plus[iPid]->Fill( pt );
                //                        else if (charge<0)  h_pt_rec_PID_minus[iPid]->Fill( pt );

                //                        // dedx:
                //                        h_dEdx_vs_mom_PID[iPid]->Fill( mom_to_charge, dedx );
                //                        if (charge>0)       h_dEdx_vs_mom_PID_plus[iPid]->Fill( mom_to_charge, dedx );
                //                        else if (charge<0)  h_dEdx_vs_mom_PID_minus[iPid]->Fill( mom_to_charge, dedx );

                //                        break;
                //                    }
            } // end of pidFlag == 0

//            cout << "test3"<< endl;





            // find id of corresponding sim track in sim array
            if(1)
            {
                int simId = fTrackSimID[iP];

                int idInSimArr = -1;
                for ( Int_t iP = 0; iP < f_SIM_NumberOfTracks; iP++ )
                    if ( simId == f_SIM_TrackId[iP] )
                    {
                        idInSimArr = iP;
                        break;
                    }

                if ( idInSimArr != -1 )
                {
                    int chargeSim = f_SIM_TrackCharge[idInSimArr];
                    double etaSim = f_SIM_TrackEta[idInSimArr];
                    double phiSim = f_SIM_TrackPhi[idInSimArr];
                    double ptSim = f_SIM_TrackPt[idInSimArr];

                    //
                    if ( charge == chargeSim )
                        counterTracksRecSimChargeMatch++;

                    h2D_chargeSimVsRec->Fill( charge, chargeSim );
                    h2D_etaSimVsRec->Fill( eta, etaSim );
                    h2D_ptSimVsRec->Fill( pt, ptSim );

                    h1D_RESOLUTION_etaSimVsRec->Fill( eta-etaSim );
                    h1D_RESOLUTION_phiSimVsRec->Fill( phi-phiSim );
                    if(ptSim>0.01) h1D_RESOLUTION_ptSimVsRec->Fill( (pt-ptSim)/ptSim );
                }
                else
                    counterNoSimFound++;
            } // end of matching sim-rec





        } // end of REC track

        h_nTracksInEtaCuts_REC_ACCEPTED->Fill( nAccTracks );
        //        delete [] flag_array;



//        cout << "test_after"<< endl;

        // finish event for winPairs SIM & REC:
        for ( int iType = 0; iType < nPartTypes; iType++)
            for ( int iCW = 0; iCW < nCW; iCW++)
                for ( int iPt = 0; iPt < nPtBins; ++iPt )
                    for ( int iEta = 0; iEta < nEtaBins; ++iEta )
                    {
                        winPairs_SIM[subsampleId][iType][iCW][ cBinId[iCW] ][iPt][iEta].finishEvent();
                        winPairs_REC[subsampleId][iType][iCW][ cBinId[iCW] ][iPt][iEta].finishEvent();
                        winPairs_REC_PRIM[subsampleId][iType][iCW][ cBinId[iCW] ][iPt][iEta].finishEvent();

                    }




//        cout << "test_after2"<< endl;

//        continue;



//        // ### now match REC tracks to SIM tracks
//        for (Int_t trackId = 0; trackId < fNumberOfTracks; trackId++)
//        {
//            int simId = fTrackSimID[trackId];
//            int chargeRec = fTrackCharge[trackId];
//            double etaRec = fTrackEta[trackId];
//            double ptRec = fTrackPt[trackId];

//            if ( fabs(etaRec) > 1.0 ) //1.4 )
//                continue;
//            if ( ptRec < 0.15 )//0.1 )
//                continue;



//            // find id of corresponding sim track in sim array
//            int idInSimArr = -1;
//            if(1)
//                for ( Int_t iP = 0; iP < f_SIM_NumberOfTracks; iP++ )
//                    if ( simId == f_SIM_TrackId[iP] )
//                    {
//                        idInSimArr = iP;
//                        break;
//                    }

//            int chargeSim = f_SIM_TrackCharge[idInSimArr];
//            double etaSim = f_SIM_TrackEta[idInSimArr];
//            double ptSim = f_SIM_TrackPt[idInSimArr];


//            //
//            nAccTracks++;
//            counterTracksAllRec++;
//            if ( chargeRec == chargeSim )
//                counterTracksRecSimChargeMatch++;

//            h2D_chargeSimVsRec->Fill( chargeSim, chargeRec);
//            h2D_etaSimVsRec->Fill( etaSim, etaRec );
//            h2D_ptSimVsRec->Fill( ptSim, ptRec );


//        }




    } // ### end of event loop


    cout << endl;
    cout << endl;
    cout << endl;
    cout << "QA_HaveSimPID = " << QA_HaveSimPID << endl;
    cout << "QA_NoSimPID = " << QA_NoSimPID << endl;
    cout << "QA_mom_outside = " << QA_mom_outside << endl;

    // aftermath:
    h_QA_meanMultInFiles->SetBinContent( runCounter+1, meanMultInRun/nEvInRun );


//return;

    // ### write output
    TFile out(outname,"recreate");

    // ##### IA:
    h_ImpactPar->Write();

    h_nTracksInEtaCuts_SIM->Write();
    h_nTracksInEtaCuts_REC->Write();
    h_nTracksInEtaCuts_REC_ACCEPTED->Write();
    h_nTracksForMultCentrality->Write();
    h_QA_centr_bins->Write();


    h_QA_meanMultInFiles->Write();

    h_VertexXY->Write();
    h_VertexZ->Write();
    h_VertexZ_accepted->Write();
//    h_VertexZ_accepted_excl_0->Write();

    h1D_vertexZspec->Write();
    h2D_eta_vs_vZ->Write();
    h2D_eta_vs_vZ_if_nClustAbove30->Write();

    h2D_chi2_perCluster_vs_vZ->Write();

    h2D_hits_vs_vZ->Write();
    h2D_hits_vs_vZ_if_etaWithin1->Write();
    h2D_hits_vs_eta_if_vZ_above_20->Write();
    h2D_hits_vs_eta_if_vZ_below_minus20->Write();
    h2D_hits_vs_eta_if_vZ_WITHIN_20->Write();
    h2D_hits_vs_eta_if_vZ_WITHIN_10->Write();
    h2D_hits_vs_eta_if_vZ_WITHIN_2_excl_0->Write();

    // if STRICT DCAz
    h2D_STRICT_DCAz_hits_vs_eta_if_vZ_above_20->Write();
    h2D_STRICT_DCAz_hits_vs_eta_if_vZ_below_minus20->Write();
    h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_20->Write();
    h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_10->Write();
    h2D_STRICT_DCAz_hits_vs_eta_if_vZ_WITHIN_2_excl_0->Write();




    h2D_chargeSimVsRec->Write();
    h2D_etaSimVsRec->Write();
    h2D_ptSimVsRec->Write();

    h1D_RESOLUTION_etaSimVsRec->Write();
    h1D_RESOLUTION_phiSimVsRec->Write();
    h1D_RESOLUTION_ptSimVsRec ->Write();

    // ### SIM:
    h_SIM_Pt->Write();
    h_SIM_Eta->Write();
    h_SIM_Rapidity->Write();
    h_SIM_Phi->Write();

    h_SIM_Pt_in_eta1->Write();
    h_SIM_Eta_in_pt015_2->Write();


    // ### REC raw:
    h_Pt->Write();
    h_Eta->Write();
    h_Phi->Write();
    h_Charge->Write();

    h_Pt_in_eta1->Write();
    h_Eta_in_pt015_2->Write();

    h_RecPrim_Pt_in_eta1->Write();
    h_RecPrim_Eta_in_pt015_2->Write();

    h1D_chi2_perCluster->Write();
    h1D_nClusters->Write();





    // ### REC within vZ range:
    h_ACCEPTED_Pt->Write();
    h_ACCEPTED_Eta->Write();
    h_ACCEPTED_Phi->Write();
    h_ACCEPTED_Charge->Write();

    h_ACCEPTED_Pt_in_eta1->Write();
    h_ACCEPTED_Eta_in_pt015_2->Write();

    h_ACCEPTED_RecPrim_Pt_in_eta1->Write();
    h_ACCEPTED_RecPrim_Eta_in_pt015_2->Write();

    h1D_BEFORE_CUT_chi2_perCluster->Write();
    h1D_ACCEPTED_chi2_perCluster->Write();
    h1D_ACCEPTED_nClusters->Write();




    // efficiencies (all rec):
    TH1D *clone_h_ACCEPTED_Pt_in_eta1 = (TH1D*)h_ACCEPTED_Pt_in_eta1->Clone("clone_h_ACCEPTED_Pt_in_eta1");
    clone_h_ACCEPTED_Pt_in_eta1->Divide( h_SIM_Pt_in_eta1 );
    clone_h_ACCEPTED_Pt_in_eta1->Write( "h_eff_vs_pt_in_eta1");

    TH1D *clone_h_ACCEPTED_Eta_in_pt015_2 = (TH1D*)h_ACCEPTED_Eta_in_pt015_2->Clone("clone_h_ACCEPTED_Eta_in_pt015_2");
    clone_h_ACCEPTED_Eta_in_pt015_2->Divide( h_SIM_Eta_in_pt015_2 );
    clone_h_ACCEPTED_Eta_in_pt015_2->Write( "h_eff_vs_eta_in_pt015_2");

    // efficiencies (rec prim):
    TH1D *clone_h_ACCEPTED_RecPrim_Pt_in_eta1 = (TH1D*)h_ACCEPTED_RecPrim_Pt_in_eta1->Clone("clone_h_ACCEPTED_RecPrim_Pt_in_eta1");
    clone_h_ACCEPTED_RecPrim_Pt_in_eta1->Divide( h_SIM_Pt_in_eta1 );
    clone_h_ACCEPTED_RecPrim_Pt_in_eta1->Write( "h_eff_PRIMARY_vs_pt_in_eta1");

    TH1D *clone_h_ACCEPTED_RecPrim_Eta_in_pt015_2 = (TH1D*)h_ACCEPTED_RecPrim_Eta_in_pt015_2->Clone("clone_h_ACCEPTED_RecPrim_Eta_in_pt015_2");
    clone_h_ACCEPTED_RecPrim_Eta_in_pt015_2->Divide( h_SIM_Eta_in_pt015_2 );
    clone_h_ACCEPTED_RecPrim_Eta_in_pt015_2->Write( "h_eff_PRIMARY_vs_eta_in_pt015_2");


    // dE/dx
    h_dEdx->Write();
    //    h_dEdx_vs_mom->Write();
    //    h_dEdx_vs_mom_pions_plus->Write();

    //    h_dEdx_vs_mom_pions_plus ->Write();
    //    h_dEdx_vs_mom_kaons_plus ->Write();
    //    h_dEdx_vs_mom_protons_plus->Write();
    //    h_dEdx_vs_mom_electrons_plus ->Write();

    //    h_dEdx_vs_mom_pions_minus->Write();
    //    h_dEdx_vs_mom_kaons_minus->Write();
    //    h_dEdx_vs_mom_protons_minus ->Write();
    //    h_dEdx_vs_mom_electrons_minus ->Write();




    // rapidity:
    h_y_rec_all->Write();
    h_eta_rec_all->Write();
    h_pt_rec_all->Write();
    h_dEdx_vs_mom_all->Write();

    if(0)for( int iPid = 0; iPid < nPIDs; iPid++)
    {
        h_y_rec_PID[iPid]->Write();
        h_y_rec_PID_plus[iPid]->Write();
        h_y_rec_PID_minus[iPid]->Write();

        h_eta_rec_PID[iPid]->Write();
        h_eta_rec_PID_plus[iPid]->Write();
        h_eta_rec_PID_minus[iPid]->Write();

        h_pt_rec_PID[iPid]->Write();
        h_pt_rec_PID_plus[iPid]->Write();
        h_pt_rec_PID_minus[iPid]->Write();

        h_dEdx_vs_mom_PID[iPid]->Write();
        h_dEdx_vs_mom_PID_plus[iPid]->Write();
        h_dEdx_vs_mom_PID_minus[iPid]->Write();

        if(0)
        {
            h_nSigmaPion   ->Write();
            h_nSigmaKaon   ->Write();
            h_nSigmaProton ->Write();

            h_2D_nSigma_Pion_Kaon   ->Write();
            h_2D_nSigma_Pion_Proton ->Write();
            h_2D_nSigma_Kaon_Proton ->Write();

        }
    }


    h_zdc_total_energy->Write();


//    fTreeIM->Write();


    // ### write DCA:
    for ( int iType = 0; iType < nLabels; iType++ )
    {
        // all pt
        fHist1D_DCA_XY[iType] ->Write();
        fHist1D_DCA_Z[iType] ->Write();
        fHist1D_DCA_Z_SPEC_QA[iType] ->Write();
        fHist2D_DCA_XYvsZ[iType] ->Write();

        fHist1D_DCA_XYZ[iType] ->Write();
        fHist1D_DCA_XYZ_Artif_FB_Sign[iType] ->Write();

        // vs pt:
        fHist2D_DCA_XY_vs_pt[iType] ->Write();
        fHist2D_DCA_Z_vs_pt[iType]  ->Write();
        fHist2D_DCA_XYZ_vs_pt[iType]  ->Write();
        fHist2D_DCA_XYZ_vs_pt_Artif_FB_Sign[iType]  ->Write();
        // //        fOutputList->Add( fHist2D_DCA_Z_SPEC_QA_vs_pt[iType] );
        fHist2D_DCA_XYZ_vs_pt_FB_02_06_Artif_FB_Sign[iType]  ->Write();



        // vs pt, eta:
        fHist3D_DCA_XY_vs_PtEta[iType] ->Write();
        fHist3D_DCA_Z_vs_PtEta[iType] ->Write();
        fHist3D_DCA_XYZ_vs_PtEta[iType] ->Write();
        fHist3D_DCA_XYZ_vs_PtEta_Artif_FB_Sign[iType] ->Write();

    } // end of particle type loop







    // ### WRITE FB:
    for ( int iSub = 0; iSub < nWinPairWithSubsamples; iSub++)
        for ( int iType = 0; iType < nPartTypes; iType++)
            for ( int iCW = 0; iCW < nCW; iCW++)
                for ( int cBin = 0; cBin < nCentrBins[iCW]; ++cBin )
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                        for ( int iEta = 0; iEta < nEtaBins; ++iEta )
                        {
                            WinPairFilling *wp = &winPairs_SIM[iSub][iType][iCW][cBin][iPt][iEta];
                            wp->histAccumulatedValues->Write();

                            wp = &winPairs_REC[iSub][iType][iCW][cBin][iPt][iEta];
                            wp->histAccumulatedValues->Write();

                            wp = &winPairs_REC_PRIM[iSub][iType][iCW][cBin][iPt][iEta];
                            wp->histAccumulatedValues->Write();
                        }




    out.Close();








    // new:
    cout << "counterTracksRecSimChargeMatch = " <<  counterTracksRecSimChargeMatch << endl;
    cout << "counterTracksAllRec = " <<  counterTracksAllRec << endl;
    cout << "counterNoSimFound = " <<  counterNoSimFound << endl;
    cout << "fraction of wrong charge in rec: " <<  1- (double)counterTracksRecSimChargeMatch / counterTracksAllRec << endl;





    TCanvas *canvQA = new TCanvas("canvQA","canvQA",20,20,1200,800);
    canvQA->Divide(2,2);
    canvQA->cd(1)->SetLogz();
    h2D_chargeSimVsRec->DrawCopy( "lego2" );
    canvQA->cd(2)->SetLogz();
    h2D_etaSimVsRec->DrawCopy( "colz" );
    canvQA->cd(3)->SetLogz();
    h2D_ptSimVsRec->DrawCopy( "colz" );


    //
    TCanvas *canvEventQA = new TCanvas( "canvEventQA","canvEventQA", 120,20, 1200,800);
    canvEventQA->Divide(2,2);

    canvEventQA->cd(1);
    //    h1D_vertexZ->DrawCopy();
    canvEventQA->cd(2);
    h1D_vertexZspec->DrawCopy();

}







