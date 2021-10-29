// =======================================================================================================================================
// This code extracts various kinds of information and an event tree from the MPD DST files with simulations.
// Files for analysis are merged into a chain from the strDirInputDST directory, filenames are taken from inputFilesByHand.txt.
// Author: Igor Altsybeev (SPbSU), Igor.Altsybeev@cern.ch, 10.2020
// =======================================================================================================================================


#define OLD_MPDROOT 0

#include "TH1D.h"
#include "TH2D.h"
#include "THn.h"

#include <TMath.h>
#include <TSystem.h>
#include <TVector3.h>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TInterpreter.h"

#include <iostream>
#include <fstream>

using namespace std;

#include "TDatabasePDG.h"


int get_tree_MPD_task(
        int start_with_file_id = 0
        , int max_file_id_to_read = 3
        , const char * strDirInputDST = "/eos/nica/mpd/sim/data/exp/dst-BiBi-09.5GeV-mp07-20-pwg1-500ev/BiBi/09.5GeV-mb/Smash/BiBi-09.5GeV-mp07-20-pwg1-500ev-1"
        )
{
    Int_t nevents = 0; // number of events to analyse, use 0 to all events

    TString strOutputDir = Form( "/lhep/users/altsybee/test_March_2021_before_commit/output_tree_2020_09_02_SMASH_%d_%d", start_with_file_id, max_file_id_to_read );
    gInterpreter->ProcessLine( Form( ".! mkdir %s", strOutputDir.Data() ) );

    TString outname = Form( "%s/output_with_tree_fId_%d_%d.root", strOutputDir.Data(), start_with_file_id, /*start_with_file_id+max_files_to_read*/max_file_id_to_read );
    cout << "output file name: " << outname << endl;

    TChain chain("mpdsim");

    // ### files to take
//    ifstream base( "/lhep/users/altsybee/analysis_IA/inputFilesByHand.txt" );
    ifstream base( "/lhep/users/altsybee/test_March_2021_before_commit/inputFilesByHand.txt" );


    const int MAX_N_FILES = 8000;
    int file_id_in_folder = 0;

    char *str = new char [1000];

    int nQA_nFiles_taken = 0;
    while (!base.eof())
    {
        base.getline(str, 1000, '\n');
        //        cout << str << endl;
        TString strFull = Form( "%s/%s", strDirInputDST, str );
        if ( TString(strFull).Length() < 3 ) // last line reached
            break;


        if ( file_id_in_folder >= start_with_file_id )
        {
            cout << TString(strFull) << endl;
            chain.Add( TString(strFull) );
            nQA_nFiles_taken++;
        }

        file_id_in_folder++;

        if ( file_id_in_folder == MAX_N_FILES )
            break;
        if ( file_id_in_folder > max_file_id_to_read ) //start_with_file_id + max_files_to_read )
            break;
    }

    cout << "nQA_nFiles_taken = " << nQA_nFiles_taken << endl;

    //    return 0;


    double bImpRanges[] = { 4.69658, 6.61771, 8.08484, 9.32845, 10.4377, 11.427, 12.3386, 13.194, 19.912, };
    int nCbins = sizeof(bImpRanges)/sizeof(*bImpRanges);


    // ##### HISTOS Feb 2020:
    TH1D *h_ImpactPar = new TH1D("h_ImpactPar",";b (fm);Entries",2000,-2,20); //440,-2,20);

    TH1D *h_nTracksInEtaCuts_SIM = new TH1D("h_nTracksInEtaCuts_SIM",";n tracks;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);
    TH1D *h_nTracksInEtaCuts_REC = new TH1D("h_nTracksInEtaCuts_REC",";n tracks;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);

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


    Double_t sigM = 4.0, sigE = 4.0, energy = 7., koef = 1.; /// n-sigma bands for PID selection
    TString Generator = "URQMD", Tracking = "CF";
    MpdPid *pid = new MpdPid(sigM, sigE, energy, koef, Generator, Tracking, "pikapr");


    // ##### IA:
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

    TH1D *h_SIM_Pt_in_eta16 = new TH1D("h_SIM_Pt_in_eta16", "h_SIM_Pt_in_eta16", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_SIM_Eta_in_pt01_10 = new TH1D("h_SIM_Eta_in_pt01_10", "h_SIM_Eta_in_pt01_10", _nYBins, etaMin, etaMax );

    // ### REC:
    TH1D* h_Pt = new TH1D("h_Pt", "h_Pt", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D* h_Eta = new TH1D("h_Eta", "h_Eta", _nYBins, etaMin, etaMax );
    TH1D* h_Phi = new TH1D("h_Phi", "h_Phi", 400, -2*TMath::TwoPi(), 2*TMath::TwoPi());
    TH1D* h_Charge = new TH1D("h_Charge", "h_Charge", 7, -3.5, 3.5);
    //    TH1D* h_Y = new TH1D("h_Y", "h_Y", _nYBins, etaMin, etaMax ); //500, -1, 1);
    TH1D* h_chi2 = new TH1D("h_chi2", "h_chi2", 400, 0, 10 ); //500, -1, 1);
    TH2D* h_VertexXY = new TH2D("h_VertexXY", "h_VertexXY",200, -1, 1, 200, -1, 1 );
    TH1D* h_VertexZ = new TH1D("h_VertexZ", "h_VertexZ", 800, -80, 80);
    TH1D* h_VertexZ_accepted = new TH1D("h_VertexZ_accepted", "h_VertexZ_accepted", 800, -80, 80);

    TH1D *h_Pt_in_eta16 = new TH1D("h_Pt_in_eta16", "h_Pt_in_eta16", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_Eta_in_pt01_10 = new TH1D("h_Eta_in_pt01_10", "h_Eta_in_pt01_10", _nYBins, etaMin, etaMax );

    TH1D* h_nNumberOfHitsTPC = new TH1D("h_nNumberOfHitsTPC", "h_nNumberOfHitsTPC", 201, -0.5, 200.5);
    TH1D* h_nKalmanTrackHits = new TH1D("h_nKalmanTrackHits", "h_nKalmanTrackHits", 201, -0.5, 200.5);
    TH2D* h_2D_nMpdTrack_vs_nKalmanTrackHits = new TH2D("h_2D_nMpdTrack_vs_nKalmanTrackHits", "h_2D_nMpdTrack_vs_nKalmanTrackHits", 201, -0.5, 200.5, 201, -0.5, 200.5);

    TH1D* h_howManyNsimForGivenRecTrack = new TH1D("h_howManyNsimForGivenRecTrack", "h_howManyNsimForGivenRecTrack", 21, -0.5, 20.5);

    // N SIGMA plots:
    TH1D* h_nSigmaPion   = new TH1D("h_nSigmaPion", "h_nSigmaPion", 200, 0, 10 );
    TH1D* h_nSigmaKaon   = new TH1D("h_nSigmaKaon", "h_nSigmaKaon", 200, 0, 10 );
    TH1D* h_nSigmaProton   = new TH1D("h_nSigmaProton", "h_nSigmaProton", 200, 0, 10 );

    TH2D* h_2D_nSigma_Pion_Kaon   = new TH2D("h_2D_nSigma_Pion_Kaon", "h_2D_nSigma_Pion_Kaon", 20, 0, 10, 20, 0, 10 );
    TH2D* h_2D_nSigma_Pion_Proton   = new TH2D("h_2D_nSigma_Pion_Proton", "h_2D_nSigma_Pion_Proton", 20, 0, 10, 20, 0, 10 );
    TH2D* h_2D_nSigma_Kaon_Proton   = new TH2D("h_2D_nSigma_Kaon_Proton", "h_2D_nSigma_Kaon_Proton", 20, 0, 10, 20, 0, 10 );


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




    // ##### EVENT TREE:
    static const Int_t NMaxTrack = 12000;

    Int_t fEventId; //!
    Float_t fImpactPar;   //!
    Float_t fVertexX,fVertexY,fVertexZ; //!
    Float_t f_zdc_total_energy; //!

    TTree       *fEventTree = new TTree("tree_sim_rec","tree_sim_rec"); //! output Tree
    fEventTree->Branch("fEventId",&fEventId,"eventId/I");
    fEventTree->Branch("fImpactPar",&fImpactPar,"impPar/F");

    fEventTree->Branch("vertexX",&fVertexX,"vertexX/F");
    fEventTree->Branch("vertexY",&fVertexY,"vertexY/F");
    fEventTree->Branch("vertexZ",&fVertexZ,"vertexZ/F");

    fEventTree->Branch("zdc_total_energy", &f_zdc_total_energy,"zdc_total_energy/F");

    // SIM:
    Int_t f_SIM_NumberOfTracks; //!
    Int_t f_SIM_TrackId[NMaxTrack];
    Float_t f_SIM_TrackPt[NMaxTrack],f_SIM_TrackEta[NMaxTrack],f_SIM_TrackPhi[NMaxTrack], f_SIM_TrackY[NMaxTrack]; //
    Int_t f_SIM_TrackCharge[NMaxTrack];
    Int_t f_SIM_pid[NMaxTrack];
    Int_t f_SIM_MotherId[NMaxTrack];

    fEventTree->Branch("SIM_NumberOfTracks",&f_SIM_NumberOfTracks,"SIM_NumberOfTracks/I");
    fEventTree->Branch("SIM_TrackId",f_SIM_TrackId,"SIM_TrackId[SIM_NumberOfTracks]/I");

    fEventTree->Branch("SIM_TrackPt",f_SIM_TrackPt,"SIM_TrackPt[SIM_NumberOfTracks]/F");
    fEventTree->Branch("SIM_TrackPhi",f_SIM_TrackPhi,"SIM_TrackPhi[SIM_NumberOfTracks]/F");
    fEventTree->Branch("SIM_TrackEta",f_SIM_TrackEta,"SIM_TrackEta[SIM_NumberOfTracks]/F");
    fEventTree->Branch("SIM_TrackY",f_SIM_TrackY,"SIM_TrackY[SIM_NumberOfTracks]/F");
    fEventTree->Branch("SIM_TrackCharge",f_SIM_TrackCharge,"SIM_TrackCharge[SIM_NumberOfTracks]/I");
    fEventTree->Branch("SIM_pid",f_SIM_pid,"SIM_pid[SIM_NumberOfTracks]/I");
    fEventTree->Branch("SIM_MotherId",f_SIM_MotherId,"SIM_MotherId[SIM_NumberOfTracks]/I");

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
    Int_t fnKalmanHits[NMaxTrack];
    Float_t fTrackChi2[NMaxTrack];

    fEventTree->Branch("nTracks",&fNumberOfTracks,"nTracks/I");
    //    fEventTree->Branch("fSimID",fSimID,"fSimID[nTracks]/I");
    fEventTree->Branch("trackSimID",fTrackSimID,"trackSimID[nTracks]/I");
    fEventTree->Branch("trackSimMotherId",fTrackSimMotherId,"trackSimMotherId[nTracks]/I");

    fEventTree->Branch("trackPt",fTrackPt,"trackPt[nTracks]/F");
    fEventTree->Branch("trackPhi",fTrackPhi,"trackPhi[nTracks]/F");
    fEventTree->Branch("trackEta",fTrackEta,"trackEta[nTracks]/F");
    fEventTree->Branch("trackPtBeforeAbs",fTrackPtBeforeAbs,"trackPtBeforeAbs[nTracks]/F");

    fEventTree->Branch("trackCharge",fTrackCharge,"trackCharge[nTracks]/I");
    //    fEventTree->Branch("trackChargeUsingPt",fTrackChargeUsingPt,"trackChargeUsingPt[nTracks]/I");
    fEventTree->Branch("trackMom",fTrackMom,"trackMom[nTracks]/F");
    fEventTree->Branch("dEdx",fdEdx,"dEdx[nTracks]/F");

    fEventTree->Branch("DCA_x",fDCA_x,"DCA_x[nTracks]/F");
    fEventTree->Branch("DCA_y",fDCA_y,"DCA_y[nTracks]/F");
    fEventTree->Branch("DCA_z",fDCA_z,"DCA_z[nTracks]/F");

    fEventTree->Branch( "nTPChits",fnTPChits,"nTPChits[nTracks]/I" );
    //    fEventTree->Branch("fnKalmanHits",fnKalmanHits,"fnKalmanHits[nTracks]/I" );
    fEventTree->Branch( "trackChi2",fTrackChi2,"trackChi2[nTracks]/F" );



    Float_t fPreviousImpactPar = -1;   // NEED TO REJECT REPEATING EVENTS! (SMASH, first try 100k, Aug 2020)

    //  ### end of IA

    TDatabasePDG *db_PDG = new TDatabasePDG;

    Int_t ev = chain.GetEntries();

    if (nevents == 0) nevents = ev;
    cout << " Number of events in DST file = " << nevents << endl;
    cout << endl;


    //    Bool_t found;
    TClonesArray *mcTracks = NULL;
    TClonesArray *mpdTracks = NULL;

    //  reco event, tracks
    MpdEvent *event = new MpdEvent();
    chain.SetBranchAddress("MCTrack", &mcTracks);
    chain.SetBranchAddress("MPDEvent.", &event);

    TClonesArray *mpdKalmanTracks = NULL;
    chain.SetBranchAddress("TpcKalmanTrack",&mpdKalmanTracks);


    TClonesArray *vertices = NULL;
    chain.SetBranchAddress("Vertex",&vertices);


    // mc particles
#if OLD_MPDROOT == 1
    FairMCTrack *mcTrack = 0x0;//new FairMCTrack;
    FairMCEventHeader *fmcHeader = NULL;
#else
    MpdMCTrack *mcTrack = 0x0;//new MpdMCTrack;
    //MpdMCEventHeader *fmcHeader = NULL;
    FairMCEventHeader *fmcHeader = NULL;
#endif

    MpdTrack *mpdTrack = new MpdTrack;

    chain.SetBranchAddress("MCEventHeader.", &fmcHeader);

    // zdc hits
    TClonesArray *ZDCHits = 0x0;
    MpdZdcDigi* ZDCHit;

    //    cout << "ZDCHits pointer = " << ZDCHits << endl;

    chain.SetBranchAddress("ZdcDigi",&ZDCHits);
    //    cout << "ZDCHits pointer = " << ZDCHits << endl;

    //    return 0;

    double meanMultInRun = 0;
    int runCounter = 0;



    // ####################
    // ### event loop
    for ( int iEv = 0; iEv < nevents; iEv++ )
    {
        if( iEv%100 == 0 )
            cout << "analysing " << (int)iEv << "/" << nevents << " event... \r"; cout.flush();

        chain.GetEntry(iEv);


        double bImp = fmcHeader->GetB();

        // !!! REJECT REPEATING EVENTS! (SMASH, first try 100k, Aug 2020)
        if ( fabs( bImp - fPreviousImpactPar ) < 0.0001 )
            continue;
        fPreviousImpactPar = bImp;


        h_ImpactPar->Fill( bImp );
        //        cout << "fmcHeader->GetB()  = " << fmcHeader->GetB()  << ", fmcHeader->GetEventID()=" << fmcHeader->GetEventID() << endl;

        int cBin = -1;
        for( int j = 0; j < nCbins; j++ )
            if ( bImp < bImpRanges[j] )
            {
                cBin = j;
                break;
            }
        if ( cBin == -1 )
        {
            cout << "AHTUNG!!! cBin == -1" << endl;
            continue;
        }
        int cBinId[10];
        cBinId[0] = cBin; // tmp


        // primary vertex:
        MpdVertex *vertex = (MpdVertex*) vertices->First();
        TVector3 primaryVertex;
        vertex->Position(primaryVertex);

        float x_vertex_mpd = primaryVertex.X();
        float y_vertex_mpd = primaryVertex.Y();
        float z_vertex_mpd = primaryVertex.Z();

        h_VertexXY->Fill( x_vertex_mpd, y_vertex_mpd );
        h_VertexZ->Fill( z_vertex_mpd );

        if ( fabs( z_vertex_mpd ) > 50 ) //40 )
            continue;

        h_VertexZ_accepted->Fill( z_vertex_mpd );



        // ### ZDC info
        //        cout << "ZDCHits pointer = " << ZDCHits << endl;
        Int_t number_of_zdchits = ZDCHits->GetEntries();
        //        cout << "number_of_zdchits = " << number_of_zdchits << endl;
        double zdc_total_energy = 0;
        for (Int_t zdchit_number = 0; zdchit_number < number_of_zdchits; ++zdchit_number)
        {
            ZDCHit = (MpdZdcDigi*) ZDCHits->At(zdchit_number);
            //          Int_t detector_ID = ZDCHit->GetDetectorID();//1,2
            //          Int_t module_ID = ZDCHit->GetModuleID();//1-45
            Double_t energy_deposit_per_hit = ZDCHit->GetELoss();

            zdc_total_energy += energy_deposit_per_hit;
            //          cout << energy_deposit_per_hit << " ";
            //ZDC_energy_mpd[ (detector_ID - 1)*45 + module_ID - 1] += energy_deposit_per_hit;
        }
        //        cout << endl;
        //        cout << "zdc_total_energy = " << zdc_total_energy << endl;
        h_zdc_total_energy->Fill( zdc_total_energy );

        int nSimInCuts = 0;


        // ### FOR TREE:
        fEventId = iEv;
        fImpactPar = bImp;

        fVertexX = x_vertex_mpd;
        fVertexY = y_vertex_mpd;
        fVertexZ = z_vertex_mpd;

        f_zdc_total_energy = zdc_total_energy;

        //        continue;

        int mcNtracks = mcTracks->GetEntries();
        mpdTracks = (TClonesArray*) event->GetGlobalTracks();
        int fNtracks = mpdTracks->GetEntriesFast();

        // arrays of mc flags for reduce mpd matching overhead
        Int_t *flag_array = new Int_t[mcNtracks+1];
        for(Int_t nt = 0; nt <= mcNtracks; nt++)
            flag_array[nt] = 0;

        //        int subsampleId = gRandom->Integer(nWinPairWithSubsamples);


        // ### SIM particle loop
        f_SIM_NumberOfTracks = 0;
        for ( Int_t iP = 0; iP < mcNtracks; iP++ )
        {
#if OLD_MPDROOT == 1
            mcTrack = (FairMCTrack*) mcTracks -> UncheckedAt(iP);
#else
            mcTrack = (MpdMCTrack*) mcTracks -> UncheckedAt(iP);
#endif



            int motherId = mcTrack->GetMotherId();
            if ( motherId != -1) continue;

            //            bool is_pi_K_p = false;
            int pdg = mcTrack->GetPdgCode();
            //            int abs_pdg = abs( pdg );
            //            if ( abs_pdg == 211 || abs_pdg == 321 || abs_pdg == 2212 )
            //                    is_pi_K_p = true;
            //            if (!is_pi_K_p)
            //                continue;

            // ### GET CHARGE:
            //            int charge = (int) partPDG->Charge()/3;
            int charge = 0;

            TParticlePDG *partPDG = db_PDG->GetParticle(pdg);
            if ( partPDG==0 )
            {
                //                cout << "motherId = " << motherId << ", pdg = " << pdg << ", partPDG==0!!!" << endl;
                //                continue;
                // TMP!!!
                if ( pdg == 1000010020 ) // deuteron
                {
                    //                    cout << "deutron" << endl;
                    charge = 1;
                }
                else if ( pdg == -1000010020 ) // anti-deuteron
                    charge = -1;
                else if ( pdg == 1000020040 ) // alpha
                    charge = 2;
                else if ( pdg == -1000020040 ) // anti-alpha
                    charge = -2;
                else if ( pdg == 1000010030 ) // triton
                    charge = 2;
                else if ( pdg == -1000010030 ) // anti-triton
                    charge = -2;
                else // BY HAND!!!
                    charge = 1;
            }
            else
            {
                //                cout << "motherId = " << motherId << ", pdg = " << pdg <<  endl;
                charge = (int) partPDG->Charge()/3;
            }

            //            int charge = (int) partPDG->Charge()/3;
            if ( charge==0 )
                continue;




            double p = mcTrack->GetP(); //if ( (p < 0.1) || (p > 1.6) ) continue;
            double pz = mcTrack->GetPz();
            double eta = 0.5 * TMath::Log( ( p + pz )/( p - pz ) );

            double pt = mcTrack->GetPt();
            //            double eta = mcTrack->GetEta();
            double rapidity = mcTrack->GetRapidity();
            //            double phi = 0;//mcTrack->GetPhi();

            double px = mcTrack->GetPx();
            double py = mcTrack->GetPy();
            double phi = asin( py/pt );
            if (px<0)    phi = TMath::Pi()-phi;
            if ( phi<0 ) phi+=TMath::TwoPi();


            h_SIM_Pt->Fill( pt );
            h_SIM_Eta->Fill( eta );
            h_SIM_Rapidity->Fill( rapidity );
            h_SIM_Phi->Fill( phi );

            // FOR EFFICIENCY PLOTS:
            if( fabs(eta)<1.6 )
                h_SIM_Pt_in_eta16->Fill( pt );
            if( pt > 0.1 && pt < 10 )
                h_SIM_Eta_in_pt01_10->Fill( eta );

            if( fabs(eta)<1.6 && (pt > 0.1 && pt < 10) )
                nSimInCuts++;


            if ( fabs(eta) > 2 && fabs(rapidity) > 2 )
                continue;


            // ### FOR TREE:
            f_SIM_TrackId[f_SIM_NumberOfTracks] = iP;
            f_SIM_TrackPt[f_SIM_NumberOfTracks] = pt;
            f_SIM_TrackPhi[f_SIM_NumberOfTracks] = phi;
            f_SIM_TrackEta[f_SIM_NumberOfTracks] = eta;
            f_SIM_TrackY[f_SIM_NumberOfTracks] = rapidity;
            f_SIM_TrackCharge[f_SIM_NumberOfTracks] = charge;
            f_SIM_pid[f_SIM_NumberOfTracks] = pdg;
            f_SIM_MotherId[f_SIM_NumberOfTracks] = mcTrack->GetMotherId();
            f_SIM_NumberOfTracks++;

            if ( f_SIM_NumberOfTracks == NMaxTrack )
            {
                cout << "AHTUNG! f_SIM_NumberOfTracks = " << f_SIM_NumberOfTracks << "!!! break the loop!" << endl;
                break;
            }



        }
        h_nTracksInEtaCuts_SIM->Fill( nSimInCuts );

        meanMultInRun += nSimInCuts;

        // QA: write av mult in this file (if 1000 events per file!)
        if ( (iEv+1)%1000 == 0 ) // new file will start next event!
        {
            h_QA_meanMultInFiles->SetBinContent( runCounter+1, meanMultInRun/1000 );
            cout << "meanMultInRun/1000 = " << meanMultInRun/1000 << endl;
            meanMultInRun = 0;
            runCounter++;
        }


        int nRecInCuts = 0;

        // ### REC track loop
        fNumberOfTracks = 0;
        for (Int_t trackId = 0; trackId < fNtracks; trackId++)
        {
            mpdTrack = (MpdTrack*) mpdTracks->UncheckedAt(trackId);
            MpdKalmanTrack *kalmanTrack = (MpdKalmanTrack*) mpdKalmanTracks->UncheckedAt(trackId);

            int ID = mpdTrack->GetID();
            flag_array[ID]++;
            if ( flag_array[ID] >= 2 ) // if 2 or more sim for a given rec - skip
                continue;


#if OLD_MPDROOT == 1
            mcTrack = (FairMCTrack*) mcTracks->UncheckedAt(ID);
#else
            mcTrack = (MpdMCTrack*) mcTracks->UncheckedAt(ID);
#endif
            int pdgc = mcTrack->GetPdgCode();
            int mother = mcTrack->GetMotherId();
            double rapidity_MC = mcTrack->GetRapidity();

            double eta = mpdTrack->GetEta();
            double AbsEta = TMath::Abs(eta);
            double phi = mpdTrack->GetPhi();
            if ( phi < 0 )
                phi += TMath::TwoPi();
            double px = mpdTrack->GetPx();
            double py = mpdTrack->GetPy();
            double pz = mpdTrack->GetPz();

            double chi2 = mpdTrack->GetChi2();



            double mom = TMath::Sqrt(px*px + py*py + pz*pz);

            double pt = mpdTrack->GetPt();
            int chargeUsingPt = ( pt >= 0 ? -1 : 1 );
            pt = TMath::Abs(pt);
            double ptBeforeAbs = mpdTrack->GetPt();

            int charge = mpdTrack->GetCharge();

            double dedx = mpdTrack->GetdEdXTPC();
            int tofFlag = mpdTrack->GetTofFlag();

            int nTPCclusters = mpdTrack->GetNofHits();
            h_nNumberOfHitsTPC->Fill( nTPCclusters );


            // from Kalman track:
            int nKalmanHits = kalmanTrack->GetNofHits();
            h_nKalmanTrackHits->Fill( nKalmanHits );
            h_2D_nMpdTrack_vs_nKalmanTrackHits->Fill( nTPCclusters, nKalmanHits );


            int pidFlag = 0; double m2 = 0.;
            if ( (tofFlag == 2) || (tofFlag == 6) ) { pidFlag = 2; m2 = mpdTrack->GetTofMass2(); }

            // Setting cuts
            //            if (mother != -1) continue;
            if (mpdTrack->GetNofHits() < 10 /*20*/ ) continue;
            if (AbsEta > 2.0 /*1.8*//*1.6*/) continue;
            //if (kftrack->GetRecoQuality()) continue;

            // ### IA:
            h_Pt->Fill(pt);                     // plot the pt value of the track in a histogram
            h_Eta->Fill(eta);
            h_Phi->Fill(phi);
            h_Charge->Fill(charge);
            h_chi2->Fill( chi2 );
            //            h_Y->Fill(trackMC->Y());

            // FOR EFFICIENCY PLOTS:
            if( fabs(eta)<1.6 )
                h_Pt_in_eta16->Fill( pt );
            if( pt > 0.1 && pt < 10 )
                h_Eta_in_pt01_10->Fill( eta );


            if( fabs(eta)<1.6 && (pt > 0.1 && pt < 10) )
                nRecInCuts++;


            // ### FOR TREE:
            fTrackSimID[fNumberOfTracks] = ID;
            fTrackSimMotherId[fNumberOfTracks] = mother;

            fTrackPt[fNumberOfTracks] = pt;
            fTrackPhi[fNumberOfTracks] = phi;
            fTrackEta[fNumberOfTracks] = eta;
            fTrackPtBeforeAbs[fNumberOfTracks] = ptBeforeAbs;


            fTrackCharge[fNumberOfTracks] = charge;
            fTrackChargeUsingPt[fNumberOfTracks] = chargeUsingPt;
            fTrackMom[fNumberOfTracks] = mom;
            fdEdx[fNumberOfTracks] = dedx;

            fDCA_x[fNumberOfTracks] = mpdTrack->GetDCAX();
            fDCA_y[fNumberOfTracks] = mpdTrack->GetDCAY();
            fDCA_z[fNumberOfTracks] = mpdTrack->GetDCAZ();

            fnTPChits[fNumberOfTracks] = nTPCclusters;
            fnKalmanHits[fNumberOfTracks] = nKalmanHits;
            fTrackChi2[fNumberOfTracks] = chi2;



            fNumberOfTracks++;


            if ( fNumberOfTracks == NMaxTrack )
            {
                cout << "AHTUNG! fNumberOfTracks = " << fNumberOfTracks << "!!! break the loop!" << endl;
                break;
            }



            bool ret = kFALSE;

            if (( pidFlag == 0) || (pidFlag == 2))
            {
                // FOR nSigma
                ret = pid->FillProbs(mom, dedx, charge); /// only dE/dx available
                if(ret==kTRUE)
                {
                    double nSigmaPion   = pid->GetNsigmaToBetheBloch("pi");
                    double nSigmaKaon   = pid->GetNsigmaToBetheBloch("ka");
                    double nSigmaProton = pid->GetNsigmaToBetheBloch("pr");
                    //            cout << "n sigma: " << nSigmaPion << " " << nSigmaKaon << " " << nSigmaProton << endl;


                    h_nSigmaPion   ->Fill( nSigmaPion   );
                    h_nSigmaKaon   ->Fill( nSigmaKaon   );
                    h_nSigmaProton ->Fill( nSigmaProton );

                    h_2D_nSigma_Pion_Kaon   ->Fill( nSigmaPion,  nSigmaKaon  );
                    h_2D_nSigma_Pion_Proton   ->Fill( nSigmaPion,  nSigmaProton  );
                    h_2D_nSigma_Kaon_Proton   ->Fill( nSigmaKaon,  nSigmaProton  );
                }



                // other fillers:
                h_dEdx->Fill(dedx);
                //                h_dEdx_vs_mom->Fill( mom/charge, dedx );
                //                cout << "dedx = " << dedx << endl;

                h_y_rec_all->Fill( rapidity_MC );
                h_eta_rec_all->Fill( eta );
                h_pt_rec_all->Fill( pt );

                double mom_to_charge = fabs(mom/charge);

                h_dEdx_vs_mom_all->Fill( mom_to_charge, dedx );

                // loop over PIDs:
                if(0)for( int iPid = 0; iPid < nPIDs; iPid++)
                    if ( fabs(pdgc) == arrPidIds[iPid] )
                    {
                        // rapidity
                        h_y_rec_PID[iPid]->Fill( rapidity_MC );
                        if (charge>0)       h_y_rec_PID_plus[iPid]->Fill( rapidity_MC );
                        else if (charge<0)  h_y_rec_PID_minus[iPid]->Fill( rapidity_MC );

                        // eta
                        h_eta_rec_PID[iPid]->Fill( eta );
                        if (charge>0)       h_eta_rec_PID_plus[iPid]->Fill( eta );
                        else if (charge<0)  h_eta_rec_PID_minus[iPid]->Fill( eta );

                        // pt
                        h_pt_rec_PID[iPid]->Fill( pt );
                        if (charge>0)       h_pt_rec_PID_plus[iPid]->Fill( pt );
                        else if (charge<0)  h_pt_rec_PID_minus[iPid]->Fill( pt );

                        // dedx:
                        h_dEdx_vs_mom_PID[iPid]->Fill( mom_to_charge, dedx );
                        if (charge>0)       h_dEdx_vs_mom_PID_plus[iPid]->Fill( mom_to_charge, dedx );
                        else if (charge<0)  h_dEdx_vs_mom_PID_minus[iPid]->Fill( mom_to_charge, dedx );

                        break;
                    }
            } // end of pidFlag == 0

        } // end of kalmanTrack

        // 6.04.2020: count how many times a given recTrack points to a given simTrack
        for(Int_t nt = 0; nt <= mcNtracks; nt++)
            h_howManyNsimForGivenRecTrack->Fill( flag_array[nt] );

        delete [] flag_array;

        h_nTracksInEtaCuts_REC->Fill( nRecInCuts );

        fEventTree->Fill();
    } // end of event loop

    cout << endl;





    // ### write output
    TFile out(outname,"recreate");

    // ##### IA:
    h_ImpactPar->Write();

    h_nTracksInEtaCuts_SIM->Write();
    h_nTracksInEtaCuts_REC->Write();

    h_QA_meanMultInFiles->Write();

    h_VertexXY->Write();
    h_VertexZ->Write();
    h_VertexZ_accepted->Write();

    // ### SIM:
    h_SIM_Pt->Write();
    h_SIM_Eta->Write();
    h_SIM_Rapidity->Write();
    h_SIM_Phi->Write();

    h_SIM_Pt_in_eta16->Write();
    h_SIM_Eta_in_pt01_10->Write();


    // ### REC:
    h_Pt->Write();
    h_Eta->Write();
    h_Phi->Write();
    h_Charge->Write();

    h_Pt_in_eta16->Write();
    h_Eta_in_pt01_10->Write();

    // efficiencies:
    h_Pt_in_eta16->Divide( h_SIM_Pt_in_eta16 );
    h_Pt_in_eta16->Write( "h_eff_vs_pt_in_eta16");

    h_Eta_in_pt01_10->Divide( h_SIM_Eta_in_pt01_10 );
    h_Eta_in_pt01_10->Write( "h_eff_vs_eta_in_pt01_10");



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

    h_nNumberOfHitsTPC->Write();
    h_nKalmanTrackHits->Write();
    h_2D_nMpdTrack_vs_nKalmanTrackHits->Write();

    h_howManyNsimForGivenRecTrack->Write();

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
    }

    h_nSigmaPion   ->Write();
    h_nSigmaKaon   ->Write();
    h_nSigmaProton ->Write();

    h_2D_nSigma_Pion_Kaon   ->Write();
    h_2D_nSigma_Pion_Proton ->Write();
    h_2D_nSigma_Kaon_Proton ->Write();


    h_zdc_total_energy->Write();


    fEventTree->Write();


    out.Close();

    gROOT->ProcessLine(".q");

    return 0;
}

