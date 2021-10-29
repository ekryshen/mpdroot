// =======================================================================================================================================
// This analysis code:
// 1) extracts the dEdx vs pTPC distributions and put them in file_dEdx_in_bins.root
// 2) calculates true (i.e. primary-particle level) moments (first, second, cross-moments) and put them in file_truth_moments_IM.root.
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

#include "utils.C"


const int nSubs = 20;//15;//8;//30;



//
const int nSpecies = 6;

// ### vars:
double _subsN[nSpecies][nSubs] = {};
double _subsN2[nSpecies][nSubs] = {};
double _subsCross[nSpecies*(nSpecies-1)/2][nSubs] = {};



// this ev:
int _n_thisEv[nSpecies] = {};


int _subsNev[nSubs] = {};

void finishEvent(int subsId)
{
    // calc:
    int counterCross = 0;
    for( int i = 0; i < nSpecies; i++ )
    {
        _subsN[i][subsId] += _n_thisEv[i];
        _subsN2[i][subsId] += _n_thisEv[i] * _n_thisEv[i];

        for( int j = i+1; j < nSpecies; j++ )
            _subsCross[counterCross++][subsId] += _n_thisEv[i]*_n_thisEv[j];
    }

    _subsNev[subsId]++;

    // reset:
    for( int i = 0; i < nSpecies; i++ )
        _n_thisEv[i] = 0;
}


double calc_bcorr( double cov, double f, double b, double f2, double b2 )
{
    return (cov-f*b)/sqrt( (f2-f*f) * (b2-b*b) );

}



void calcStdDev( double *arr, double mean, double &std_dev )
{
    std_dev = 0;

    for ( int iSub = 0; iSub < nSubs; iSub++)
    {
        float diff = arr[iSub] - mean;
        std_dev += diff*diff;
    }

    // stdDev of the mean:
    std_dev /= (nSubs-1); // -1 is important. see also Marek's lecture - https://indico.cern.ch/event/722562/
    std_dev = sqrt(std_dev);
    std_dev = std_dev / sqrt(nSubs);
}


// #############################################################
void get_dEdx_vs_mom_and_true_pi_K_p_moments_from_tree()
{
    TFile *fileInput = new TFile( "an_results_40_file_with_IM_tree_dEdx.root" );
    TTree       *fTreeIM = (TTree*) fileInput->Get( "treeIM" );
    cout << "n entries: " << fTreeIM->GetEntries() << endl;

    ULong64_t IM_event; //!
    //    Float_t IM_vZ; //!
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

    fTreeIM->SetBranchAddress("event",        &IM_event );
    fTreeIM->SetBranchAddress("centrality",   &IM_centrality );
    fTreeIM->SetBranchAddress("pt",   &IM_pt );
    fTreeIM->SetBranchAddress("mom",   &IM_mom );
    fTreeIM->SetBranchAddress("dEdx",   &IM_dEdx );
    fTreeIM->SetBranchAddress("eta",   &IM_eta );

    fTreeIM->SetBranchAddress("sim_PID",   &IM_sim_PID );


    //
    int nBins_mom = 130;//500;
    double momMin = 0.2;// 0.;
    double momMax = 1.5;//2;
    int nBins_dEdx = 300*5;//280*5; //980*2;
    double dEdxMin = 0;//2000;//20.;
    double dEdxMax = 30000;//1000;
    TH2D *fDeDx_AllTracks = new TH2D( "fDeDx_AllTracks","all tracks;p_{TPC} (GeV/c);dE/dx (a.u.)", nBins_mom, momMin, momMax, nBins_dEdx, dEdxMin, dEdxMax);
    TH2D *fDeDx_el = new TH2D( "fDeDx_el",  "electrons;p_{TPC} (GeV/c);dE/dx (a.u.)", nBins_mom, momMin, momMax, nBins_dEdx, dEdxMin, dEdxMax);
    TH2D *fDeDx_pi = new TH2D( "fDeDx_pi",  "pions;p_{TPC} (GeV/c);dE/dx (a.u.)", nBins_mom, momMin, momMax, nBins_dEdx, dEdxMin, dEdxMax);
    TH2D *fDeDx_K  = new TH2D( "fDeDx_K",   "kaons;p_{TPC} (GeV/c);dE/dx (a.u.)",  nBins_mom, momMin, momMax, nBins_dEdx, dEdxMin, dEdxMax);
    TH2D *fDeDx_p  = new TH2D( "fDeDx_p",   "protons;p_{TPC} (GeV/c);dE/dx (a.u.)",  nBins_mom, momMin, momMax, nBins_dEdx, dEdxMin, dEdxMax);
    //    TH2D *fDeDx_p_SPEC_ONLY_TPC_nSigma  = new TH2D( "fDeDx_p_SPEC_ONLY_TPC_nSigma",   "fDeDx_p_SPEC_ONLY_TPC_nSigma",  nBins_mom, momMin, momMax, nBins_dEdx, dEdxMin, dEdxMax);

    // pos
    TH2D *fDeDx_pos_AllTracks = new TH2D( "fDeDx_pos_AllTracks","all pos tracks;p_{TPC} (GeV/c);dE/dx (a.u.)", nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);
    TH2D *fDeDx_pos_el = new TH2D( "fDeDx_pos_el",  "electrons;p_{TPC} (GeV/c);dE/dx (a.u.)",      nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);
    TH2D *fDeDx_pos_pi = new TH2D( "fDeDx_pos_pi",  "pions;p_{TPC} (GeV/c);dE/dx (a.u.)",    nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);
    TH2D *fDeDx_pos_K  = new TH2D( "fDeDx_pos_K",   "kaons;p_{TPC} (GeV/c);dE/dx (a.u.)",    nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);
    TH2D *fDeDx_pos_p  = new TH2D( "fDeDx_pos_p",   "protons;p_{TPC} (GeV/c);dE/dx (a.u.)",  nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);

    // neg
    TH2D *fDeDx_neg_AllTracks = new TH2D( "fDeDx_neg_AllTracks","all neg tracks;p_{TPC} (GeV/c);dE/dx (a.u.)", nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);
    TH2D *fDeDx_neg_el = new TH2D( "fDeDx_neg_el",  "anti electrons;p_{TPC} (GeV/c);dE/dx (a.u.)", nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);
    TH2D *fDeDx_neg_pi = new TH2D( "fDeDx_neg_pi",  "anti pions;p_{TPC} (GeV/c);dE/dx (a.u.)",     nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);
    TH2D *fDeDx_neg_K  = new TH2D( "fDeDx_neg_K",   "anti kaons;p_{TPC} (GeV/c);dE/dx (a.u.)",     nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);
    TH2D *fDeDx_neg_p  = new TH2D( "fDeDx_neg_p",   "anti protons;p_{TPC} (GeV/c);dE/dx (a.u.)",   nBins_mom, momMin, momMax, 2*nBins_dEdx, -dEdxMax, dEdxMax);


    // #### event loop
    ULong64_t prevEvId = -1;
    int subsampleId = gRandom->Integer(nSubs);
    for( int iEv = 0; iEv < fTreeIM->GetEntries(); iEv++ )
    {
        fTreeIM->GetEntry(iEv);

        if(iEv == 0)
            prevEvId = IM_event;

        if( IM_event != prevEvId )
        {
            finishEvent(subsampleId);
            // assign new subsample id
            subsampleId = gRandom->Integer(nSubs);
            prevEvId = IM_event;
        }


        if( fabs(IM_eta) > 0.5 )//0.5 )
            continue;



        // pt cuts:
        if ( fabs(IM_mom) > 0.3 && fabs(IM_mom) < 1.5 )
        {
            // FB analysis:
            if( IM_eta > 0.1 )
            {
                if( IM_sim_PID == 211        ) _n_thisEv[0]++;
                else if( IM_sim_PID == 321   ) _n_thisEv[1]++;
                else if( IM_sim_PID == 2212  ) _n_thisEv[2]++;
            }
            else if( IM_eta < -0.1 )
            {
                if( IM_sim_PID == -211  ) _n_thisEv[3]++;
                else if( IM_sim_PID == -321  ) _n_thisEv[4]++;
                else if( IM_sim_PID == -2212 ) _n_thisEv[5]++;
            }
        }


        fDeDx_AllTracks->Fill( IM_mom, IM_dEdx );
        if ( IM_sim_PID>0 ) fDeDx_pos_AllTracks->Fill( IM_mom, IM_dEdx );
        else fDeDx_neg_AllTracks->Fill( IM_mom, IM_dEdx );


        if ( fabs(IM_sim_PID) == 11 )
        {
            fDeDx_el->Fill( IM_mom, IM_dEdx );
            if ( IM_sim_PID == 11 ) fDeDx_pos_el->Fill( IM_mom, IM_dEdx );
            else fDeDx_neg_el->Fill( IM_mom, IM_dEdx );
        }

        if ( fabs(IM_sim_PID) == 211 )
        {
            fDeDx_pi->Fill( IM_mom, IM_dEdx );
            if ( IM_sim_PID == 211 ) fDeDx_pos_pi->Fill( IM_mom, IM_dEdx );
            else fDeDx_neg_pi->Fill( IM_mom, IM_dEdx );
        }

        if ( fabs(IM_sim_PID) == 321 )
        {
            fDeDx_K->Fill( IM_mom, IM_dEdx );
            if ( IM_sim_PID == 321 ) fDeDx_pos_K->Fill( IM_mom, IM_dEdx );
            else fDeDx_neg_K->Fill( IM_mom, IM_dEdx );

        }

        if ( fabs(IM_sim_PID) == 2212 )
        {
            fDeDx_p->Fill( IM_mom, IM_dEdx );
            if ( IM_sim_PID == 2212 ) fDeDx_pos_p->Fill( IM_mom, IM_dEdx );
            else fDeDx_neg_p->Fill( IM_mom, IM_dEdx );
        }
    }
    // finish last event:
    finishEvent(subsampleId);


    TCanvas *canv_dEdx = new TCanvas("canv_dEdx","canv dE/dx",20,20,1200,600 );
    canv_dEdx->Divide(3,2);

    canv_dEdx->cd(1);
    fDeDx_AllTracks->DrawCopy("colz");
    canv_dEdx->cd(2);
    fDeDx_el->DrawCopy("colz");
    canv_dEdx->cd(3);
    fDeDx_pi->DrawCopy("colz");
    canv_dEdx->cd(4);
    fDeDx_K->DrawCopy("colz");
    canv_dEdx->cd(5);
    fDeDx_p->DrawCopy("colz");

    // pos
    TCanvas *canv_dEdx_pos = new TCanvas("canv_dEdx_pos","canv dE/dx pos",80,40,800,800 );
    canv_dEdx_pos->Divide(2,2);

    canv_dEdx_pos->cd(1)->SetLogz();
    fDeDx_pos_AllTracks->DrawCopy("colz");
    //    canv_dEdx_pos->cd(2);
    //    fDeDx_pos_el->DrawCopy("colz");
    canv_dEdx_pos->cd(2)->SetLogz();
    fDeDx_pos_pi->DrawCopy("colz");
    canv_dEdx_pos->cd(3)->SetLogz();
    fDeDx_pos_K->DrawCopy("colz");
    canv_dEdx_pos->cd(4)->SetLogz();
    fDeDx_pos_p->DrawCopy("colz");

    // neg
    TCanvas *canv_dEdx_neg = new TCanvas("canv_dEdx_neg","canv dE/dx neg",120,60,800,800 );
    canv_dEdx_neg->Divide(2,2);

    canv_dEdx_neg->cd(1)->SetLogz();
    fDeDx_neg_AllTracks->DrawCopy("colz");
    //    canv_dEdx_neg->cd(2);
    //    fDeDx_neg_el->DrawCopy("colz");
    canv_dEdx_neg->cd(2)->SetLogz();
    fDeDx_neg_pi->DrawCopy("colz");
    canv_dEdx_neg->cd(3)->SetLogz();
    fDeDx_neg_K->DrawCopy("colz");
    canv_dEdx_neg->cd(4)->SetLogz();
    fDeDx_neg_p->DrawCopy("colz");


    // ### write
    TFile *file_dEdx_in_bins = new TFile( "file_dEdx_in_bins.root", "recreate" );

    fDeDx_AllTracks->Write();
    fDeDx_el->Write();
    fDeDx_pi->Write();
    fDeDx_K->Write();
    fDeDx_p->Write();

    // pos
    fDeDx_pos_AllTracks->Write("fDeDx_pos_all");
    fDeDx_pos_el->Write();
    fDeDx_pos_pi->Write("fDeDx_0");
    fDeDx_pos_K->Write("fDeDx_1");
    fDeDx_pos_p->Write("fDeDx_2");

    // neg
    fDeDx_neg_AllTracks->Write("fDeDx_neg_all");
    fDeDx_neg_el->Write();
    fDeDx_neg_pi->Write("fDeDx_3");
    fDeDx_neg_K->Write("fDeDx_4");
    fDeDx_neg_p->Write("fDeDx_5");


    // averaging:
    double avN[nSpecies] = {};
    double avN2[nSpecies] = {};
    double avCross[nSpecies*(nSpecies-1)/2] = {};


    int nEvents = 0;
    for( int subId = 0; subId < nSubs; subId++ )
        nEvents += _subsNev[subId];

    int counterCross = 0;
    for( int i = 0; i < nSpecies; i++ )
    {
        for( int subId = 0; subId < nSubs; subId++ )
        {
            avN[i] += _subsN[i][subId];
            avN2[i] += _subsN2[i][subId];

            _subsN[i][subId] /= _subsNev[subId];
            _subsN2[i][subId] /= _subsNev[subId];
        }
        avN[i] /= nEvents;
        avN2[i] /= nEvents;

        // cross moments:
        for( int j = i+1; j < nSpecies; j++ )
        {
            for( int subId = 0; subId < nSubs; subId++ )
            {
                avCross[counterCross] += _subsCross[counterCross][subId];
                _subsCross[counterCross][subId] /= _subsNev[subId];
            }
            avCross[counterCross] /= nEvents;
            counterCross++;
        }
    }



    // ### bcorr
    double av_bcorr[nSpecies*(nSpecies-1)/2] = {};
    double _subs_bcorr[nSpecies*(nSpecies-1)/2][nSubs] = {};

    counterCross = 0;
    for( int i = 0; i < nSpecies; i++ )
        for( int j = i+1; j < nSpecies; j++ )
        {
            for( int subId = 0; subId < nSubs; subId++ )
                _subs_bcorr[counterCross][subId] = calc_bcorr( _subsCross[counterCross][subId],  _subsN[i][subId], _subsN[j][subId], _subsN2[i][subId], _subsN2[j][subId] );

            av_bcorr[counterCross] = calc_bcorr( avCross[counterCross],  avN[i], avN[j], avN2[i], avN2[j] );
            counterCross++;
        }





    double std_dev_N[nSpecies] = {};
    double std_dev_N2[nSpecies] = {};
    double std_dev_Cross[nSpecies*(nSpecies-1)/2] = {};
    double std_dev_bcorr[nSpecies*(nSpecies-1)/2] = {};
    counterCross = 0;
    for( int i = 0; i < nSpecies; i++ )
    {
        calcStdDev( _subsN[i], avN[i], std_dev_N[i] );
        calcStdDev( _subsN2[i], avN2[i], std_dev_N2[i] );
        for( int j = i+1; j < nSpecies; j++ )
        {
            calcStdDev( _subsCross[counterCross], avCross[counterCross], std_dev_Cross[counterCross] );
            calcStdDev( _subs_bcorr[counterCross], av_bcorr[counterCross], std_dev_bcorr[counterCross] );

            counterCross++;
        }
    }



    // histos with moments
    TH1D *histFirstMom = new TH1D( "histFirstMom","; ;Entries",nSpecies,-0.5,nSpecies-0.5);
    TString arrFirstMomNames[nSpecies] = {"#LT #pi^{+} #GT", "#LT K^{+} #GT", "#LT p #GT", "#LT #pi^{#minus} #GT", "#LT K^{#minus} #GT", "#LT #bar{p} #GT" };
    for(Int_t i = 0; i < nSpecies; i++)
        histFirstMom->GetXaxis()->SetBinLabel( i+1, arrFirstMomNames[i].Data() );
    tuneHist1D(histFirstMom);

    TH1D *histSecondMom = new TH1D( "histSecondMom","; ;Entries", nSpecies, -0.5, nSpecies-0.5 );
    TString arrSecMomNames[nSpecies] = {"#LT #pi^{+ 2} #GT", "#LT K^{+ 2} #GT", "#LT p^{2} #GT", "#LT #pi^{#minus 2} #GT", "#LT K^{#minus 2} #GT", "#LT #bar{p}^{2} #GT" };
    for(Int_t i = 0; i < nSpecies; i++)
        histSecondMom->GetXaxis()->SetBinLabel( i+1, arrSecMomNames[i].Data() );
    tuneHist1D(histSecondMom);

    const int nComb = nSpecies*(nSpecies-1)/2;
    TH1D *histMixedMom = new TH1D( "histMixedMom","; ;Entries", nComb, -0.5, nComb-0.5 );
    TString arrMixedMomNames[nComb] = { "#LT #pi^{+}K^{+} #GT", "#LT #pi^{+}p #GT", "#LT #pi^{+}#pi^{#minus} #GT", "#LT #pi^{+}K^{#minus} #GT", "#LT #pi^{+}#bar{p} #GT",
                                        "#LT K^{+}p #GT",   "#LT K^{+}#pi^{#minus} #GT", "#LT K^{+}K^{#minus} #GT", "#LT K^{+}#bar{p} #GT",
                                        "#LT p#pi^{#minus} #GT", "#LT pK^{#minus} #GT", "#LT p#bar{p} #GT",
                                        "#LT #pi^{#minus}K^{#minus} #GT", "#LT #pi^{#minus}#bar{p} #GT",
                                        "#LT K^{#minus}#bar{p} #GT",
                                      };



    for(Int_t i = 0; i < nComb; i++)
        histMixedMom->GetXaxis()->SetBinLabel( i+1, arrMixedMomNames[i].Data() );
    tuneHist1D(histMixedMom);



    // bcorr:
    TH1D *hist_bcorr = new TH1D( "hist_bcorr","; ;Entries", nComb, -0.5, nComb-0.5 );
    TString arrBcorrBinsNames[nComb] = { "#pi^{+}K^{+}", "#pi^{+}p", "#pi^{+}#pi^{#minus}", "#pi^{+}K^{#minus}", "#pi^{+}#bar{p}",
                                         "K^{+}p",   "K^{+}#pi^{#minus}", "K^{+}K^{#minus}", "K^{+}#bar{p}",
                                         "p#pi^{#minus}", "pK^{#minus}", "p#bar{p}",
                                         "#pi^{#minus}K^{#minus}", "#pi^{#minus}#bar{p}",
                                         "K^{#minus}#bar{p}",
                                       };
    for(Int_t i = 0; i < nComb; i++)
        hist_bcorr->GetXaxis()->SetBinLabel( i+1, arrBcorrBinsNames[i].Data() );
    tuneHist1D(hist_bcorr);



    histFirstMom->SetMarkerStyle(25);
    histFirstMom->SetMarkerColor(kBlue);
    histFirstMom->SetLineColor(kBlue);

    histSecondMom->SetMarkerStyle(25);
    histSecondMom->SetMarkerColor(kBlue);
    histSecondMom->SetLineColor(kBlue);

    histMixedMom->SetMarkerStyle(25);
    histMixedMom->SetMarkerColor(kBlue);
    histMixedMom->SetLineColor(kBlue);


    hist_bcorr->SetMarkerStyle(25);
    hist_bcorr->SetMarkerColor(kRed);
    hist_bcorr->SetLineColor(kRed);




    // ### Print & Draw:

    // ##### 1st Moments
    for( int i = 0; i < nSpecies; i++ )
    {
        cout << " First  --> " << i << " = " <<  avN[i] << " +/- " << std_dev_N[i] << endl;
        histFirstMom->SetBinContent(i+1, avN[i]);
        histFirstMom->SetBinError(i+1, std_dev_N[i]);
    }
    TCanvas *canv_FirstMom = new TCanvas("canv_FirstMom","FirstMom",120,80,800,600 );
    tuneCanvas( canv_FirstMom );
    histFirstMom->DrawCopy();




    // ##### 2nd Moments
    for( int i = 0; i < nSpecies; i++ )
    {
        cout << " Second --> " << i << " = " <<  avN2[i] << " +/- " << std_dev_N2[i] << endl;
        histSecondMom->SetBinContent(i+1, avN2[i]);
        histSecondMom->SetBinError(i+1, std_dev_N2[i]);
    }
    TCanvas *canv_SecMom = new TCanvas("canv_SecMom","SecMom",150,110,800,600 );
    tuneCanvas( canv_SecMom );
    histSecondMom->DrawCopy();


    // ##### Mixed Moments
    counterCross = 0;
    for( int i = 0; i < nSpecies; i++ )
        for( int j = i+1; j < nSpecies; j++ )
        {
            cout << " Mixed --> " << i << "-" << j << " = " <<  avCross[counterCross] << " +/- " << std_dev_Cross[counterCross] << endl;
            histMixedMom->SetBinContent(counterCross+1, avCross[counterCross] );
            histMixedMom->SetBinError(counterCross+1, std_dev_Cross[counterCross] );

            counterCross++;
        }
    TCanvas *canv_MixedMom = new TCanvas("canv_MixedMom","MixedMom",170,120,800,600 );
    tuneCanvas( canv_MixedMom );
    histMixedMom->DrawCopy();



    counterCross = 0;
    for (Int_t i=0; i<nSpecies; i++){
        for (Int_t j=0; j<nSpecies; j++){
            if (i>=j) continue;
            cout << " b_corr --> " << i << "-" << j << " = " <<  av_bcorr[counterCross] << " +/- " <<  std_dev_bcorr[counterCross] << endl;
            hist_bcorr->SetBinContent(counterCross+1, av_bcorr[counterCross] );
            hist_bcorr->SetBinError(counterCross+1, std_dev_bcorr[counterCross] );
            counterCross++;
        }
    }


    cout << "#########" << endl;

    for (Int_t i=0; i<nSpecies; i++)
    {
        double Var = avN2[i] - avN[i]*avN[i];
        cout << "PID " << i << ": Var = " << Var << endl;
    }

    counterCross = 0;
    for (Int_t i=0; i<nSpecies; i++)
    {
        for (Int_t j=i+1; j<nSpecies; j++)
        {
            double VarA = avN2[i] - avN[i]*avN[i];
            double VarB = avN2[j] - avN[j]*avN[j];
            double COV_FB = avCross[counterCross] - avN[i] * avN[j];
            double bcorr = COV_FB / sqrt(VarA*VarB);
            cout << "COV_FB (" << i << ", " << j << ") = " << COV_FB << ", bcorr = " << bcorr << ", bcorr (recalc) = " << av_bcorr[counterCross] << " +/- " << std_dev_bcorr[counterCross] << endl;

            counterCross++;
        }

    }

    cout << endl;
    cout << "nEvents = " << nEvents << endl;



    // ### write results (moments)
    TFile *file_truth_moments_IM = new TFile( "file_truth_moments_IM.root", "recreate" );
    histFirstMom->Write();
    histSecondMom->Write();
    histMixedMom->Write();
    hist_bcorr->Write();
    file_truth_moments_IM->Close();
}














