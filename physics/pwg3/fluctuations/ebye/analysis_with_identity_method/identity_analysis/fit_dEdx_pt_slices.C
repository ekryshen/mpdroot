// ===========================================================================================================================================================================
// This code takes the file_dEdx_in_bins.root (obtained after get_dEdx_vs_mom_and_true_pi_K_p_moments_from_tree.C)
// and fit narrow momentum slices for each particle species with a function (i.e. Gaus). The results for all slices are then stored into fitted_dEdx_in_pt_bins.root
// Author: Igor Altsybeev (SPbSU), 10.2020
// ===========================================================================================================================================================================


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

#include "../utils.C"

const Int_t fnParticleBins = 6;//4;//2;

const Int_t fnMomBins      = 130;

Int_t fMindEdx = -30000;
Int_t fMaxdEdx = 30000;

int markers[] = { 20, 20, 20, 24, 24, 24 };
int colors[] = { kGreen+2, kBlue, kRed, kGreen+2, kBlue, kRed, };


Double_t Ggauss1(Double_t *x, Double_t *par)
{
    return par[0] * exp ( -pow( fabs( x[0]-par[1])/par[2] , par[3] ) ) * (1+erf(par[4]*(x[0]-par[1])/(par[2]*sqrt(2))))   ;
}


Double_t sumGgauss3(Double_t *x, Double_t *par)
{
    return Ggauss1( x, &par[0] ) + Ggauss1( x, &par[5] ) + Ggauss1( x, &par[10] )   ;
}





static TH1D *****hParticles;


void fit_dEdx_pt_slices()
{
    TFile *fLineShapesLookUpTable = new TFile( "file_dEdx_in_bins.root" );

    TFile *fileWithFittedHistos = new TFile( "fitted_dEdx_in_pt_bins.root", "recreate" );

    //    funcGenGaus->FixParameter(3, 2 );
    //    funcGenGaus->SetParLimits(4, -2, 2);


    TH1D *h_gg_skewness = new TH1D("h_gg_skewness",";Skewness;Entries", 400, -5, 5 );
    TH1D *h_gg_kurtosis = new TH1D("h_gg_kurtosis",";Kurtosis;Entries", 400, -5, 5 );

    TGraphErrors *gr_ampl_vs_p[fnParticleBins];
    TGraphErrors *gr_mean_vs_p[fnParticleBins];
    TGraphErrors *gr_sigma_vs_p[fnParticleBins];
    TGraphErrors *gr_skew_vs_p[fnParticleBins];
    TGraphErrors *gr_kurt_vs_p[fnParticleBins];


    //
    int cBin = 0;
    int eBin = 0;

    int counterGoodFit = 0;
    int counterBadFit = 0;
    int counterNoFit = 0;

    TF1 *funcGenGaus[fnMomBins][fnParticleBins] = {};

    for (Int_t pid = 0; pid<fnParticleBins+2; pid++)
    {
        TH2D *h2D_Particles = 0x0;
        if ( pid < fnParticleBins ) // PIDs
            h2D_Particles = (TH2D*) fLineShapesLookUpTable->Get( Form( "fDeDx_%d", pid ) );
        else if (pid==fnParticleBins) //pos tracks
            h2D_Particles = (TH2D*) fLineShapesLookUpTable->Get( "fDeDx_pos_all" );
        else if (pid==fnParticleBins+1) //neg tracks
            h2D_Particles = (TH2D*) fLineShapesLookUpTable->Get( "fDeDx_neg_all" );


        if( pid < fnParticleBins ) // PID (not All pos and All neg)
        {
            gr_ampl_vs_p[pid]  = new TGraphErrors;
            gr_mean_vs_p[pid]  = new TGraphErrors;
            gr_sigma_vs_p[pid] = new TGraphErrors;
            gr_skew_vs_p[pid]  = new TGraphErrors;
            gr_kurt_vs_p[pid]  = new TGraphErrors;
        }


        //        cout << h2D_Particles << endl;
        for (Int_t ptBin = 0; ptBin<fnMomBins; ptBin++)
        {
            //                                     cout << eBin << " " << cBin << " " << ptBin << " " << pid << endl;
            //                    hParticles[eBin][cBin][ptBin][pid] = (TH1D*) fLineShapesLookUpTable->Get( Form( "particle_eta%d_centr%d_pt%d_pid%d", eBin, cBin, ptBin, pid ) );
            TH1D *proj = h2D_Particles->ProjectionX();
            //                                     cout << proj << endl;
            double lowEdge = proj->GetBinLowEdge(ptBin+1);
            double upEdge = lowEdge + proj->GetBinWidth(ptBin+1);
            //            cout << lowEdge << " " << upEdge << endl;
            TH1D *h = h2D_Particles->ProjectionY( Form( "fDeDx_eta%d_centr%d_pt%d_range_%.2f_%.2f_pid%d", eBin, cBin, ptBin, lowEdge, upEdge, pid ), ptBin+1, ptBin+1 );

            if( pid < fnParticleBins ) // PID (not All pos and All neg):
            {

                // FIT AND REWRITE BINS!
                if ( h->GetEntries() > 10 )//100 )
                {
                    //            cout << h << endl;
                    h->Fit("gaus","Q");
                    TF1 *func = h->GetFunction("gaus");
                    //            fileWithFittedHistos->WriteObject( fileWithFittedHistos, h->GetName() );
                    //            fileWithFittedHistos->cd();




                    if(1 && func ) // have gaus
                    {
                        counterGoodFit++;

                        double *gausParams = func->GetParameters();
                        double A = gausParams[0];
                        double mean = gausParams[1];
                        double sigma = gausParams[2];
                        cout << A << " " << mean << " " << sigma << endl;

                        // Gen gaus:
                        funcGenGaus[ptBin][pid] = new TF1( Form("fGenGaus_%d", pid), Ggauss1, /*fMindEdx, fMaxdEdx,*/ mean-4*sigma, mean+4*sigma, 5 );
                        TF1 *fitFunc = funcGenGaus[ptBin][pid];
                        fitFunc->SetNpx(1500);
                        //                fitFunc->FixParameter(3, 2 );
                        fitFunc->SetParLimits(3, 1, 2 ); // kurtosis
                        fitFunc->SetParLimits(4, -2, 2); // skewness

                        fitFunc->SetParameters( A, mean, sigma, 2, 0 );
                        h->Fit( fitFunc, "Q" );//, "", gausParams[1]-3*gausParams[2], gausParams[1]+3*gausParams[2] );
                        h->Fit( fitFunc, "Q" );//, "", gausParams[1]-3*gausParams[2], gausParams[1]+3*gausParams[2] );
                        //                h->Fit( fitFunc, "Q" );//, "", gausParams[1]-3*gausParams[2], gausParams[1]+3*gausParams[2] );
                        double *ggParams = fitFunc->GetParameters();
                        cout << ggParams[0] << " " << ggParams[1] << " " << ggParams[2] << " " << ggParams[3] << " " << ggParams[4] << endl;

                        h_gg_skewness->Fill( ggParams[4] );
                        h_gg_kurtosis->Fill( ggParams[3] );

                        double momBinCenter = (upEdge+lowEdge)/2;
                        int nP = gr_ampl_vs_p[pid]->GetN();
                        gr_ampl_vs_p[pid]  ->SetPoint( nP, momBinCenter, ggParams[0] );
                        gr_mean_vs_p[pid]  ->SetPoint( nP, momBinCenter, ggParams[1] );
                        gr_sigma_vs_p[pid] ->SetPoint( nP, momBinCenter, ggParams[2] );
                        gr_skew_vs_p[pid]  ->SetPoint( nP, momBinCenter, ggParams[3] );
                        gr_kurt_vs_p[pid]  ->SetPoint( nP, momBinCenter, ggParams[4] );

                        TF1 *func = h->GetFunction( Form("fGenGaus_%d", pid) );
                        if(0 && func)
                            for ( int i = 0; i < h->GetNbinsX(); i++ )
                            {
                                double value = func->Eval( h->GetBinCenter(i+1) );
                                h->SetBinContent( i+1, value );
                            }
                        h->SetTitle( Form( "Skewness %.2f, Kurtosis %.2f", ggParams[4], ggParams[3] ) );

                        //                        delete funcGenGaus;


                        funcGenGaus[ptBin][pid]->SetName( Form( "func_%s", h->GetName() ) );
                    } // end of gen gaus
                    else
                    {
                        counterBadFit++;
                    }

                } // end of gaus
                else
                    counterNoFit++;
            } // end of if ( pid < fnParticleBins )
            else // sum of gen gausses
            {
                cout << "sum of gen gausses for momBin " << ptBin << " and pid " << pid << endl;
                double minForFunc = h->GetBinLowEdge( h->FindFirstBinAbove( 1 ) );
                double maxForFunc = h->GetBinLowEdge( h->FindLastBinAbove( 1 ) );
                //                cout << "minForFunc " << minForFunc << " maxForFunc " << maxForFunc << endl;
                TF1 *funcSum = new TF1( "funcSum", sumGgauss3, /*fMindEdx, fMaxdEdx,*/ minForFunc, maxForFunc, 15 );
                funcSum->SetNpx(2500);
                //                cout << "funcSum " << funcSum << endl;


                int addIfNeg = ( pid == fnParticleBins) ? 0 : fnParticleBins/2;

                TF1 *fPID[3] = {};
                fPID[0] = funcGenGaus[ptBin][0+addIfNeg];
                fPID[1] = funcGenGaus[ptBin][1+addIfNeg];
                fPID[2] = funcGenGaus[ptBin][2+addIfNeg];

                //                cout << " fPID[0] " << fPID[0] << " fPID[1] " << fPID[1]<< " fPID[2] " << fPID[2] << endl;

                //                cout << "fPID[2]->GetParameter(0) = " << fPID[2]->GetParameter(0) << endl;

                if(fPID[0])
                    for (Int_t iPar = 0; iPar<5; iPar++)
                    { funcSum->SetParameter( iPar, fPID[0]->GetParameter(iPar) ); }//cout << "iPar = " << iPar << ", GetParameter(iPar) = " << fPID[0]->GetParameter(iPar) << endl; }
                if(fPID[1])
                    for (Int_t iPar = 0; iPar<5; iPar++)
                    { funcSum->SetParameter( 5+iPar, fPID[1]->GetParameter(iPar) ); }//cout << "iPar = " << iPar << ", GetParameter(iPar) = " << fPID[1]->GetParameter(iPar) << endl; }
                if(fPID[2])
                    for (Int_t iPar = 0; iPar<5; iPar++)
                    { funcSum->SetParameter( 10+iPar, fPID[2]->GetParameter(iPar) ); }//cout << "iPar = " << iPar << ", GetParameter(iPar) = " << fPID[2]->GetParameter(iPar) << endl; }

                //                cout << "before fixing params..." << endl;

                // fix skew:
                if(fPID[0]) funcSum->FixParameter( 4,    fPID[0]->GetParameter(4) );
                if(fPID[1]) funcSum->FixParameter( 5+4,  fPID[1]->GetParameter(4) );
                if(fPID[2]) funcSum->FixParameter( 10+4, fPID[2]->GetParameter(4) );

                // fix kurt:
                if(fPID[0]) funcSum->FixParameter( 3,    fPID[0]->GetParameter(3) );
                if(fPID[1]) funcSum->FixParameter( 5+3,  fPID[1]->GetParameter(3) );
                if(fPID[2]) funcSum->FixParameter( 10+3, fPID[2]->GetParameter(3) );

                //                cout << "before destructing ampl BY HAND..." << endl;

                // destruct ampl BY HAND:
                if(fPID[0]) funcSum->SetParameter( 0,    fPID[0]->GetParameter(0) * gRandom->Gaus(1,0.1) );
                if(fPID[1]) funcSum->SetParameter( 5+0,  fPID[1]->GetParameter(0) * gRandom->Gaus(1,0.1) );
                if(fPID[2]) funcSum->SetParameter( 10+0, fPID[2]->GetParameter(0) * gRandom->Gaus(1,0.1)  );

                // fix ampl = 0 if no entries (=no func) for this PID in this mom slice!
                if(!fPID[0]) funcSum->FixParameter( 0, 0  );
                if(!fPID[1]) funcSum->FixParameter( 5+0, 0  );
                if(!fPID[2]) funcSum->FixParameter( 10+0, 0  );

                //                //                fitFunc->FixParameter(3, 2 );
                //                funcSum->SetParLimits(3, 1, 2 ); // kurtosis
                //                funcSum->SetParLimits(4, -2, 2); // skewness

                //                funcSum->SetParameters( A, mean, sigma, 2, 0 );
                //                cout << "before fit ..." << endl;
                h->Fit( funcSum, "Q" );//, "", gausParams[1]-3*gausParams[2], gausParams[1]+3*gausParams[2] );
                //                cout << "after fit ..." << endl;

                delete funcSum;

            }



            h->Write();

            cout << " end of momBin " << ptBin << " and pid " << pid << endl;
        } // end of mom bins
    } // end of pids

    cout << "test after loops!" << endl;

    //    return;

    h_gg_skewness->Write();
    h_gg_kurtosis->Write();

    for (Int_t pid = 0; pid<fnParticleBins+2; pid++)
        for (Int_t ptBin = 0; ptBin<fnMomBins; ptBin++)
            if(funcGenGaus[ptBin][pid]) funcGenGaus[ptBin][pid]->Write();


    cout << "counterGoodFit = " << counterGoodFit << ", counterBadFit = " << counterBadFit << ", counterNoFit = " << counterNoFit << endl;



    TCanvas *canvQA = new TCanvas("canvQA","canvQA",20,20,1200,800);
    canvQA->Divide(2,2);
    canvQA->cd(1)->SetLogz();
    h_gg_skewness->DrawCopy();
    canvQA->cd(2)->SetLogz();
    h_gg_kurtosis->DrawCopy();


    TCanvas *canv_params = new TCanvas("canv_params","canv_params",120,20,1200,800);
    canv_params->Divide(3,2);
    canv_params->cd(1);
    for (Int_t pid = 0; pid<fnParticleBins; pid++)
        drawGraph( gr_ampl_vs_p[pid], markers[pid], colors[pid], pid == 0 ? "APL" : "PL", 1, 1 );

    canv_params->cd(2);
    for (Int_t pid = 0; pid<fnParticleBins; pid++)
        drawGraph( gr_mean_vs_p[pid], markers[pid], colors[pid], pid == 0 ? "APL" : "PL", 1, 1 );

    canv_params->cd(3);
    for (Int_t pid = 0; pid<fnParticleBins; pid++)
        drawGraph( gr_sigma_vs_p[pid], markers[pid], colors[pid], pid == 0 ? "APL" : "PL", 1, 1 );

    canv_params->cd(4);
    for (Int_t pid = 0; pid<fnParticleBins; pid++)
        drawGraph( gr_skew_vs_p[pid], markers[pid], colors[pid], pid == 0 ? "APL" : "PL", 1, 1 );

    canv_params->cd(5);
    for (Int_t pid = 0; pid<fnParticleBins; pid++)
        drawGraph( gr_kurt_vs_p[pid], markers[pid], colors[pid], pid == 0 ? "APL" : "PL", 1, 1 );
}












