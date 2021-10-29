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


const int nSubs = 20;//15;//8;//30;
const int nSpecies = 6;

void copyErrors( TH1D *hBase, TH1D *hToWrite )
{
    for ( int i = 0; i < hBase->GetNbinsX(); i++ )
    {
        double err = hBase->GetBinError(i+1);
        hToWrite->SetBinError( i+1, err );
    }
}


void draw_comparison_truth_IM()
{

    TFile *fileTruth = new TFile( "file_truth_moments_IM.root" );
    TFile *fileExtracted = new TFile( "file_extracted_moments_IM.root" );

    TH1D *histFirstMom_truth = (TH1D*)fileTruth->Get( "histFirstMom" );
    TH1D *histSecondMom_truth = (TH1D*)fileTruth->Get( "histSecondMom" );
    TH1D *histMixedMom_truth = (TH1D*)fileTruth->Get( "histMixedMom" );

    TH1D *histFirstMom_IM = (TH1D*)fileExtracted->Get( "histFirstMom" );
    TH1D *histSecondMom_IM = (TH1D*)fileExtracted->Get( "histSecondMom" );
    TH1D *histMixedMom_IM = (TH1D*)fileExtracted->Get( "histMixedMom" );

    gStyle->SetOptStat(false);

    histFirstMom_truth->SetMarkerStyle(25);
    histFirstMom_truth->SetMarkerColor(kBlue+1);
    histFirstMom_truth->SetLineColor(kBlue+1);
    histFirstMom_truth->SetLineWidth(2);
    histFirstMom_truth->SetMarkerSize(1.8);


    histFirstMom_IM->SetMarkerStyle(20);
    histFirstMom_IM->SetMarkerColor(kRed);
    histFirstMom_IM->SetLineColor(kRed);
    histFirstMom_IM->SetLineWidth(2);
    histFirstMom_IM->SetMarkerSize(1.8);
    histFirstMom_IM->GetYaxis()->SetTitleSize(0.07);
    histFirstMom_IM->GetYaxis()->SetTitleOffset(1);
    histFirstMom_IM->GetYaxis()->CenterTitle();
    histFirstMom_IM->GetYaxis()->SetLabelSize(0.055);
    histFirstMom_IM->GetXaxis()->SetLabelSize(0.078);


    //
    histSecondMom_truth->SetMarkerStyle(25);
    histSecondMom_truth->SetMarkerColor(kBlue+1);
    histSecondMom_truth->SetLineColor(kBlue+1);
    histSecondMom_truth->SetLineWidth(2);
    histSecondMom_truth->SetMarkerSize(1.8);

    histSecondMom_IM->SetMarkerStyle(20);
    histSecondMom_IM->SetMarkerColor(kRed);
    histSecondMom_IM->SetLineColor(kRed);
    histSecondMom_IM->SetLineWidth(2);
    histSecondMom_IM->SetMarkerSize(1.8);
    histSecondMom_IM->GetYaxis()->SetTitleSize(0.07);
    histSecondMom_IM->GetYaxis()->SetTitleOffset(1);
    histSecondMom_IM->GetYaxis()->CenterTitle();
    histSecondMom_IM->GetYaxis()->SetLabelSize(0.055);
    histSecondMom_IM->GetXaxis()->SetLabelSize(0.078);

    //
    histMixedMom_truth->SetMarkerStyle(25);
    histMixedMom_truth->SetMarkerColor(kBlue+1);
    histMixedMom_truth->SetLineColor(kBlue+1);
    histMixedMom_truth->SetLineWidth(2);
    histMixedMom_truth->SetMarkerSize(1.8);
    histMixedMom_IM->SetMarkerStyle(20);
    histMixedMom_IM->SetMarkerColor(kRed);
    histMixedMom_IM->SetLineColor(kRed);
    histMixedMom_IM->SetLineWidth(2);
    histMixedMom_IM->SetMarkerSize(1.8);
    histMixedMom_IM->GetYaxis()->SetTitleSize(0.07);
    histMixedMom_IM->GetYaxis()->SetTitleOffset(1);
    histMixedMom_IM->GetYaxis()->CenterTitle();
    histMixedMom_IM->GetYaxis()->SetLabelSize(0.055);
    histMixedMom_IM->GetXaxis()->SetLabelSize(0.053);




    TCanvas *canv_FirstMom = new TCanvas("canv_FirstMom","FirstMom",120,80,800,600 );
    tuneCanvas( canv_FirstMom );
    histFirstMom_IM->GetYaxis()->SetTitle( "#LTN#GT" );
    histFirstMom_IM->DrawCopy("P");
    histFirstMom_truth->DrawCopy("same");


    TLegend *leg_n = new TLegend(0.51,0.4,0.9,0.67);
    tuneLegend( leg_n );
    leg_n->AddEntry( histFirstMom_IM, "Identity Method", "p");
    leg_n->AddEntry( histFirstMom_truth, "truth", "p");
    leg_n->Draw();


    TCanvas *canv_SecMom = new TCanvas("canv_SecMom","SecMom",150,110,800,600 );
    tuneCanvas( canv_SecMom );
    histSecondMom_IM->GetYaxis()->SetTitle("#LT N^{2} #GT");
    histSecondMom_IM->DrawCopy("P");
    histSecondMom_truth->DrawCopy("same");

    TLegend *leg_n2 = new TLegend(0.51,0.4,0.9,0.67);
    tuneLegend( leg_n2 );
    leg_n2->AddEntry( histSecondMom_IM, "Identity Method", "p");
    leg_n2->AddEntry( histSecondMom_truth, "truth", "p");
    leg_n2->Draw();



    TCanvas *canv_MixedMom = new TCanvas("canv_MixedMom","MixedMom",170,120,800,600 );
    tuneCanvas( canv_MixedMom );
    histMixedMom_IM->DrawCopy("P");
    histMixedMom_truth->DrawCopy("same");



    // ### RATIOS:
    copyErrors( histFirstMom_truth, histFirstMom_IM );
    copyErrors( histSecondMom_truth, histSecondMom_IM );
    copyErrors( histMixedMom_truth, histMixedMom_IM );


    // #####
    TCanvas *canv_FirstMom_RATIO = new TCanvas("canv_FirstMom_RATIO","FirstMom ratio",120,80,800,600 );
    tuneCanvas( canv_FirstMom_RATIO );
    canv_FirstMom_RATIO->SetGrid();
    histFirstMom_IM->Divide( histFirstMom_truth );
    histFirstMom_IM->GetYaxis()->SetNdivisions(507);
    histFirstMom_IM->GetYaxis()->SetRangeUser(0.955, 1.045);
    histFirstMom_IM->GetYaxis()->SetTitle("IM / truth");
    histFirstMom_IM->DrawCopy();


    TLine *line_unity = new TLine(-0.5,1,histFirstMom_IM->GetNbinsX()-0.5,1);
    line_unity->SetLineWidth(2);
    line_unity->SetLineStyle(9);
    line_unity->SetLineColor(kGray+1);
    line_unity->DrawClone();

    histFirstMom_IM->DrawCopy("same");


    // ####
    TCanvas *canv_SecMom_RATIO = new TCanvas("canv_SecMom_RATIO","SecMom ratio",120,80,800,600 );
    tuneCanvas( canv_SecMom_RATIO );
    canv_SecMom_RATIO->SetGrid();
    histSecondMom_IM->Divide( histSecondMom_truth );
    histSecondMom_IM->GetYaxis()->SetNdivisions(507);
    histSecondMom_IM->GetYaxis()->SetRangeUser(0.87, 1.13);
    histSecondMom_IM->GetYaxis()->SetTitle("IM / truth");
    histSecondMom_IM->DrawCopy();

    line_unity = new TLine(-0.5,1,histSecondMom_IM->GetNbinsX()-0.5,1);
    line_unity->SetLineWidth(2);
    line_unity->SetLineStyle(9);
    line_unity->SetLineColor(kGray+1);
    line_unity->DrawClone();

    histSecondMom_IM->DrawCopy("same");

    // ####
    TCanvas *canv_MixedMom_RATIO = new TCanvas("canv_MixedMom_RATIO","MixedMom ratio",120,80,800,600 );
    tuneCanvas( canv_MixedMom_RATIO );
    canv_MixedMom_RATIO->SetGrid();
    histMixedMom_IM->Divide( histMixedMom_truth );
    histMixedMom_IM->GetYaxis()->SetNdivisions(507);
    histMixedMom_IM->GetYaxis()->SetRangeUser(0.87, 1.13);
    histMixedMom_IM->GetYaxis()->SetTitle("IM / truth");

    // FOR FB:
    if(0)
    {
        int binId = 1; histMixedMom_IM->SetBinContent( binId,-1);      histMixedMom_IM->GetXaxis()->SetBinLabel( binId,"");
        binId = 2;  histMixedMom_IM->SetBinContent( binId,-1);      histMixedMom_IM->GetXaxis()->SetBinLabel( binId,"");
        binId = 6;  histMixedMom_IM->SetBinContent( binId,-1);      histMixedMom_IM->GetXaxis()->SetBinLabel( binId,"");
        binId = 13;  histMixedMom_IM->SetBinContent( binId,-1);      histMixedMom_IM->GetXaxis()->SetBinLabel( binId,"");
        binId = 14;  histMixedMom_IM->SetBinContent( binId,-1);      histMixedMom_IM->GetXaxis()->SetBinLabel( binId,"");
        binId = 15;  histMixedMom_IM->SetBinContent( binId,-1);      histMixedMom_IM->GetXaxis()->SetBinLabel( binId,"");


        for (int i=5; i<13; i++ )
        {
            double binCont = histMixedMom_IM->GetBinContent( i+1 +1 );
            double err = histMixedMom_IM->GetBinError( i+1 +1 );

            histMixedMom_IM->SetBinContent( i+1, binCont);
            histMixedMom_IM->SetBinError( i+1, err);
            histMixedMom_IM->GetXaxis()->SetBinLabel( i+1, histMixedMom_IM->GetXaxis()->GetBinLabel(i+1+1) );
        }

    }



    histMixedMom_IM->DrawCopy();


//    line_unity = new TLine(-0.5,1,histMixedMom_IM->GetNbinsX()-0.5,1);
    line_unity = new TLine(-0.5+2,1,histMixedMom_IM->GetNbinsX()-0.5-4,1);
    line_unity->SetLineWidth(2);
    line_unity->SetLineStyle(9);
    line_unity->SetLineColor(kGray+1);
    line_unity->DrawClone();

    histMixedMom_IM->DrawCopy("same");


    // ### write
    TFile *file_results = new TFile( "file_results_canvases.root", "recreate" );

    canv_FirstMom->Write();
    canv_SecMom->Write();
    canv_MixedMom->Write();

    canv_FirstMom_RATIO->Write();
    canv_SecMom_RATIO->Write();
    canv_MixedMom_RATIO->Write();

    file_results->Close();
}














