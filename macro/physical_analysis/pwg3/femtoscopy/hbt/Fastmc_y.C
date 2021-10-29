#include "stdlib.h"
#include "stdio.h"

// the macros to create the distribution on the y 
// (only FASTMC, pions)
// and comparison with data from the experiment E895.

void Fastmc_y()
{
 
//  exp. fata file (3-column figures)
   std::ifstream out1("mt_E895/8gev_pim_dndy"); //20-30
   Double_t mtstarpi[100],dNdmtpi[100],edNdmtpi[100],edNdmtpi1[100], epts[100];
   Int_t i=0;
   
   while(out1)
   {
   out1>>mtstarpi[i]>>dNdmtpi[i]>>edNdmtpi[i];
   epts[i]= 0.01;//(-ptlow[i]+pthight[i])/ptstar[i];   
        
  // printf("%f %f %3f\n",ptstar[i], v2star[i],ev2star[i]);
   i++;       
   }
   out1.close();
     
   Int_t npart = 0;
   Float_t v2, phi;
   
// FASTMC Output file  
   TFile *f = new TFile("RunOutput.root");
   const Int_t kMax = 20000; 

   Int_t   pdg[kMax];
   Int_t   Mpdg[kMax];
   Float_t Px[kMax];
   Float_t Py[kMax];
   Float_t Pz[kMax];
   Float_t E[kMax];   
   Float_t X[kMax];
   Float_t Y[kMax];
   Float_t Z[kMax];
   Float_t T[kMax]; 

  TTree *td = (TTree*)f->Get("td");
  Int_t nevents = td->GetEntries();  
  Info("mult.C", "Nevents %d ", nevents);  
  
  td->SetBranchAddress("npart",&npart);
  td->SetBranchAddress("Px",Px);
  td->SetBranchAddress("Py",Py);
  td->SetBranchAddress("Pz",Pz);
  td->SetBranchAddress("E",E);
  td->SetBranchAddress("X",X);
  td->SetBranchAddress("Y",Y);
  td->SetBranchAddress("Z",Z);
  td->SetBranchAddress("T",T);
  td->SetBranchAddress("pdg",pdg);
  td->SetBranchAddress("Mpdg",Mpdg);
 
   TH1D *hy = new TH1D("hy", "hy", 100, -1.2, 1.2);
   
 
   td->Draw("0.5*log((E+Pz)/(E-Pz))>>hy", "pdg==-211");

   hy->Scale(100.0/ (nevents*2.4));
   
     
   gStyle->SetOptStat(10000000);
   gStyle->SetStatBorderSize(0);

   TCanvas *c2 = new TCanvas("c2", "c2",364,44,800,600);
   gStyle->SetOptStat(0);
   c2->Range(-24.9362,1.71228,25.0213,4.77534);
   c2->SetFillColor(10);
   c2->SetFillStyle(4000);
   c2->SetBorderSize(2);
   c2->SetFrameFillColor(0);
   c2->SetFrameFillStyle(0);

    
          TGraph *gr1= new TGraphErrors(25,mtstarpi,dNdmtpi,epts,edNdmtpi); 
          gr1->SetMarkerColor(1);
          gr1->SetLineColor(1);
          gr1->SetMarkerStyle(20);
          
	  TMultiGraph *mg=new TMultiGraph();
         
          mg->Draw("AP");
          mg->Add(gr1);     
          mg->GetXaxis()->SetTitleSize(0.05);
          mg->GetYaxis()->SetTitleSize(0.06); 
          mg->GetXaxis()->SetTitleOffset(1.);
	  mg->GetYaxis()->SetTitleOffset(0.7); 
          mg->GetXaxis()->SetTitle("y = 0.5 ln((E+P_{z})/(E-P_{z}))");
          mg->GetYaxis()->SetTitle("dN/dy");

          TLegend *legend=new TLegend(0.6, 0.6, 0.9, 0.9);
          legend->AddEntry(gr1, "E895(8) 0-5%  #pi^{-}", "p");
          legend->AddEntry(hy, "FASTMC", "l");
          legend->Draw();
         
          hy->Draw("same::hist");

 }
