#include "stdlib.h"
#include "stdio.h"

// the macros to create the distribution on the transverse mass mt 
// (only FASTMC)
// and comparison with data from the experiment NA49.

void Fastmc_mt()
{
 
// exp. fata file (4-column figures)
   std::ifstream out1("mt_E895/4gev_pim_y0"); //central 4gev_pim_y0 b1: 4gev_pim_b1 
   Double_t mtstarpi[100],dNdmtpi[100],edNdmtpi[100],edNdmtpi1[100], epts[100];
   Int_t i=0;
   
   while(out1)
   {
   out1>>mtstarpi[i]>>dNdmtpi[i]>>edNdmtpi[i]>>edNdmtpi1[i];
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
 
   TH1D *hpt1 = new TH1D("hpt1", "hpt1", 100, 0., 1.);
   TH1D *hpt2 = new TH1D("hpt2", "hpt2", 100, 0., 1.);
   TH1D *hpt3 = new TH1D("hpt3", "hpt3", 100, 0., 1.);
   
 
  td->Draw("sqrt(Px*Px+Py*Py+0.13957*0.13957)-0.13957>>hpt1", "(1.0/(sqrt(Px*Px+Py*Py+0.13957*0.13957)))*(pdg==-211&&0.5*log((E+Pz)/(E-Pz))>=-0.05&&0.5*log((E+Pz)/(E-Pz))<=0.05)");
 //  td->Draw("sqrt(Px*Px+Py*Py+0.13957*0.13957)-0.13957>>hpt1", "(1.0/(sqrt(Px*Px+Py*Py+0.13957*0.13957)))*(pdg==-211&&0.5*log((E+Pz)/(E-Pz))>=-0.15&&0.5*log((E+Pz)/(E-Pz))<=-0.05)");

//b2 td->Draw("sqrt(Px*Px+Py*Py+0.13957*0.13957)-0.13957>>hpt1","(1.0/(sqrt(Px*Px+Py*Py+0.13957*0.13957)))
//*(pdg==-211&&log((E+Pz)/(E-Pz))>=0.15&&log((E+Pz)/(E-Pz))<=0.25)");
//test td->Draw("sqrt(Px*Px+Py*Py+0.13957*0.13957)-0.13957>>hpt1","(1.0/(sqrt(Px*Px+Py*Py+0.13957*0.13957)))*
//(pdg==-211&&log((E+Pz)/(E-Pz))>=-0.25&&Mpdg==-1&&log((E+Pz)/(E-Pz))<=-0.15)");

   hpt1->Scale(100.0/ (2.0*TMath::Pi()*nevents*1.*0.1)); // dp=1.8/18; dy=[-1:1]=2
   
        //  hptnw1->SetLineStyle(2);
       //   hptnw2->SetLineStyle(2);
        //  hptnw3->SetLineStyle(2);
    
   gStyle->SetOptStat(10000000);
   gStyle->SetStatBorderSize(0);

   TCanvas *c2 = new TCanvas("c2", "c2",364,44,800,600);
   gStyle->SetOptStat(0);
   c2->Range(-24.9362,1.71228,25.0213,4.77534);
   c2->SetFillColor(10);
   c2->SetFillStyle(4000);
   c2->SetBorderSize(2);
   c2->SetLogy();
   c2->SetFrameFillColor(0);
   c2->SetFrameFillStyle(0);

    
 //   c1->Divide(1,3);
                                                                                              
          TGraph *gr1= new TGraphErrors(40,mtstarpi,dNdmtpi,epts,edNdmtpi); //27 b1 2AGeV
          gr1->SetMarkerColor(1);
          gr1->SetLineColor(1);
          gr1->SetMarkerStyle(20);
          
           

	  TMultiGraph *mg=new TMultiGraph();
	 //mg->SetMaximum(0.2);
    //     mg->SetMinimum(0.);
         
          mg->Draw("AP");
          mg->Add(gr1);     
//          mg->Add(gr2);
//          mg->Add(gr3);
          mg->GetXaxis()->SetTitleSize(0.05);
          mg->GetYaxis()->SetTitleSize(0.05); 
          mg->GetXaxis()->SetTitleOffset(0.9);
	       mg->GetYaxis()->SetTitleOffset(0.9); 
          mg->GetXaxis()->SetTitle("m_{t} (GeV/c)");
          mg->GetYaxis()->SetTitle("1/(2 #pi) d^{2}N/ N m_{t} dm_{t} dY, c^{2}/GeV^{2}");

          TLegend *legend=new TLegend(0.6, 0.6, 0.9, 0.9);
          legend->AddEntry(gr1, "E895 0-5%  #pi^{-}", "p");
//          legend->AddEntry(gr2, "             K ", "p");
//          legend->AddEntry(gr3, "             p", "p");
          legend->AddEntry(hpt1, "FASTMC", "l");
	  // legend->AddEntry(hpt2, " K", "l");
          // legend->AddEntry(hpt3, "p", "l");
           // legend->AddEntry(hv2res3_165, "FASTMC: T= 165 MeV", "l");
            legend->Draw();


         
          hpt1->Draw("same::hist");

 }
