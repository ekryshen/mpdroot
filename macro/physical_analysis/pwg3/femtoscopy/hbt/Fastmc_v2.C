 #include "stdlib.h"
#include "stdio.h"

// the macros to create the distribution on the v2 
// (only FASTMC)
// and comparison with data from the experiment STAR.

void fig_v2()
{
  Double_t v2pi[100], v1pi[100], v2K[100], v2p[100], ev2pi[100], ev2K[100], ev2p[100],pt, ptt[100], ept[100];
  Int_t n2pi[100], n2K[100], n2p[100];
 
// exp. fata file (7-column figures)
   std::ifstream out1("v2_STAR/v2star20-30"); //20-30
   Double_t ptstar[100],ptlow[100],pthight[100],v2star[100],ev2star[100],epts[100],x1,x2;
   Int_t i=0;
   
   while(out1)
   {
   out1>>ptstar[i]>>ptlow[i]>>pthight[i]>>v2star[i]>>ev2star[i]>>x1>>x2;
   epts[i]= 0.01;//(-ptlow[i]+pthight[i])/ptstar[i];   
        
   printf("%f %f %3f\n",ptstar[i], v2star[i],ev2star[i]);
   i++;       
   }
   out1.close();
 
   TH1D *hv2 = new TH1D("hv2", "hv2", 100, 0.01, 2.);
   TH1D *hv0 = new TH1D("hv0", "hv0", 100, 0.01, 2.);
   TH1D *hv2res1 = new TH1D("hv2res1", "hv2res1", 100, 0.01, 2.);
   
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
 
  
  for (Int_t k=0;k<nevents;k++) 
  {
            
      td->GetEntry(k);  
              
     for (Int_t i=0;i<npart;i++) 
	   {

      if(((pdg[i]==211)||(abs(pdg[i])==321)||(abs(pdg[i])==2212)) 
      && (abs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))<1.0)){

      pt = TMath::Sqrt(Px[i]*Px[i]+Py[i]*Py[i]);
      
      phi = TMath::Abs(TMath::ATan(Py[i]/Px[i]));

      if(Px[i]<0&&Py[i]>0)phi = TMath::Pi()-phi;
      if(Px[i]<0&&Py[i]<0)phi = TMath::Pi()+phi;
      if(Px[i]>0&&Py[i]<0)phi = 2.*TMath::Pi()-phi;
                                                                          
      v2 = TMath::Cos(2*phi);  
      double ev2=0.1;
  // printf("%f %f %3f\n",pt, v2,ev2);

      hv2->Fill(pt,v2);
      hv0->Fill(pt,1.);
      
    } //if    
  } // for
} // for


    TCanvas *c1= new TCanvas("c1OQ","efficiencie and contaminations  ",200,10,800,600);
    

          TGraph *gr5= new TGraphErrors(10,ptstar,v2star,epts,ev2star);
          gr5->SetMarkerColor(4);
          gr5->SetLineColor(4);
          gr5->SetMarkerStyle(20);

	  TMultiGraph *mg=new TMultiGraph();
          mg->SetMaximum(0.2);
          mg->SetMinimum(0.);
          mg->Add(gr5);          
          mg->Draw("AP");
          mg->GetXaxis()->SetTitle("pt, GeV/c");
          mg->GetYaxis()->SetTitle(" v2 ");

         hv2res1->Divide(hv2,hv0,1.,1.);
         hv2res1->SetLineWidth(1);
         hv2res1->Draw("histo:same");

 }
