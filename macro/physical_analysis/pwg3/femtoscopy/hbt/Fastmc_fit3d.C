

void Fastmc_fit3d()
{
// FASTMC, pions
//macros to fit 3D CFs (Fastmc_CF3d.C) and to compare with STAR data from
//J. Adams et al., STAR Collaboration,  Phys. Rev. C71, 044906 (2005).


 const int NPT = 3 ;

 Double_t rout [NPT] ;
 Double_t rside [NPT] ;
 Double_t rlong [NPT] ;
 Double_t lambda [NPT] ;
 Double_t erout [NPT] ;
 Double_t erside [NPT] ;
 Double_t erlong [NPT] ;
 Double_t elambda [NPT] ;

//star0.2,0.3,0.4,0.525

// intervals of transverse momentum pt
 Double_t pt_mean[NPT] = {0.05, 0.15, 0.3};

 Double_t ept[NPT] ={0.0005, 0.0005, 0.0005};

 TString binnames[NPT] = {"0-01","01-02","02-04"};

 //STAR central 0--5%:
//pi+
// Malinina
/*
 Double_t ptSTAR[3]={0.05,0.15,0.3};
 Double_t RoSTAR[3]={6.1,6.0,4.5};
 Double_t RsSTAR[3]={6.0,5.6,4.5};
 Double_t RlSTAR[3]={5.6,5.5,4.0};
 */
 // Dryablov (8AGeV) 
 Double_t ptSTAR[3]={0.05,0.15,0.3};
 Double_t RoSTAR[3]={5.236,6.206,4.609};
 Double_t RsSTAR[3]={5.41,5.03,4.136};
 Double_t RlSTAR[3]={5.17,5.3,3.247};
 Double_t lamSTAR[3]={0.513,0.722,0.535};
 
  Double_t eRoSTAR[3]={0.35,0.38,0.424};
 Double_t eRsSTAR[3]={0.35,0.22,0.304};
 Double_t eRlSTAR[3]={0.34,0.21,0.26};
 Double_t eSTARl[3]={0.05,0.06,0.07};
 // Musulmabekov(8AGeV):
 /*
 Double_t ptSTAR[3]={0.069,0.149,0.299};
 Double_t RoSTAR[3]={5.153,6.061,4.444};
 Double_t RsSTAR[3]={5.333,4.875,4.19};
 Double_t RlSTAR[3]={5.165,5.275,3.187};
 */
// Double_t lamSTAR[3]={0.7,0.9,0.9};
// Double_t eSTAR[3]={0.1,0.3,0.5};
// Double_t eSTARl[3]={0.03,0.05,0.1};
 
 Double_t eptSTAR[3]={0.01,0.01,0.01};

 Double_t ptPion[3]={0.2,0.425,0.75};
 Double_t eptPion[3]={0.001,0.001,0.001};
 Double_t RoPion[3]={6.37,6.0,6.36};
 Double_t eRoPion[3]={0.37,0.61,1.67};
 Double_t RsPion[3]={4.87,4.04,4.4};
 Double_t eRsPion[3]={0.25,0.44,1.57};
 Double_t RlPion[3]={8.55,6.60,10.1};
 Double_t eRlPion[3]={0.5,0.67,4.55};
 Double_t lamPion[3]={0.5,0.55,0.70};
 Double_t elamPion[3]={0.04,0.0755,0.32};


 for(int ipt=0; ipt<NPT; ipt++)   //  pt loop
 {

 TF3 *fitc = new TF3("fitc","1+[3]*exp(-25.76578*(x*x*[0]*[0]+y*y*[1]*[1]+z*z*[2]*[2]))");
 fitc->SetParameters(3.0, 3.0, 3.0, 1.0);

//Fastmc_CF3d.C Output file
 TFile *f = new TFile("CFs.root");

 TH3F *hCF = (TH3F*)f->Get("hCF"+binnames[ipt]);

 hCF->Fit(fitc,"Q0");
 fitc = (TF3*)(hCF->GetFunction("fitc"));
 rout[ipt] = fitc->GetParameter(0);
 rside[ipt] = fitc->GetParameter(1);
 rlong[ipt] = fitc->GetParameter(2);
 lambda[ipt] = fitc->GetParameter(3);
 erout[ipt] = fitc->GetParError(0);
 erside[ipt] = fitc->GetParError(1);
 erlong[ipt] = fitc->GetParError(2);
 elambda[ipt] = fitc->GetParError(3);

 printf("pt = %f  Ro = %f +/- %f  Rs = %f +/- %f Rl = %f +/- %f  lambda =  %f +/- %f\n", pt_mean[ipt],
        rout[ipt],erout[ipt], rside[ipt],erside[ipt], rlong[ipt],erlong[ipt], lambda[ipt],elambda[ipt]);

 } // end pt loop



 Int_t NPT1=NPT;


 TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
                                                                                                  
  gPad->SetFillStyle(4000);
  gPad->SetFillColor(0);
  c1->Divide(2,2,0.001,0.001);
                                                                                             
 
 //c1->GetFrame()->SetBorderMode(-1);
 c1->cd(1);
 gPad->SetTopMargin(0.02);
 gPad->SetRightMargin(0.1);
 gPad->SetLeftMargin(0.2);
 gPad->SetBottomMargin(0.2);
 TMultiGraph *mg =new TMultiGraph();

 TGraph *gr = new TGraphErrors(3,ptSTAR,RoSTAR,eptSTAR,eRoSTAR);
 gr->SetMarkerStyle(24);
 gr->SetLineColor(1); 


 gr->SetMarkerColor(1);
 
 gr->SetMarkerSize(0.8);
 mg->Add(gr);


 TGraph *gr1 = new TGraphErrors(2,ptPion,RoPion,eptPion,eRoPion);
 gr1->SetMarkerStyle(20);
 gr1->SetLineColor(1);
 gr1->SetMarkerColor(1);
 gr1->SetMarkerSize(0.8);

// mg->Add(gr1);

 TGraph *out = new TGraphErrors(3, pt_mean, rout,ept, erout) ;
 out->SetMarkerStyle(22);
 out->SetLineColor(1);
 out->SetMarkerColor(1);
 out->SetMarkerSize(0.8);

 mg->Add(out);
 mg->SetMaximum(10);
 mg->SetMinimum(1);
 mg->Draw("AP");
 mg->GetXaxis()->SetTitle("k_{t} (GeV/c)");
 mg->GetXaxis()->SetLabelSize(0.05);
 mg->GetXaxis()->SetTitleSize(0.06);
 mg->GetXaxis()->SetTitleOffset(0.8);
 mg->GetYaxis()->SetTitle("R_{out} (fm) ");
 mg->GetYaxis()->SetLabelSize(0.05);
 mg->GetYaxis()->SetTitleSize(0.06);
 mg->GetYaxis()->SetTitleOffset(0.7);

 c1->cd(2);
 gPad->SetTopMargin(0.02);
 gPad->SetRightMargin(0.1);
 gPad->SetLeftMargin(0.2);
 gPad->SetBottomMargin(0.2);
 TMultiGraph *mg =new TMultiGraph();
 TGraph *gr = new TGraphErrors(3,ptSTAR,RsSTAR,eptSTAR,eRsSTAR);
 gr->SetMarkerStyle(24);
 gr->SetLineColor(1);
 gr->SetMarkerColor(1);
 gr->SetMarkerSize(0.8);
// gr->SetMarkerSize(0.8);

 mg->Add(gr);


 TGraph *gr1 = new TGraphErrors(2,ptPion,RsPion,eptPion,eRsPion);
 gr1->SetMarkerStyle(20);
 gr1->SetLineColor(1);
 gr1->SetMarkerColor(1);
 gr1->SetMarkerSize(0.8);

// mg->Add(gr1);

 TGraph *side = new TGraphErrors(3, pt_mean, rside, ept, erside) ;
 side->SetMarkerStyle(22);
 side->SetLineColor(1);
 side->SetMarkerColor(1);
 side->SetMarkerSize(0.8);

 mg->Add(side);
 mg->SetMaximum(10);
 mg->SetMinimum(1);
 mg->Draw("AP");
 mg->GetXaxis()->SetTitle("k_{t} (GeV/c)");
 mg->GetXaxis()->SetLabelSize(0.05);
 mg->GetXaxis()->SetTitleSize(0.06);
 mg->GetXaxis()->SetTitleOffset(0.8);
 mg->GetYaxis()->SetTitle("R_{side} (fm) ");
 mg->GetYaxis()->SetLabelSize(0.05);
 mg->GetYaxis()->SetTitleSize(0.06);
 mg->GetYaxis()->SetTitleOffset(0.7);

 TLegend *legend=new TLegend(0.6, 0.6, 0.9, 0.9);
 legend->AddEntry(gr, "E895(8) experiment", "p");
// legend->AddEntry(gr1, "UrQMD, direct #pi", "p");
 legend->AddEntry(side, "FASTMC,  #pi^{+}", "p");
 legend->Draw();

 c1->cd(1);
 //TLatex *latex=new TLatex();
// latex->SetNDC();
 //latex->SetTextSize(0.07);
 //latex->DrawLatex(0.2, 0.7, "Au+Au #sqrt{s} = 200 AGeV, 6\% centrality");

 c1->cd(3);
 gPad->SetTopMargin(0.02);
 gPad->SetRightMargin(0.1);
 gPad->SetLeftMargin(0.2);
 gPad->SetBottomMargin(0.2);
 TMultiGraph *mg =new TMultiGraph();
 TGraph *gr = new TGraphErrors(3,ptSTAR,RlSTAR,eptSTAR,eRlSTAR);
 gr->SetMarkerStyle(24);
 gr->SetLineColor(1);
 gr->SetMarkerSize(0.8);

// gr->SetMarkerSize(2);
 gr->SetMarkerColor(1);
 mg->Add(gr);

 TGraph *gr1 = new TGraphErrors(2,ptPion,RlPion,eptPion,eRlPion);
 gr1->SetMarkerStyle(20);
 gr1->SetLineColor(1);
 gr1->SetMarkerColor(1);
 gr1->SetMarkerSize(0.8);

 //mg->Add(gr1);

 TGraph *llong = new TGraphErrors(3, pt_mean, rlong, ept, erlong) ;
 llong->SetMarkerStyle(22);
 llong->SetLineColor(1);
 llong->SetMarkerColor(1);
 llong->SetMarkerSize(0.8);

 mg->Add(llong);
 mg->SetMaximum(10);
 mg->SetMinimum(1);
 mg->Draw("AP");
 mg->GetXaxis()->SetTitle("k_{t} (GeV/c)");
 mg->GetXaxis()->SetLabelSize(0.05);
 mg->GetXaxis()->SetTitleSize(0.06);
 mg->GetXaxis()->SetTitleOffset(0.8);
 mg->GetYaxis()->SetTitle("R_{long} (fm) ");
 mg->GetYaxis()->SetLabelSize(0.05);
 mg->GetYaxis()->SetTitleSize(0.06);
 mg->GetYaxis()->SetTitleOffset(0.7);

  c1->cd(4);
 gPad->SetTopMargin(0.02);
 gPad->SetRightMargin(0.1);
 gPad->SetLeftMargin(0.2);
 gPad->SetBottomMargin(0.2);

 TMultiGraph *mg =new TMultiGraph();
 TGraph *gr = new TGraphErrors(3,ptSTAR,lamSTAR,eptSTAR,eSTARl);
 gr->SetMarkerStyle(24);
 gr->SetLineColor(1);
 gr->SetMarkerSize(0.8);
 gr->SetMarkerColor(1);
 mg->Add(gr);

 TGraph *gr1 = new TGraphErrors(2,ptPion,lamPion,eptPion,elamPion);
 gr1->SetMarkerStyle(20);
 gr1->SetLineColor(1);
 gr1->SetMarkerColor(1);
 gr1->SetMarkerSize(0.8);

// mg->Add(gr1);

 TGraph *glambda = new TGraphErrors(3, pt_mean, lambda, ept, elambda) ;
 glambda->SetMarkerStyle(22);
 glambda->SetLineColor(1);
 glambda->SetMarkerColor(1);
 glambda->SetMarkerSize(0.8);

 mg->Add(glambda);
 mg->SetMaximum(2);
 mg->SetMinimum(0);
 mg->Draw("AP");
 mg->GetXaxis()->SetTitle("k_{t} (GeV/c)");
 mg->GetXaxis()->SetLabelSize(0.05);
 mg->GetXaxis()->SetTitleSize(0.06);
 mg->GetXaxis()->SetTitleOffset(0.8);
 mg->GetYaxis()->SetTitle("#lambda");
 mg->GetYaxis()->SetLabelSize(0.05);
 mg->GetYaxis()->SetTitleSize(0.06);
 mg->GetYaxis()->SetTitleOffset(0.7);



}
