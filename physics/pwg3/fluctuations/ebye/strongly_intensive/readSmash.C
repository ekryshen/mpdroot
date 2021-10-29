// E. Andronov 2020
// this code can be used to extract fluctuation observables from the SMASH output (ROOT version)

#include "TTree.h"
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "TMath.h"
#include "TString.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>

using namespace::std;

void readSmash(char * filelist_chunk = " ",char * outname = " "){
  cout << "start of script" << endl;
  ifstream base(filelist_chunk);

  //ifstream base("/lhep/users/evandron/myFilelist_urqmd_11gev_fromParfenov_full.txt");
  //ifstream base("myFilelist_urqmd_7point7gev_fromParfenov_full.txt");
  //ifstream base("myFilelist_vhlle_7point7gev_fromParfenov_full.txt");

  char *str = new char [1000];

  TChain chain("particles");
  int num_of_files_in_folder = 0;
  while (!base.eof())
  {
    base.getline(str, 1000, '\n');
    cout << str << endl;
    if ( TString(str).Length() < 3 ) // last line reached (typically it's empty)
    break;

    num_of_files_in_folder++;

    chain.Add( TString(str) );

  }

  cout << "num_of_files_in_folder = " << num_of_files_in_folder << endl;
  Int_t nevents = 0;

  Int_t ev = chain.GetEntries();

  if (nevents == 0) nevents = ev;
  cout << " Number of events in DST file = " << nevents << endl;

  const int MAX_N_PART = 10000;
  Int_t npart;
  Double_t impact_b;
  Double_t px[MAX_N_PART];
  Double_t py[MAX_N_PART];
  Double_t pz[MAX_N_PART];
  Int_t charge[MAX_N_PART];

  TH2D* pTMultHist;
  pTMultHist = new TH2D("pTMultHist","pTMultHist;n;mean p_{T};entries",20,-0.5,19.5,500,-0.5,2.5);

  chain.SetBranchAddress("npart",&npart);
  chain.SetBranchAddress("impact_b",&impact_b);
  chain.SetBranchAddress("px",px);
  chain.SetBranchAddress("py",py);
  chain.SetBranchAddress("pz",pz);
  chain.SetBranchAddress("charge",charge);

  long double fullSumPt = 0;
  long double fullEBEPt = 0;
  long double fullEBEPt_2 = 0;
  long double fullEBEPt_3 = 0;

  long double fullEBEdptdpt = 0;
  long double fullEBEdptdpt_2 = 0;
  long double fullEBEdptdptdpt = 0;

  long double fullSumPt2 = 0;
  long double fullSumPtpow2 = 0;
  long double fullSumPt3 = 0;
  long double fullSumPtpow3 = 0;
  long double fullSumMult = 0;
  long double fullSumMultMult = 0;
  long double fullSumMultMultMult = 0;
  long double fullSumMultPt = 0;
  long double fullSumMultMultPt = 0;
  long double fullSumPt2Pt = 0;
  long double fullSumPt2Mult = 0;
  long double fullSumPtpow2Mult = 0;

  long double meanSumPt = 0;
  long double meanEBEPt = 0;
  long double meanEBEPt_2 = 0;
  long double meanEBEPt_3 = 0;
  long double meanEBESumPt = 0;

  long double meanEBEdptdpt = 0;
  long double meanEBEdptdpt_2 = 0;
  long double meanEBEdptdptdpt = 0;

  long double meanSumPt2 = 0;
  long double meanSumPtcomma2 = 0;
  long double meanSumPt3 = 0;
  long double meanSumPtcomma3 = 0;
  long double meanSumPtcomma2Pt = 0;
  long double meanMultPt2 = 0;
  long double meanMultPtcomma2 = 0;



  long double meanPt = 0;
  long double meanPt2 = 0;
  long double meanMultPt = 0;
  long double meanMultMultPt = 0;
  long double meanX3 = 0;
  long double meanX2 = 0;
  long double meanX = 0;
  long double meanEBEX2 = 0;
  long double meanEBEX = 0;
  long double varX = 0;
  long double varX3 = 0;
  long double ptOmega = 0;
  long double sumPtOmega = 0;

  long double nEventsIn;
  nEventsIn = 0;
  cout << "nEventsIn = " << nEventsIn << endl;

  long double nEventsInSubsample_D_nonzeroTotalMult;
  nEventsInSubsample_D_nonzeroTotalMult=0;

  long double nEventsInSubsample_D_nononeTotalMult;
  nEventsInSubsample_D_nononeTotalMult=0;

  long double nEventsInSubsample_D_nontwoTotalMult;
  nEventsInSubsample_D_nontwoTotalMult=0;

  Int_t totalMult;
  Double_t totalPt;
  Double_t totalPt2;
  Double_t totalPt3;

  for (Int_t i = 0; i < nevents; i++){
    //cout << "event number = " << i << endl;


    chain.GetEntry(i);
    if(i%10000==0){cout << "event number = " << i << endl; cout << "impact_b = " << impact_b << endl;}

    if (impact_b>1.54) {
      continue;
    }
    //cout << "npart = " << npart << endl;
    totalMult = 0;
    totalPt = 0;
    totalPt2 = 0;
    totalPt3 = 0;

    double p, pt, eta;

    for (Int_t i = 0; i < npart; i++){
      //    cout << "px[" << i << "] = " << px[i] << endl;
      p = TMath::Sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
      pt = TMath::Sqrt(px[i]*px[i] + py[i]*py[i]);
      eta = 0.5 * TMath::Log( ( p + pz[i] )/( p - pz[i] ) );
      if (charge[i]!=0 && eta<1 && eta>-1 && pt>0.15 && pt<2) {
        totalMult = totalMult + 1;
        totalPt = totalPt + pt;
        totalPt2 = totalPt2 + pt*pt;
        totalPt3 = totalPt3 + pt*pt*pt;
      }
    }

    nEventsIn = nEventsIn+1;
    if (totalMult>0) {
      nEventsInSubsample_D_nonzeroTotalMult = nEventsInSubsample_D_nonzeroTotalMult + 1;
    }
    if (totalMult>1) {
      nEventsInSubsample_D_nononeTotalMult = nEventsInSubsample_D_nononeTotalMult + 1;
    }
    if (totalMult>2) {
      nEventsInSubsample_D_nontwoTotalMult = nEventsInSubsample_D_nontwoTotalMult + 1;
    }



    fullSumPt += totalPt;
    if(totalMult>0){
      fullEBEPt += totalPt/totalMult;
      pTMultHist->Fill(totalMult, totalPt/totalMult);
    }
    if(totalMult>1){
      fullEBEPt_2 += totalPt/totalMult;
      fullEBEdptdpt += (totalPt*totalPt-totalPt2)/(totalMult*totalMult-totalMult);
    }
    if(totalMult>2){
      fullEBEPt_3 += totalPt/totalMult;
      fullEBEdptdpt_2 += (totalPt*totalPt-totalPt2)/(totalMult*totalMult-totalMult);
      fullEBEdptdptdpt += (totalPt*totalPt*totalPt - 3*totalPt*totalPt2 + 2*totalPt3)/(totalMult*totalMult*totalMult-3*totalMult*totalMult + 2*totalMult);
    }

    fullSumPt2 += totalPt2;
    fullSumPt3 += totalPt3;
    fullSumPtpow2 += totalPt*totalPt;
    fullSumPtpow3 += totalPt*totalPt*totalPt;
    fullSumMult += totalMult;
    fullSumMultMult += totalMult*totalMult;
    fullSumMultMultMult += totalMult*totalMult*totalMult;
    fullSumMultPt += totalMult*totalPt;
    fullSumPt2Pt += totalPt*totalPt2;
    fullSumPt2Mult += totalPt2*totalMult;
    fullSumPtpow2Mult += totalPt*totalPt*totalMult;
    fullSumMultMultPt += totalMult*totalMult*totalPt;

  }
  //TCanvas* c = new TCanvas("c","c",800,800);
  //c->cd();
  //pTMultHist->Draw("colz");
  //c->SaveAs("ptn.pdf");

  meanSumPt = fullSumPt/nEventsIn;

  meanEBEPt = fullEBEPt / nEventsInSubsample_D_nonzeroTotalMult;
  meanEBEPt_2 = fullEBEPt_2 / nEventsInSubsample_D_nononeTotalMult;
  meanEBEPt_3 = fullEBEPt_3 / nEventsInSubsample_D_nontwoTotalMult;

  meanEBESumPt = fullSumPt / nEventsInSubsample_D_nonzeroTotalMult;
  meanEBEdptdpt = fullEBEdptdpt / nEventsInSubsample_D_nononeTotalMult;
  meanEBEdptdpt_2 = fullEBEdptdpt_2 / nEventsInSubsample_D_nontwoTotalMult;
  meanEBEdptdptdpt = fullEBEdptdptdpt / nEventsInSubsample_D_nontwoTotalMult;

  meanSumPt2 = fullSumPtpow2 / nEventsIn;
  meanPt = fullSumPt / fullSumMult;
  meanPt2 = fullSumPt2 / fullSumMult;
  meanSumPtcomma2 = fullSumPt2 / nEventsIn;
  meanMultPt = fullSumMultPt / nEventsIn;
  meanMultMultPt = fullSumMultMultPt / nEventsIn;

  meanSumPt3 = fullSumPtpow3 / nEventsIn;
  meanSumPtcomma3 = fullSumPt3 / nEventsIn;
  meanSumPtcomma2Pt = fullSumPt2Pt / nEventsIn;
  meanMultPt2 = fullSumPtpow2Mult / nEventsIn;
  meanMultPtcomma2 = fullSumPt2Mult / nEventsIn;

  meanX3 = fullSumMultMultMult / nEventsIn;
  meanX2 = fullSumMultMult / nEventsIn;
  meanX = fullSumMult / nEventsIn;
  meanEBEX2 = fullSumMultMult / nEventsInSubsample_D_nonzeroTotalMult;
  meanEBEX = fullSumMult / nEventsInSubsample_D_nonzeroTotalMult;
  long double omegaEBE=0;
  omegaEBE = (meanEBEX2 - meanEBEX*meanEBEX)/meanEBEX;

  varX = meanX2 - meanX*meanX;
  varX3 = meanX3 - 3*meanX2*meanX2 + 2*meanX2;

  ptOmega = (meanPt2 - meanPt*meanPt) / meanPt;
  //cout << "ptOmega " << ptOmega << endl;
  sumPtOmega = (meanSumPt2 - meanSumPt*meanSumPt) / meanSumPt;


  double mult, meanpt_star, meanpt_alice, dptdpt, dptdpt_star, dptdpt_star_2, gammapt, gammapt_star, delta, sigma, corr;

  mult = meanX;

  meanpt_star = meanEBEPt;

  meanpt_alice = meanPt;
  dptdpt = ((meanSumPt2-meanSumPtcomma2) - 2*meanPt*(meanMultPt-meanSumPt) + meanPt*meanPt*(meanX2 - meanX))/(meanX2 - meanX);

  dptdpt_star = meanEBEdptdpt - meanEBEPt_2*meanEBEPt_2;

  dptdpt_star_2 = meanEBEdptdpt_2 - meanEBEPt_3*meanEBEPt_3;
  gammapt_star = meanEBEdptdptdpt - 3*meanEBEdptdpt_2*meanEBEPt_3 + 2*meanEBEPt_3*meanEBEPt_3*meanEBEPt_3;
  gammapt_star = meanEBEPt_3*gammapt_star/(dptdpt_star_2*dptdpt_star_2);

  delta = (meanX*sumPtOmega - meanSumPt*varX/meanX) / (meanX*ptOmega);
  sigma = (meanX*sumPtOmega + meanSumPt*varX/meanX - 2*(meanMultPt-meanX*meanSumPt)) / (meanX*ptOmega);
  long double tmp_1, tmp_2, tmp_3;
  tmp_1 = (meanSumPt3 - 3*meanSumPtcomma2Pt +2*meanSumPtcomma3)/varX3;
  tmp_2 = meanPt*(meanMultPt2 - 2*meanSumPt2 - meanMultPtcomma2 + 2*meanSumPtcomma2)/varX3;
  tmp_3 = meanPt*meanPt*(meanMultMultPt - 3*meanMultPt +2*meanSumPt)/varX3;
  gammapt = tmp_1 - 3*tmp_2 + 3*tmp_3 - meanPt*meanPt*meanPt;
  gammapt = meanPt*gammapt/(dptdpt*dptdpt);


  corr = ((meanEBESumPt/meanEBEPt)-meanEBEX)/omegaEBE;




  std::ofstream out;
  out.open(outname);
  if (out.is_open()) {
    out << "mult = " << mult << endl;
    out << "meanpt_star = " << meanpt_star << endl;
    out << "meanpt_alice = " << meanpt_alice << endl;
    out << "dptdpt_star = " << dptdpt_star << endl;
    out << "dptdpt = " << dptdpt << endl;
    out << "gammapt_star = " << gammapt_star << endl;
    out << "gammapt = " << gammapt << endl;
    out << "delta = " << delta << endl;
    out << "sigma = " << sigma << endl;
    out << "corr = " << corr << endl;
  }


}
