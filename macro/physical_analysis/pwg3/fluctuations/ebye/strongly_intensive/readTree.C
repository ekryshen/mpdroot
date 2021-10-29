// E.Andronov 2020
// this code allows to analyze tree extracted from the miniDST files
// in order to get numerical values of the fluctuation observables
// such as Sigma[N,PT], Delta[N,PT] etc.

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TCanvas.h"
#include <iostream>

static const Int_t NMaxTrack = 12000;
using namespace std;



void analyzeTree(TTree* extractedTree, Int_t nSubSamples, Int_t nEventsInSubSample);

void analyzeTree(TTree* extractedTree, Int_t nSubSamples, Int_t nEventsInSubSample){

  char strLowCentrality[35];
  int lowCentrality = 0;
  sprintf(strLowCentrality, "%d", lowCentrality);
  TString strrLowCentrality(strLowCentrality);

  char strHighCentrality[35];
  int highCentrality = 10;
  sprintf(strHighCentrality, "%d", highCentrality);
  TString strrHighCentrality(strHighCentrality);

  char strNSub[35];
  sprintf(strNSub, "%d", nSubSamples);
  TString strrNSub(strNSub);


  Double_t nSubSamplesDbl;
  nSubSamplesDbl = (Double_t)nSubSamples;

  Double_t dtmp2 = sqrt(nSubSamplesDbl*(nSubSamplesDbl-1)); //for subsample

  Int_t nOneWindows;
  nOneWindows = 10;

  Int_t totalMult[nOneWindows];
  Double_t totalPt[nOneWindows];
  Double_t totalPt2[nOneWindows];
  Double_t totalPt3[nOneWindows];

  Int_t totalMultPrimary[nOneWindows];
  Double_t totalPtPrimary[nOneWindows];
  Double_t totalPt2Primary[nOneWindows];
  Double_t totalPt3Primary[nOneWindows];

  Int_t totalMultSim[nOneWindows];
  Double_t totalPtSim[nOneWindows];
  Double_t totalPt2Sim[nOneWindows];
  Double_t totalPt3Sim[nOneWindows];


  Float_t impactPar;
  Int_t f_SIM_NumberOfTracks; //!
  Float_t f_SIM_TrackPt[NMaxTrack],f_SIM_TrackEta[NMaxTrack];
  Int_t f_SIM_TrackCharge[NMaxTrack];
  Int_t f_SIM_fromGen[NMaxTrack];
  Int_t fNumberOfTracks; //!
  Float_t fTrackPt[NMaxTrack],fTrackEta[NMaxTrack];
  Float_t fTrackPtPrimary[NMaxTrack],fTrackEtaPrimary[NMaxTrack];
  Int_t fTrackCharge[NMaxTrack];
  Float_t fDCA_x[NMaxTrack];
  Float_t fDCA_y[NMaxTrack];
  Float_t fDCA_z[NMaxTrack];
  Int_t fHasPrimaryTrack[NMaxTrack];
  extractedTree->SetBranchAddress("fImpactPar",&impactPar);
  extractedTree->SetBranchAddress("SIM_NumberOfTracks",&f_SIM_NumberOfTracks);
  extractedTree->SetBranchAddress("SIM_TrackPt",f_SIM_TrackPt);
  extractedTree->SetBranchAddress("SIM_TrackEta",f_SIM_TrackEta);
  extractedTree->SetBranchAddress("SIM_TrackCharge",f_SIM_TrackCharge);
  extractedTree->SetBranchAddress("SIM_fromGen",f_SIM_fromGen);
  extractedTree->SetBranchAddress("nTracks",&fNumberOfTracks);
  extractedTree->SetBranchAddress("trackPt",fTrackPt);
  extractedTree->SetBranchAddress("trackEta",fTrackEta);
  extractedTree->SetBranchAddress("trackPtPrimary",fTrackPtPrimary);
  extractedTree->SetBranchAddress("trackEtaPrimary",fTrackEtaPrimary);
  extractedTree->SetBranchAddress("trackCharge",fTrackCharge);
  extractedTree->SetBranchAddress("DCA_x",fDCA_x);
  extractedTree->SetBranchAddress("DCA_y",fDCA_y);
  extractedTree->SetBranchAddress("DCA_z",fDCA_z);
  extractedTree->SetBranchAddress("hasPrimaryTrack",fHasPrimaryTrack);

  TH1D* pTHist[16][nSubSamples];
  TH1D* meanpTHist[16][nSubSamples];
  TH2D* pTMultHist[16][nSubSamples];

  for (int i=0; i<10; i++) {
    char str[15];
    sprintf(str, "%d", i);
    TString strr(str);

    for (Int_t h=0; h<nSubSamples; h++) {
      char strSub[15];
      sprintf(strSub, "%d", h);
      TString strrSub(strSub);

      pTHist[i][h] = new TH1D(strrNSub+"x"+strr+"x"+strrSub+"pTHist",strrNSub+"x"+strr+"x"+strrSub+"pTHist;p_{T};entries",500,-0.5,2.5);
      meanpTHist[i][h] = new TH1D(strrNSub+"x"+strr+"x"+strrSub+"meanpTHist",strrNSub+"x"+strr+"x"+strrSub+"meanpTHist;mean p_{T};entries",500,-0.5,2.5);
      pTMultHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"pTMultHist",strrNSub+"x"+strr+"x"+strrSub+"pTMultHist;n;mean p_{T};entries",500,-0.5,499.5,500,-0.5,2.5);
    }
  }

  Int_t nStartingEvent, nEndingEvent;
  nStartingEvent = 0;
  nEndingEvent = nStartingEvent+nEventsInSubSample;

  long double mult_sim_sub[10][nSubSamples];
  long double mult_global_sub[10][nSubSamples];
  long double mult_primary_sub[10][nSubSamples];
  long double mult_sim_mean[10] = {0};
  long double mult_global_mean[10] = {0};
  long double mult_primary_mean[10] = {0};
  long double mult_sim_error[10] = {0};
  long double mult_global_error[10] = {0};
  long double mult_primary_error[10] = {0};


  long double meanpt_star_sim_sub[10][nSubSamples];
  long double meanpt_star_global_sub[10][nSubSamples];
  long double meanpt_star_primary_sub[10][nSubSamples];
  long double meanpt_star_sim_mean[10] = {0};
  long double meanpt_star_global_mean[10] = {0};
  long double meanpt_star_primary_mean[10] = {0};
  long double meanpt_star_sim_error[10] = {0};
  long double meanpt_star_global_error[10] = {0};
  long double meanpt_star_primary_error[10] = {0};

  long double meanpt_alice_sim_sub[10][nSubSamples];
  long double meanpt_alice_global_sub[10][nSubSamples];
  long double meanpt_alice_primary_sub[10][nSubSamples];
  long double meanpt_alice_sim_mean[10] = {0};
  long double meanpt_alice_global_mean[10] = {0};
  long double meanpt_alice_primary_mean[10] = {0};
  long double meanpt_alice_sim_error[10] = {0};
  long double meanpt_alice_global_error[10] = {0};
  long double meanpt_alice_primary_error[10] = {0};

  long double delta_sim_sub[10][nSubSamples];
  long double delta_global_sub[10][nSubSamples];
  long double delta_primary_sub[10][nSubSamples];
  long double delta_sim_mean[10] = {0};
  long double delta_global_mean[10] = {0};
  long double delta_primary_mean[10] = {0};
  long double delta_sim_error[10] = {0};
  long double delta_global_error[10] = {0};
  long double delta_primary_error[10] = {0};

  long double sigma_sim_sub[10][nSubSamples];
  long double sigma_global_sub[10][nSubSamples];
  long double sigma_primary_sub[10][nSubSamples];
  long double sigma_sim_mean[10] = {0};
  long double sigma_global_mean[10] = {0};
  long double sigma_primary_mean[10] = {0};
  long double sigma_sim_error[10] = {0};
  long double sigma_global_error[10] = {0};
  long double sigma_primary_error[10] = {0};

  long double corr_sim_sub[10][nSubSamples];
  long double corr_global_sub[10][nSubSamples];
  long double corr_primary_sub[10][nSubSamples];
  long double corr_sim_mean[10] = {0};
  long double corr_global_mean[10] = {0};
  long double corr_primary_mean[10] = {0};
  long double corr_sim_error[10] = {0};
  long double corr_global_error[10] = {0};
  long double corr_primary_error[10] = {0};

  long double dptdpt_sim_sub[10][nSubSamples];
  long double dptdpt_global_sub[10][nSubSamples];
  long double dptdpt_primary_sub[10][nSubSamples];
  long double dptdpt_sim_mean[10] = {0};
  long double dptdpt_global_mean[10] = {0};
  long double dptdpt_primary_mean[10] = {0};
  long double dptdpt_sim_error[10] = {0};
  long double dptdpt_global_error[10] = {0};
  long double dptdpt_primary_error[10] = {0};


  long double dptdpt_star_sim_sub[10][nSubSamples];
  long double dptdpt_star_global_sub[10][nSubSamples];
  long double dptdpt_star_primary_sub[10][nSubSamples];
  long double dptdpt_star_sim_mean[10] = {0};
  long double dptdpt_star_global_mean[10] = {0};
  long double dptdpt_star_primary_mean[10] = {0};
  long double dptdpt_star_sim_error[10] = {0};
  long double dptdpt_star_global_error[10] = {0};
  long double dptdpt_star_primary_error[10] = {0};

  long double gammapt_sim_sub[10][nSubSamples];
  long double gammapt_global_sub[10][nSubSamples];
  long double gammapt_primary_sub[10][nSubSamples];
  long double gammapt_sim_mean[10] = {0};
  long double gammapt_global_mean[10] = {0};
  long double gammapt_primary_mean[10] = {0};
  long double gammapt_sim_error[10] = {0};
  long double gammapt_global_error[10] = {0};
  long double gammapt_primary_error[10] = {0};

  long double gammapt_star_sim_sub[10][nSubSamples];
  long double gammapt_star_global_sub[10][nSubSamples];
  long double gammapt_star_primary_sub[10][nSubSamples];
  long double gammapt_star_sim_mean[10] = {0};
  long double gammapt_star_global_mean[10] = {0};
  long double gammapt_star_primary_mean[10] = {0};
  long double gammapt_star_sim_error[10] = {0};
  long double gammapt_star_global_error[10] = {0};
  long double gammapt_star_primary_error[10] = {0};


  for (int h=0; h<nSubSamples; h++) {
    cout << "subsample #" << h << endl;
    long double nEventsInSubSample_D;
    nEventsInSubSample_D=0;

    long double nEventsInSubsample_D_nonzeroTotalMultSim;
    nEventsInSubsample_D_nonzeroTotalMultSim=0;

    long double nEventsInSubsample_D_nonzeroTotalMult;
    nEventsInSubsample_D_nonzeroTotalMult=0;

    long double nEventsInSubsample_D_nonzeroTotalMultPrimary;
    nEventsInSubsample_D_nonzeroTotalMultPrimary=0;

    long double fullSumPt[10] = {0};
    long double fullEBEPt[10] = {0};

    long double fullEBEdptdpt[10] = {0};
    long double fullEBEdptdptdpt[10] = {0};

    long double fullSumPt2[10] = {0};
    long double fullSumPtpow2[10] = {0};
    long double fullSumPt3[10] = {0};
    long double fullSumPtpow3[10] = {0};
    long double fullSumMult[10] = {0};
    long double fullSumMultMult[10] = {0};
    long double fullSumMultMultMult[10] = {0};
    long double fullSumMultPt[10] = {0};
    long double fullSumMultMultPt[10] = {0};
    long double fullSumPt2Pt[10] = {0};
    long double fullSumPt2Mult[10] = {0};
    long double fullSumPtpow2Mult[10] = {0};

    long double meanSumPt[10]={0};
    long double meanEBEPt[10] = {0};
    long double meanEBESumPt[10] = {0};

    long double meanEBEdptdpt[10] = {0};
    long double meanEBEdptdptdpt[10] = {0};

    long double meanSumPt2[10]={0};
    long double meanSumPtcomma2[10]={0};
    long double meanSumPt3[10]={0};
    long double meanSumPtcomma3[10]={0};
    long double meanSumPtcomma2Pt[10]={0};
    long double meanMultPt2[10]={0};
    long double meanMultPtcomma2[10]={0};



    long double meanPt[10]={0};
    long double meanPt2[10]={0};
    long double meanMultPt[10]={0};
    long double meanMultMultPt[10]={0};
    long double meanX3[10]={0};
    long double meanX2[10]={0};
    long double meanX[10]={0};
    long double meanEBEX2[10]={0};
    long double meanEBEX[10]={0};
    long double varX[10]={0};
    long double varX3[10]={0};
    long double ptOmega[10]={0};
    long double sumPtOmega[10]={0};




    long double fullSumPtPrimary[10] = {0};
    long double fullEBEPtPrimary[10] = {0};

    long double fullEBEdptdptPrimary[10] = {0};
    long double fullEBEdptdptdptPrimary[10] = {0};

    long double fullSumPt2Primary[10] = {0};
    long double fullSumPtpow2Primary[10] = {0};
    long double fullSumPt3Primary[10] = {0};
    long double fullSumPtpow3Primary[10] = {0};
    long double fullSumMultPrimary[10] = {0};
    long double fullSumMultMultPrimary[10] = {0};
    long double fullSumMultMultMultPrimary[10] = {0};
    long double fullSumMultPtPrimary[10] = {0};
    long double fullSumMultMultPtPrimary[10] = {0};
    long double fullSumPt2PtPrimary[10] = {0};
    long double fullSumPt2MultPrimary[10] = {0};
    long double fullSumPtpow2MultPrimary[10] = {0};

    long double meanSumPtPrimary[10]={0};
    long double meanEBEPtPrimary[10] = {0};
    long double meanEBESumPtPrimary[10] = {0};

    long double meanEBEdptdptPrimary[10] = {0};
    long double meanEBEdptdptdptPrimary[10] = {0};

    long double meanSumPt2Primary[10]={0};
    long double meanSumPtcomma2Primary[10]={0};
    long double meanSumPt3Primary[10]={0};
    long double meanSumPtcomma3Primary[10]={0};
    long double meanSumPtcomma2PtPrimary[10]={0};
    long double meanMultPt2Primary[10]={0};
    long double meanMultPtcomma2Primary[10]={0};

    long double meanPtPrimary[10]={0};
    long double meanPt2Primary[10]={0};
    long double meanMultPtPrimary[10]={0};
    long double meanMultMultPtPrimary[10]={0};
    long double meanX3Primary[10]={0};
    long double meanX2Primary[10]={0};
    long double meanXPrimary[10]={0};
    long double meanEBEX2Primary[10]={0};
    long double meanEBEXPrimary[10]={0};
    long double varXPrimary[10]={0};
    long double varX3Primary[10]={0};
    long double ptOmegaPrimary[10]={0};
    long double sumPtOmegaPrimary[10]={0};



    long double fullSumPtSim[10] = {0};
    long double fullEBEPtSim[10] = {0};

    long double fullEBEdptdptSim[10] = {0};
    long double fullEBEdptdptdptSim[10] = {0};

    long double fullSumPt2Sim[10] = {0};
    long double fullSumPtpow2Sim[10] = {0};
    long double fullSumPt3Sim[10] = {0};
    long double fullSumPtpow3Sim[10] = {0};
    long double fullSumMultSim[10] = {0};
    long double fullSumMultMultSim[10] = {0};
    long double fullSumMultMultMultSim[10] = {0};
    long double fullSumMultPtSim[10] = {0};
    long double fullSumMultMultPtSim[10] = {0};
    long double fullSumPt2PtSim[10] = {0};
    long double fullSumPt2MultSim[10] = {0};
    long double fullSumPtpow2MultSim[10] = {0};

    long double meanSumPtSim[10]={0};
    long double meanEBEPtSim[10] = {0};
    long double meanEBESumPtSim[10] = {0};

    long double meanEBEdptdptSim[10] = {0};
    long double meanEBEdptdptdptSim[10] = {0};

    long double meanSumPt2Sim[10]={0};
    long double meanSumPt3Sim[10]={0};
    long double meanSumPtcomma2Sim[10]={0};
    long double meanSumPtcomma3Sim[10]={0};
    long double meanSumPtcomma2PtSim[10]={0};
    long double meanMultPt2Sim[10]={0};
    long double meanMultPtcomma2Sim[10]={0};

    long double meanPtSim[10]={0};
    long double meanPt2Sim[10]={0};
    long double meanMultPtSim[10]={0};
    long double meanMultMultPtSim[10]={0};
    long double meanX3Sim[10]={0};
    long double meanX2Sim[10]={0};
    long double meanXSim[10]={0};
    long double meanEBEX2Sim[10]={0};
    long double meanEBEXSim[10]={0};
    long double varXSim[10]={0};
    long double varX3Sim[10]={0};
    long double ptOmegaSim[10]={0};
    long double sumPtOmegaSim[10]={0};


    for (int j=nStartingEvent; j<nEndingEvent; j++) {
      extractedTree->GetEntry(j);

      for (Int_t jW=0; jW<nOneWindows; jW++) {
        totalMult[jW] = 0;
        totalPt[jW] = 0;
        totalPt2[jW] = 0;
        totalPt3[jW] = 0;

        totalMultSim[jW] = 0;
        totalPtSim[jW] = 0;
        totalPt2Sim[jW] = 0;
        totalPt3Sim[jW] = 0;

        totalMultPrimary[jW] = 0;
        totalPtPrimary[jW] = 0;
        totalPt2Primary[jW] = 0;
        totalPt3Primary[jW] = 0;
      }

      //implement event cuts below!!!

      double bForTreeLow, bForTreeHigh;
      bForTreeLow = -1;
      bForTreeHigh = 100;


      switch (lowCentrality) {
        case 0:
        bForTreeLow = -1;
        break;
        case 5:
        bForTreeLow = 3.23611;
        break;
        case 10:
        bForTreeLow = 4.59567;
        break;
        case 15:
        bForTreeLow = 5.63829;
        break;
        case 20:
        bForTreeLow = 6.51336;
        break;
        case 25:
        bForTreeLow = 7.28364;
        break;
        case 30:
        bForTreeLow = 7.97662;
        break;
        case 35:
        bForTreeLow = 8.61874;
        break;
        case 40:
        bForTreeLow = 9.21147;
        break;
        case 45:
        bForTreeLow = 9.76649;
        break;
        case 50:
        bForTreeLow = 10.2912;
        break;
        case 55:
        bForTreeLow = 10.7946;
        break;
        case 60:
        bForTreeLow = 11.2812;
        break;
        case 65:
        bForTreeLow = 11.742;
        break;
        case 70:
        bForTreeLow = 12.1956;
        break;
        case 75:
        bForTreeLow = 12.624;
        break;
        case 80:
        bForTreeLow = 13.0411;
        break;
        case 85:
        bForTreeLow = 13.466;
        break;
        case 90:
        bForTreeLow = 13.9349;
        break;
        case 95:
        bForTreeLow = 14.5656;
        break;
        default:
        bForTreeLow = -1;
        break;
      }




      switch (highCentrality) {
        case 0:
        cout << "don't do that, idiot";
        break;
        case 5:
        bForTreeHigh = 3.23611;
        break;
        case 10:
        bForTreeHigh = 4.59567;
        break;
        case 15:
        bForTreeHigh = 5.63829;
        break;
        case 20:
        bForTreeHigh = 6.51336;
        break;
        case 25:
        bForTreeHigh = 7.28364;
        break;
        case 30:
        bForTreeHigh = 7.97662;
        break;
        case 35:
        bForTreeHigh = 8.61874;
        break;
        case 40:
        bForTreeHigh = 9.21147;
        break;
        case 45:
        bForTreeHigh = 9.76649;
        break;
        case 50:
        bForTreeHigh = 10.2912;
        break;
        case 55:
        bForTreeHigh = 10.7946;
        break;
        case 60:
        bForTreeHigh = 11.2812;
        break;
        case 65:
        bForTreeHigh = 11.742;
        break;
        case 70:
        bForTreeHigh = 12.1956;
        break;
        case 75:
        bForTreeHigh = 12.624;
        break;
        case 80:
        bForTreeHigh = 13.0411;
        break;
        case 85:
        bForTreeHigh = 13.466;
        break;
        case 90:
        bForTreeHigh = 13.9349;
        break;
        case 95:
        bForTreeHigh = 14.5656;
        break;
        default:
        bForTreeHigh = 100;
        break;
      }


      if (impactPar>bForTreeHigh || impactPar<bForTreeLow) {
        continue;
      }


      for (int i=0; i<f_SIM_NumberOfTracks; i++) {
        for (Int_t jW=0; jW<nOneWindows; jW++) {
          if (f_SIM_TrackEta[i]<0.1+jW*0.1 && f_SIM_TrackEta[i]>-0.1-jW*0.1 && f_SIM_fromGen[i]==1) {
            totalMultSim[jW] = totalMultSim[jW] + 1;
            totalPtSim[jW] = totalPtSim[jW] + f_SIM_TrackPt[i];
            totalPt2Sim[jW] = totalPt2Sim[jW] + f_SIM_TrackPt[i]*f_SIM_TrackPt[i];
            totalPt3Sim[jW] = totalPt3Sim[jW] + f_SIM_TrackPt[i]*f_SIM_TrackPt[i]*f_SIM_TrackPt[i];
          }
        }
      }


      for (int i=0; i<fNumberOfTracks; i++) {
        for (Int_t jW=0; jW<nOneWindows; jW++) {
          if (fTrackEta[i]<0.1+jW*0.1 && fTrackEta[i]>-0.1-jW*0.1 && fTrackPt[i]>0.15 && fTrackPt[i]<2.0) {
            if (fDCA_x[i]*fDCA_x[i]+fDCA_y[i]*fDCA_y[i]+fDCA_z[i]*fDCA_z[i]<5) {
              totalMult[jW] = totalMult[jW] + 1;
              totalPt[jW] = totalPt[jW] + fTrackPt[i];
              totalPt2[jW] = totalPt2[jW] + fTrackPt[i]*fTrackPt[i];
              totalPt3[jW] = totalPt3[jW] + fTrackPt[i]*fTrackPt[i]*fTrackPt[i];
            }
          }
          if (fHasPrimaryTrack[i] == 1) {
            if (fTrackEtaPrimary[i]<0.1+jW*0.1 && fTrackEtaPrimary[i]>-0.1-jW*0.1 && fTrackPtPrimary[i]>0.15 && fTrackPtPrimary[i]<2.0) {
              if (fDCA_x[i]*fDCA_x[i]+fDCA_y[i]*fDCA_y[i]+fDCA_z[i]*fDCA_z[i]<5) {
                totalMultPrimary[jW] = totalMultPrimary[jW] + 1;
                totalPtPrimary[jW] = totalPtPrimary[jW] + fTrackPtPrimary[i];
                totalPt2Primary[jW] = totalPt2Primary[jW] + fTrackPtPrimary[i]*fTrackPtPrimary[i];
                totalPt3Primary[jW] = totalPt3Primary[jW] + fTrackPtPrimary[i]*fTrackPtPrimary[i]*fTrackPtPrimary[i];
                pTHist[jW][h]->Fill(fTrackPtPrimary[i]);
              }
            }
          }
        }
      }

      nEventsInSubSample_D++;

      if(totalMult[0]>0){
        nEventsInSubsample_D_nonzeroTotalMult++;
      }

      if(totalMultSim[0]>0){
        nEventsInSubsample_D_nonzeroTotalMultSim++;
      }

      if(totalMultPrimary[0]>0){
        nEventsInSubsample_D_nonzeroTotalMultPrimary++;
      }

      for (int i=0; i<10; i++) {
        fullSumPtSim[i] += totalPtSim[i];
        if(totalMultSim[i]>0){
          fullEBEPtSim[i] += totalPtSim[i]/totalMultSim[i];
        }
        fullEBEdptdptSim[i] += (totalPtSim[i]*totalPtSim[i]-totalPt2Sim[i])/(totalMultSim[i]*totalMultSim[i]-totalMultSim[i]);
        fullEBEdptdptdptSim[i] += (totalPtSim[i]*totalPtSim[i]*totalPtSim[i] - 3*totalPtSim[i]*totalPt2Sim[i] + 2*totalPt3Sim[i])/(totalMultSim[i]*totalMultSim[i]*totalMultSim[i]-3*totalMultSim[i]*totalMultSim[i] + 2*totalMultSim[i]);

        fullSumPt2Sim[i] += totalPt2Sim[i];
        fullSumPt3Sim[i] += totalPt3Sim[i];
        fullSumPtpow2Sim[i] += totalPtSim[i]*totalPtSim[i];
        fullSumPtpow3Sim[i] += totalPtSim[i]*totalPtSim[i]*totalPtSim[i];
        fullSumMultSim[i] += totalMultSim[i];
        fullSumMultMultSim[i] += totalMultSim[i]*totalMultSim[i];
        fullSumMultMultMultSim[i] += totalMultSim[i]*totalMultSim[i]*totalMultSim[i];
        fullSumMultPtSim[i] += totalMultSim[i]*totalPtSim[i];
        fullSumPt2PtSim[i] += totalPtSim[i]*totalPt2Sim[i];
        fullSumPt2MultSim[i] += totalPt2Sim[i]*totalMultSim[i];
        fullSumPtpow2MultSim[i] += totalPtSim[i]*totalPtSim[i]*totalMultSim[i];
        fullSumMultMultPtSim[i] += totalMultSim[i]*totalMultSim[i]*totalPtSim[i];

        /////////////////////////////////////////////////
        fullSumPt[i] += totalPt[i];
        if(totalMult[i]>0){
          fullEBEPt[i] += totalPt[i]/totalMult[i];
        }
        fullEBEdptdpt[i] += (totalPt[i]*totalPt[i]-totalPt2[i])/(totalMult[i]*totalMult[i]-totalMult[i]);
        fullEBEdptdptdpt[i] += (totalPt[i]*totalPt[i]*totalPt[i] - 3*totalPt[i]*totalPt2[i] + 2*totalPt3[i])/(totalMult[i]*totalMult[i]*totalMult[i]-3*totalMult[i]*totalMult[i] + 2*totalMult[i]);

        fullSumPt2[i] += totalPt2[i];
        fullSumPt3[i] += totalPt3[i];
        fullSumPtpow2[i] += totalPt[i]*totalPt[i];
        fullSumPtpow3[i] += totalPt[i]*totalPt[i]*totalPt[i];
        fullSumMult[i] += totalMult[i];
        fullSumMultMult[i] += totalMult[i]*totalMult[i];
        fullSumMultMultMult[i] += totalMult[i]*totalMult[i]*totalMult[i];
        fullSumMultPt[i] += totalMult[i]*totalPt[i];
        fullSumPt2Pt[i] += totalPt[i]*totalPt2[i];
        fullSumPt2Mult[i] += totalPt2[i]*totalMult[i];
        fullSumPtpow2Mult[i] += totalPt[i]*totalPt[i]*totalMult[i];
        fullSumMultMultPt[i] += totalMult[i]*totalMult[i]*totalPt[i];


        /////////////////////////////////////////////////
        fullSumPtPrimary[i] += totalPtPrimary[i];
        if(totalMultPrimary[i]>0){
          fullEBEPtPrimary[i] += totalPtPrimary[i]/totalMultPrimary[i];
          meanpTHist[i][h]->Fill(totalPtPrimary[i]/totalMultPrimary[i]);
          pTMultHist[i][h]->Fill(totalMultPrimary[i],totalPtPrimary[i]/totalMultPrimary[i]);
        }
        fullEBEdptdptPrimary[i] += (totalPtPrimary[i]*totalPtPrimary[i]-totalPt2Primary[i])/(totalMultPrimary[i]*totalMultPrimary[i]-totalMultPrimary[i]);
        fullEBEdptdptdptPrimary[i] += (totalPtPrimary[i]*totalPtPrimary[i]*totalPtPrimary[i] - 3*totalPtPrimary[i]*totalPt2Primary[i] + 2*totalPt3Primary[i])/(totalMultPrimary[i]*totalMultPrimary[i]*totalMultPrimary[i]-3*totalMultPrimary[i]*totalMultPrimary[i] + 2*totalMultPrimary[i]);


        fullSumPt2Primary[i] += totalPt2Primary[i];
        fullSumPt3Primary[i] += totalPt3Primary[i];
        fullSumPtpow2Primary[i] += totalPtPrimary[i]*totalPtPrimary[i];
        fullSumPtpow3Primary[i] += totalPtPrimary[i]*totalPtPrimary[i]*totalPtPrimary[i];
        fullSumMultPrimary[i] += totalMultPrimary[i];
        fullSumMultMultPrimary[i] += totalMultPrimary[i]*totalMultPrimary[i];
        fullSumMultMultMultPrimary[i] += totalMultPrimary[i]*totalMultPrimary[i]*totalMultPrimary[i];
        fullSumMultPtPrimary[i] += totalMultPrimary[i]*totalPtPrimary[i];
        fullSumPt2PtPrimary[i] += totalPtPrimary[i]*totalPt2Primary[i];
        fullSumPt2MultPrimary[i] += totalPt2Primary[i]*totalMultPrimary[i];
        fullSumPtpow2MultPrimary[i] += totalPtPrimary[i]*totalPtPrimary[i]*totalMultPrimary[i];
        fullSumMultMultPtPrimary[i] += totalMultPrimary[i]*totalMultPrimary[i]*totalPtPrimary[i];
      }


    };
    nStartingEvent = nEndingEvent;
    nEndingEvent = nEndingEvent + nEventsInSubSample;

    for (int i=0; i<10; i++) {
      meanSumPtSim[i] = fullSumPtSim[i] / nEventsInSubSample_D;
      meanEBEPtSim[i] = fullEBEPtSim[i] / nEventsInSubsample_D_nonzeroTotalMultSim;
      meanEBESumPtSim[i] = fullSumPtSim[i] / nEventsInSubsample_D_nonzeroTotalMultSim;
      meanEBEdptdptSim[i] = fullEBEdptdptSim[i] / nEventsInSubsample_D_nonzeroTotalMultSim;
      meanEBEdptdptdptSim[i] = fullEBEdptdptdptSim[i] / nEventsInSubsample_D_nonzeroTotalMultSim;

      meanSumPt2Sim[i] = fullSumPtpow2Sim[i] / nEventsInSubSample_D;
      meanPtSim[i] = fullSumPtSim[i] / fullSumMultSim[i];
      meanPt2Sim[i] = fullSumPt2Sim[i] / fullSumMultSim[i];
      meanSumPtcomma2Sim[i] = fullSumPt2Sim[i] / nEventsInSubSample_D;
      meanMultPtSim[i] = fullSumMultPtSim[i] / nEventsInSubSample_D;
      meanMultMultPtSim[i] = fullSumMultMultPtSim[i] / nEventsInSubSample_D;

      meanSumPt3Sim[i] = fullSumPtpow3Sim[i] / nEventsInSubSample_D;
      meanSumPtcomma3Sim[i] = fullSumPt3Sim[i] / nEventsInSubSample_D;
      meanSumPtcomma2PtSim[i] = fullSumPt2PtSim[i] / nEventsInSubSample_D;
      meanMultPt2Sim[i] = fullSumPtpow2MultSim[i] / nEventsInSubSample_D;
      meanMultPtcomma2Sim[i] = fullSumPt2MultSim[i] / nEventsInSubSample_D;

      meanX3Sim[i] = fullSumMultMultMultSim[i] / nEventsInSubSample_D;
      meanX2Sim[i] = fullSumMultMultSim[i] / nEventsInSubSample_D;
      meanXSim[i] = fullSumMultSim[i] / nEventsInSubSample_D;
      meanEBEX2Sim[i] = fullSumMultMultSim[i] / nEventsInSubsample_D_nonzeroTotalMultSim;
      meanEBEXSim[i] = fullSumMultSim[i] / nEventsInSubsample_D_nonzeroTotalMultSim;
      long double omegaEBE=0;
      omegaEBE = (meanEBEX2Sim[i] - meanEBEXSim[i]*meanEBEXSim[i])/meanEBEXSim[i];

      varXSim[i] = meanX2Sim[i] - meanXSim[i]*meanXSim[i];
      varX3Sim[i] = meanX3Sim[i] - 3*meanX2Sim[i]*meanX2Sim[i] + 2*meanX2Sim[i];

      ptOmegaSim[i] = (meanPt2Sim[i] - meanPtSim[i]*meanPtSim[i]) / meanPtSim[i];
      //cout << "ptOmega " << ptOmega << endl;
      sumPtOmegaSim[i] = (meanSumPt2Sim[i] - meanSumPtSim[i]*meanSumPtSim[i]) / meanSumPtSim[i];


      mult_sim_sub[i][h] = meanXSim[i];

      meanpt_star_sim_sub[i][h] = meanEBEPtSim[i];
      meanpt_alice_sim_sub[i][h] = meanPtSim[i];
      dptdpt_sim_sub[i][h] = ((meanSumPt2Sim[i]-meanSumPtcomma2Sim[i]) - 2*meanPtSim[i]*(meanMultPtSim[i]-meanSumPtSim[i]) + meanPtSim[i]*meanPtSim[i]*(meanX2Sim[i] - meanXSim[i]))/(meanX2Sim[i] - meanXSim[i]);
      dptdpt_star_sim_sub[i][h] = meanEBEdptdptSim[i] - meanEBEPtSim[i]*meanEBEPtSim[i];
      gammapt_star_sim_sub[i][h] = meanEBEdptdptdptSim[i] - 3*meanEBEdptdptSim[i]*meanEBEPtSim[i] + 2*meanEBEPtSim[i]*meanEBEPtSim[i]*meanEBEPtSim[i];
      gammapt_star_sim_sub[i][h] = meanEBEPtSim[i]*gammapt_star_sim_sub[i][h]/(dptdpt_star_sim_sub[i][h]*dptdpt_star_sim_sub[i][h]);
      //totalMultSub[i][h] = meanX[i];
      //totalOmegaSub[i][h] = varX[i]/meanX[i];
      //totalDeltaSub[i][h] = (meanX[i]*sumPtOmega[i] - meanSumPt[i]*totalOmegaSub[i][h]) / (meanX[i]*ptOmega[i]);
      //totalSigmaSub[i][h] = (meanX[i]*sumPtOmega[i] + meanSumPt[i]*totalOmegaSub[i][h] - 2*(meanMultPt[i]-meanX[i]*meanSumPt[i])) / (meanX[i]*ptOmega[i]);
      delta_sim_sub[i][h] = (meanXSim[i]*sumPtOmegaSim[i] - meanSumPtSim[i]*varXSim[i]/meanXSim[i]) / (meanXSim[i]*ptOmegaSim[i]);
      sigma_sim_sub[i][h] = (meanXSim[i]*sumPtOmegaSim[i] + meanSumPtSim[i]*varXSim[i]/meanXSim[i] - 2*(meanMultPtSim[i]-meanXSim[i]*meanSumPtSim[i])) / (meanXSim[i]*ptOmegaSim[i]);

      long double tmp_1, tmp_2, tmp_3;
      tmp_1 = (meanSumPt3Sim[i] - 3*meanSumPtcomma2PtSim[i] +2*meanSumPtcomma3Sim[i])/varX3Sim[i];
      tmp_2 = meanPtSim[i]*(meanMultPt2Sim[i] - 2*meanSumPt2Sim[i] - meanMultPtcomma2Sim[i] + 2*meanSumPtcomma2Sim[i])/varX3Sim[i];
      tmp_3 = meanPtSim[i]*meanPtSim[i]*(meanMultMultPtSim[i] - 3*meanMultPtSim[i] +2*meanSumPtSim[i])/varX3Sim[i];
      gammapt_sim_sub[i][h] = tmp_1 - 3*tmp_2 + 3*tmp_3 - meanPtSim[i]*meanPtSim[i]*meanPtSim[i];
      gammapt_sim_sub[i][h] = meanPtSim[i]*gammapt_sim_sub[i][h]/(dptdpt_sim_sub[i][h]*dptdpt_sim_sub[i][h]);
      //cout << i << " ; " << h << "gammapt_sim_sub[i][h] = " << gammapt_sim_sub[i][h] << endl;

      corr_sim_sub[i][h] = ((meanEBESumPtSim[i]/meanEBEPtSim[i])-meanEBEXSim[i])/omegaEBE;
      //cout << i << " ; " << h << "corr_sim_sub[i][h] = " << corr_sim_sub[i][h] << endl;



      /////////////////////////////////////////////////



      meanSumPt[i] = fullSumPt[i] / nEventsInSubSample_D;
      meanEBEPt[i] = fullEBEPt[i] / nEventsInSubsample_D_nonzeroTotalMult;
      meanEBESumPt[i] = fullSumPt[i] / nEventsInSubsample_D_nonzeroTotalMult;
      meanEBEdptdpt[i] = fullEBEdptdpt[i] / nEventsInSubsample_D_nonzeroTotalMult;
      meanEBEdptdptdpt[i] = fullEBEdptdptdpt[i] / nEventsInSubsample_D_nonzeroTotalMult;

      meanSumPt2[i] = fullSumPtpow2[i] / nEventsInSubSample_D;
      meanPt[i] = fullSumPt[i] / fullSumMult[i];
      meanPt2[i] = fullSumPt2[i] / fullSumMult[i];
      meanSumPtcomma2[i] = fullSumPt2[i] / nEventsInSubSample_D;
      meanMultPt[i] = fullSumMultPt[i] / nEventsInSubSample_D;
      meanMultMultPt[i] = fullSumMultMultPt[i] / nEventsInSubSample_D;

      meanSumPt3[i] = fullSumPtpow3[i] / nEventsInSubSample_D;
      meanSumPtcomma3[i] = fullSumPt3[i] / nEventsInSubSample_D;
      meanSumPtcomma2Pt[i] = fullSumPt2Pt[i] / nEventsInSubSample_D;
      meanMultPt2[i] = fullSumPtpow2Mult[i] / nEventsInSubSample_D;
      meanMultPtcomma2[i] = fullSumPt2Mult[i] / nEventsInSubSample_D;

      meanX3[i] = fullSumMultMultMult[i] / nEventsInSubSample_D;
      meanX2[i] = fullSumMultMult[i] / nEventsInSubSample_D;
      meanX[i] = fullSumMult[i] / nEventsInSubSample_D;
      meanEBEX2[i] = fullSumMultMult[i] / nEventsInSubsample_D_nonzeroTotalMult;
      meanEBEX[i] = fullSumMult[i] / nEventsInSubsample_D_nonzeroTotalMult;

      omegaEBE = (meanEBEX2[i] - meanEBEX[i]*meanEBEX[i])/meanEBEX[i];

      varX[i] = meanX2[i] - meanX[i]*meanX[i];
      varX3[i] = meanX3[i] - 3*meanX2[i]*meanX2[i] + 2*meanX2[i];

      ptOmega[i] = (meanPt2[i] - meanPt[i]*meanPt[i]) / meanPt[i];
      //cout << "ptOmega " << ptOmega << endl;
      sumPtOmega[i] = (meanSumPt2[i] - meanSumPt[i]*meanSumPt[i]) / meanSumPt[i];

      mult_global_sub[i][h] = meanX[i];

      meanpt_star_global_sub[i][h] = meanEBEPt[i];
      //cout << meanEBEPt[i] << endl;
      meanpt_alice_global_sub[i][h] = meanPt[i];
      dptdpt_global_sub[i][h] = ((meanSumPt2[i]-meanSumPtcomma2[i]) - 2*meanPt[i]*(meanMultPt[i]-meanSumPt[i]) + meanPt[i]*meanPt[i]*(meanX2[i] - meanX[i]))/(meanX2[i] - meanX[i]);
      dptdpt_star_global_sub[i][h] = meanEBEdptdpt[i] - meanEBEPt[i]*meanEBEPt[i];
      gammapt_star_global_sub[i][h] = meanEBEdptdptdpt[i] - 3*meanEBEdptdpt[i]*meanEBEPt[i] + 2*meanEBEPt[i]*meanEBEPt[i]*meanEBEPt[i];
      gammapt_star_global_sub[i][h] = meanEBEPt[i]*gammapt_star_global_sub[i][h]/(dptdpt_star_global_sub[i][h]*dptdpt_star_global_sub[i][h]);
      //totalMultSub[i][h] = meanX[i];
      //totalOmegaSub[i][h] = varX[i]/meanX[i];
      delta_global_sub[i][h] = (meanX[i]*sumPtOmega[i] - meanSumPt[i]*varX[i]/meanX[i]) / (meanX[i]*ptOmega[i]);
      sigma_global_sub[i][h] = (meanX[i]*sumPtOmega[i] + meanSumPt[i]*varX[i]/meanX[i] - 2*(meanMultPt[i]-meanX[i]*meanSumPt[i])) / (meanX[i]*ptOmega[i]);
      //totalDeltaSub[i][h] = (meanX[i]*sumPtOmega[i] - meanSumPt[i]*totalOmegaSub[i][h]) / (meanX[i]*ptOmega[i]);
      //totalSigmaSub[i][h] = (meanX[i]*sumPtOmega[i] + meanSumPt[i]*totalOmegaSub[i][h] - 2*(meanMultPt[i]-meanX[i]*meanSumPt[i])) / (meanX[i]*ptOmega[i]);
      tmp_1 = (meanSumPt3[i] - 3*meanSumPtcomma2Pt[i] +2*meanSumPtcomma3[i])/varX3[i];
      tmp_2 = meanPt[i]*(meanMultPt2[i] - 2*meanSumPt2[i] - meanMultPtcomma2[i] + 2*meanSumPtcomma2[i])/varX3[i];
      tmp_3 = meanPt[i]*meanPt[i]*(meanMultMultPt[i] - 3*meanMultPt[i] +2*meanSumPt[i])/varX3[i];
      gammapt_global_sub[i][h] = tmp_1 - 3*tmp_2 + 3*tmp_3 - meanPt[i]*meanPt[i]*meanPt[i];
      gammapt_global_sub[i][h] = meanPt[i]*gammapt_global_sub[i][h]/(dptdpt_global_sub[i][h]*dptdpt_global_sub[i][h]);


      corr_global_sub[i][h] = ((meanEBESumPt[i]/meanEBEPt[i])-meanEBEX[i])/omegaEBE;

      /////////////////////////////////////////////////

      meanSumPtPrimary[i] = fullSumPtPrimary[i] / nEventsInSubSample_D;
      meanEBEPtPrimary[i] = fullEBEPtPrimary[i] / nEventsInSubsample_D_nonzeroTotalMultPrimary;
      meanEBESumPtPrimary[i] = fullSumPtPrimary[i] / nEventsInSubsample_D_nonzeroTotalMultPrimary;
      meanEBEdptdptPrimary[i] = fullEBEdptdptPrimary[i] / nEventsInSubsample_D_nonzeroTotalMultPrimary;
      meanEBEdptdptdptPrimary[i] = fullEBEdptdptdptPrimary[i] / nEventsInSubsample_D_nonzeroTotalMultPrimary;


      meanSumPt2Primary[i] = fullSumPtpow2Primary[i] / nEventsInSubSample_D;
      meanPtPrimary[i] = fullSumPtPrimary[i] / fullSumMultPrimary[i];
      meanPt2Primary[i] = fullSumPt2Primary[i] / fullSumMultPrimary[i];
      meanSumPtcomma2Primary[i] = fullSumPt2Primary[i] / nEventsInSubSample_D;
      meanMultPtPrimary[i] = fullSumMultPtPrimary[i] / nEventsInSubSample_D;
      meanMultMultPtPrimary[i] = fullSumMultMultPtPrimary[i] / nEventsInSubSample_D;

      meanSumPt3Primary[i] = fullSumPtpow3Primary[i] / nEventsInSubSample_D;
      meanSumPtcomma3Primary[i] = fullSumPt3Primary[i] / nEventsInSubSample_D;
      meanSumPtcomma2PtPrimary[i] = fullSumPt2PtPrimary[i] / nEventsInSubSample_D;
      meanMultPt2Primary[i] = fullSumPtpow2MultPrimary[i] / nEventsInSubSample_D;
      meanMultPtcomma2Primary[i] = fullSumPt2MultPrimary[i] / nEventsInSubSample_D;

      meanX3Primary[i] = fullSumMultMultMultPrimary[i] / nEventsInSubSample_D;
      meanX2Primary[i] = fullSumMultMultPrimary[i] / nEventsInSubSample_D;
      meanXPrimary[i] = fullSumMultPrimary[i] / nEventsInSubSample_D;
      meanEBEX2Primary[i] = fullSumMultMultPrimary[i] / nEventsInSubsample_D_nonzeroTotalMultPrimary;
      meanEBEXPrimary[i] = fullSumMultPrimary[i] / nEventsInSubsample_D_nonzeroTotalMultPrimary;

      omegaEBE = (meanEBEX2Primary[i] - meanEBEXPrimary[i]*meanEBEXPrimary[i])/meanEBEXPrimary[i];

      varXPrimary[i] = meanX2Primary[i] - meanXPrimary[i]*meanXPrimary[i];
      varX3Primary[i] = meanX3Primary[i] - 3*meanX2Primary[i]*meanX2Primary[i] + 2*meanX2Primary[i];

      ptOmegaPrimary[i] = (meanPt2Primary[i] - meanPtPrimary[i]*meanPtPrimary[i]) / meanPtPrimary[i];
      //cout << "ptOmega " << ptOmega << endl;
      sumPtOmegaPrimary[i] = (meanSumPt2Primary[i] - meanSumPtPrimary[i]*meanSumPtPrimary[i]) / meanSumPtPrimary[i];


      mult_primary_sub[i][h] = meanXPrimary[i];

      meanpt_star_primary_sub[i][h] = meanEBEPtPrimary[i];
      meanpt_alice_primary_sub[i][h] = meanPtPrimary[i];
      dptdpt_primary_sub[i][h] = ((meanSumPt2Primary[i]-meanSumPtcomma2Primary[i]) - 2*meanPtPrimary[i]*(meanMultPtPrimary[i]-meanSumPtPrimary[i]) + meanPtPrimary[i]*meanPtPrimary[i]*(meanX2Primary[i] - meanXPrimary[i]))/(meanX2Primary[i] - meanXPrimary[i]);
      dptdpt_star_primary_sub[i][h] = meanEBEdptdptPrimary[i] - meanEBEPtPrimary[i]*meanEBEPtPrimary[i];
      gammapt_star_primary_sub[i][h] = meanEBEdptdptdptPrimary[i] - 3*meanEBEdptdptPrimary[i]*meanEBEPtPrimary[i] + 2*meanEBEPtPrimary[i]*meanEBEPtPrimary[i]*meanEBEPtPrimary[i];
      gammapt_star_primary_sub[i][h] = meanEBEPtPrimary[i]*gammapt_star_primary_sub[i][h]/(dptdpt_star_primary_sub[i][h]*dptdpt_star_primary_sub[i][h]);
      //totalMultSub[i][h] = meanX[i];
      //totalOmegaSub[i][h] = varX[i]/meanX[i];
      //totalDeltaSub[i][h] = (meanX[i]*sumPtOmega[i] - meanSumPt[i]*totalOmegaSub[i][h]) / (meanX[i]*ptOmega[i]);
      //totalSigmaSub[i][h] = (meanX[i]*sumPtOmega[i] + meanSumPt[i]*totalOmegaSub[i][h] - 2*(meanMultPt[i]-meanX[i]*meanSumPt[i])) / (meanX[i]*ptOmega[i]);
      delta_primary_sub[i][h] = (meanXPrimary[i]*sumPtOmegaPrimary[i] - meanSumPtPrimary[i]*varXPrimary[i]/meanXPrimary[i]) / (meanXPrimary[i]*ptOmegaPrimary[i]);
      sigma_primary_sub[i][h] = (meanXPrimary[i]*sumPtOmegaPrimary[i] + meanSumPtPrimary[i]*varXPrimary[i]/meanXPrimary[i] - 2*(meanMultPtPrimary[i]-meanXPrimary[i]*meanSumPtPrimary[i])) / (meanXPrimary[i]*ptOmegaPrimary[i]);

      tmp_1 = (meanSumPt3Primary[i] - 3*meanSumPtcomma2PtPrimary[i] +2*meanSumPtcomma3Primary[i])/varX3Primary[i];
      tmp_2 = meanPtPrimary[i]*(meanMultPt2Primary[i] - 2*meanSumPt2Primary[i] - meanMultPtcomma2Primary[i] + 2*meanSumPtcomma2Primary[i])/varX3Primary[i];
      tmp_3 = meanPtPrimary[i]*meanPtPrimary[i]*(meanMultMultPtPrimary[i] - 3*meanMultPtPrimary[i] +2*meanSumPtPrimary[i])/varX3Primary[i];
      gammapt_primary_sub[i][h] = tmp_1 - 3*tmp_2 + 3*tmp_3 - meanPtPrimary[i]*meanPtPrimary[i]*meanPtPrimary[i];
      gammapt_primary_sub[i][h] = meanPtPrimary[i]*gammapt_primary_sub[i][h]/(dptdpt_primary_sub[i][h]*dptdpt_primary_sub[i][h]);

      corr_primary_sub[i][h] = ((meanEBESumPtPrimary[i]/meanEBEPtPrimary[i])-meanEBEXPrimary[i])/omegaEBE;






    };


  };


  for (int i=0; i<10; i++) {
    for (Int_t h=0; h<nSubSamples; h++) {
      mult_sim_mean[i] = mult_sim_mean[i] + mult_sim_sub[i][h];
      meanpt_star_sim_mean[i] = meanpt_star_sim_mean[i] + meanpt_star_sim_sub[i][h];
      meanpt_alice_sim_mean[i] = meanpt_alice_sim_mean[i] + meanpt_alice_sim_sub[i][h];
      delta_sim_mean[i] = delta_sim_mean[i] + delta_sim_sub[i][h];
      sigma_sim_mean[i] = sigma_sim_mean[i] + sigma_sim_sub[i][h];
      dptdpt_sim_mean[i] = dptdpt_sim_mean[i] + dptdpt_sim_sub[i][h];
      dptdpt_star_sim_mean[i] = dptdpt_star_sim_mean[i] + dptdpt_star_sim_sub[i][h];
      gammapt_star_sim_mean[i] = gammapt_star_sim_mean[i] + gammapt_star_sim_sub[i][h];
      gammapt_sim_mean[i] = gammapt_sim_mean[i] + gammapt_sim_sub[i][h];
      corr_sim_mean[i] = corr_sim_mean[i] + corr_sim_sub[i][h];

      mult_global_mean[i] = mult_global_mean[i] + mult_global_sub[i][h];
      meanpt_star_global_mean[i] = meanpt_star_global_mean[i] + meanpt_star_global_sub[i][h];
      meanpt_alice_global_mean[i] = meanpt_alice_global_mean[i] + meanpt_alice_global_sub[i][h];
      delta_global_mean[i] = delta_global_mean[i] + delta_global_sub[i][h];
      sigma_global_mean[i] = sigma_global_mean[i] + sigma_global_sub[i][h];
      dptdpt_global_mean[i] = dptdpt_global_mean[i] + dptdpt_global_sub[i][h];
      dptdpt_star_global_mean[i] = dptdpt_star_global_mean[i] + dptdpt_star_global_sub[i][h];
      gammapt_star_global_mean[i] = gammapt_star_global_mean[i] + gammapt_star_global_sub[i][h];
      gammapt_global_mean[i] = gammapt_global_mean[i] + gammapt_global_sub[i][h];
      corr_global_mean[i] = corr_global_mean[i] + corr_global_sub[i][h];
      //cout << i << "; " << h << " = " <<  dptdpt_global_sub[i][h] << endl;

      mult_primary_mean[i] = mult_primary_mean[i] + mult_primary_sub[i][h];
      meanpt_star_primary_mean[i] = meanpt_star_primary_mean[i] + meanpt_star_primary_sub[i][h];
      meanpt_alice_primary_mean[i] = meanpt_alice_primary_mean[i] + meanpt_alice_primary_sub[i][h];
      delta_primary_mean[i] = delta_primary_mean[i] + delta_primary_sub[i][h];
      sigma_primary_mean[i] = sigma_primary_mean[i] + sigma_primary_sub[i][h];
      dptdpt_primary_mean[i] = dptdpt_primary_mean[i] + dptdpt_primary_sub[i][h];
      dptdpt_star_primary_mean[i] = dptdpt_star_primary_mean[i] + dptdpt_star_primary_sub[i][h];
      gammapt_star_primary_mean[i] = gammapt_star_primary_mean[i] + gammapt_star_primary_sub[i][h];
      gammapt_primary_mean[i] = gammapt_primary_mean[i] + gammapt_primary_sub[i][h];
      corr_primary_mean[i] = corr_primary_mean[i] + corr_primary_sub[i][h];

      //totalMultMean[i] = totalMultMean[i] + totalMultSub[i][h];
      //totalOmegaMean[i] = totalOmegaMean[i] + totalOmegaSub[i][h];
      //totalDeltaMean[i] = totalDeltaMean[i] + totalDeltaSub[i][h];
      //totalSigmaMean[i] = totalSigmaMean[i] + totalSigmaSub[i][h];
    }

    mult_sim_mean[i] = mult_sim_mean[i]/nSubSamplesDbl;
    meanpt_star_sim_mean[i] = meanpt_star_sim_mean[i]/nSubSamplesDbl;
    meanpt_alice_sim_mean[i] = meanpt_alice_sim_mean[i]/nSubSamplesDbl;
    delta_sim_mean[i] = delta_sim_mean[i]/nSubSamplesDbl;
    sigma_sim_mean[i] = sigma_sim_mean[i]/nSubSamplesDbl;
    dptdpt_sim_mean[i] = dptdpt_sim_mean[i]/nSubSamplesDbl;
    dptdpt_star_sim_mean[i] = dptdpt_star_sim_mean[i]/nSubSamplesDbl;
    gammapt_star_sim_mean[i] = gammapt_star_sim_mean[i]/nSubSamplesDbl;
    gammapt_sim_mean[i] = gammapt_sim_mean[i]/nSubSamplesDbl;
    corr_sim_mean[i] = corr_sim_mean[i]/nSubSamplesDbl;

    mult_global_mean[i] = mult_global_mean[i]/nSubSamplesDbl;
    meanpt_star_global_mean[i] = meanpt_star_global_mean[i]/nSubSamplesDbl;
    meanpt_alice_global_mean[i] = meanpt_alice_global_mean[i]/nSubSamplesDbl;
    delta_global_mean[i] = delta_global_mean[i]/nSubSamplesDbl;
    sigma_global_mean[i] = sigma_global_mean[i]/nSubSamplesDbl;
    dptdpt_global_mean[i] = dptdpt_global_mean[i]/nSubSamplesDbl;
    dptdpt_star_global_mean[i] = dptdpt_star_global_mean[i]/nSubSamplesDbl;
    gammapt_star_global_mean[i] = gammapt_star_global_mean[i]/nSubSamplesDbl;
    gammapt_global_mean[i] = gammapt_global_mean[i]/nSubSamplesDbl;
    corr_global_mean[i] = corr_global_mean[i]/nSubSamplesDbl;

    mult_primary_mean[i] = mult_primary_mean[i]/nSubSamplesDbl;
    meanpt_star_primary_mean[i] = meanpt_star_primary_mean[i]/nSubSamplesDbl;
    meanpt_alice_primary_mean[i] = meanpt_alice_primary_mean[i]/nSubSamplesDbl;
    delta_primary_mean[i] = delta_primary_mean[i]/nSubSamplesDbl;
    sigma_primary_mean[i] = sigma_primary_mean[i]/nSubSamplesDbl;
    dptdpt_primary_mean[i] = dptdpt_primary_mean[i]/nSubSamplesDbl;
    dptdpt_star_primary_mean[i] = dptdpt_star_primary_mean[i]/nSubSamplesDbl;
    gammapt_star_primary_mean[i] = gammapt_star_primary_mean[i]/nSubSamplesDbl;
    gammapt_primary_mean[i] = gammapt_primary_mean[i]/nSubSamplesDbl;
    corr_primary_mean[i] = corr_primary_mean[i]/nSubSamplesDbl;


    for (Int_t h=0; h<nSubSamples; h++) {
      mult_sim_error[i] += (mult_sim_mean[i]-mult_sim_sub[i][h])*(mult_sim_mean[i]-mult_sim_sub[i][h]);
      meanpt_star_sim_error[i] += (meanpt_star_sim_mean[i]-meanpt_star_sim_sub[i][h])*(meanpt_star_sim_mean[i]-meanpt_star_sim_sub[i][h]);
      meanpt_alice_sim_error[i] += (meanpt_alice_sim_mean[i]-meanpt_alice_sim_sub[i][h])*(meanpt_alice_sim_mean[i]-meanpt_alice_sim_sub[i][h]);
      delta_sim_error[i] += (delta_sim_mean[i]-delta_sim_sub[i][h])*(delta_sim_mean[i]-delta_sim_sub[i][h]);
      sigma_sim_error[i] += (sigma_sim_mean[i]-sigma_sim_sub[i][h])*(sigma_sim_mean[i]-sigma_sim_sub[i][h]);
      dptdpt_sim_error[i] += (dptdpt_sim_mean[i]-dptdpt_sim_sub[i][h])*(dptdpt_sim_mean[i]-dptdpt_sim_sub[i][h]);
      dptdpt_star_sim_error[i] += (dptdpt_star_sim_mean[i]-dptdpt_star_sim_sub[i][h])*(dptdpt_star_sim_mean[i]-dptdpt_star_sim_sub[i][h]);
      gammapt_star_sim_error[i] += (gammapt_star_sim_mean[i]-gammapt_star_sim_sub[i][h])*(gammapt_star_sim_mean[i]-gammapt_star_sim_sub[i][h]);
      gammapt_sim_error[i] += (gammapt_sim_mean[i]-gammapt_sim_sub[i][h])*(gammapt_sim_mean[i]-gammapt_sim_sub[i][h]);
      corr_sim_error[i] += (corr_sim_mean[i]-corr_sim_sub[i][h])*(corr_sim_mean[i]-corr_sim_sub[i][h]);



      mult_global_error[i] += (mult_global_mean[i]-mult_global_sub[i][h])*(mult_global_mean[i]-mult_global_sub[i][h]);
      meanpt_star_global_error[i] += (meanpt_star_global_mean[i]-meanpt_star_global_sub[i][h])*(meanpt_star_global_mean[i]-meanpt_star_global_sub[i][h]);
      meanpt_alice_global_error[i] += (meanpt_alice_global_mean[i]-meanpt_alice_global_sub[i][h])*(meanpt_alice_global_mean[i]-meanpt_alice_global_sub[i][h]);
      delta_global_error[i] += (delta_global_mean[i]-delta_global_sub[i][h])*(delta_global_mean[i]-delta_global_sub[i][h]);
      sigma_global_error[i] += (sigma_global_mean[i]-sigma_global_sub[i][h])*(sigma_global_mean[i]-sigma_global_sub[i][h]);
      dptdpt_global_error[i] += (dptdpt_global_mean[i]-dptdpt_global_sub[i][h])*(dptdpt_global_mean[i]-dptdpt_global_sub[i][h]);
      dptdpt_star_global_error[i] += (dptdpt_star_global_mean[i]-dptdpt_star_global_sub[i][h])*(dptdpt_star_global_mean[i]-dptdpt_star_global_sub[i][h]);
      gammapt_star_global_error[i] += (gammapt_star_global_mean[i]-gammapt_star_global_sub[i][h])*(gammapt_star_global_mean[i]-gammapt_star_global_sub[i][h]);
      gammapt_global_error[i] += (gammapt_global_mean[i]-gammapt_global_sub[i][h])*(gammapt_global_mean[i]-gammapt_global_sub[i][h]);
      corr_global_error[i] += (corr_global_mean[i]-corr_global_sub[i][h])*(corr_global_mean[i]-corr_global_sub[i][h]);

      mult_primary_error[i] += (mult_primary_mean[i]-mult_primary_sub[i][h])*(mult_primary_mean[i]-mult_primary_sub[i][h]);
      meanpt_star_primary_error[i] += (meanpt_star_primary_mean[i]-meanpt_star_primary_sub[i][h])*(meanpt_star_primary_mean[i]-meanpt_star_primary_sub[i][h]);
      meanpt_alice_primary_error[i] += (meanpt_alice_primary_mean[i]-meanpt_alice_primary_sub[i][h])*(meanpt_alice_primary_mean[i]-meanpt_alice_primary_sub[i][h]);
      delta_primary_error[i] += (delta_primary_mean[i]-delta_primary_sub[i][h])*(delta_primary_mean[i]-delta_primary_sub[i][h]);
      sigma_primary_error[i] += (sigma_primary_mean[i]-sigma_primary_sub[i][h])*(sigma_primary_mean[i]-sigma_primary_sub[i][h]);
      dptdpt_primary_error[i] += (dptdpt_primary_mean[i]-dptdpt_primary_sub[i][h])*(dptdpt_primary_mean[i]-dptdpt_primary_sub[i][h]);
      dptdpt_star_primary_error[i] += (dptdpt_star_primary_mean[i]-dptdpt_star_primary_sub[i][h])*(dptdpt_star_primary_mean[i]-dptdpt_star_primary_sub[i][h]);
      gammapt_star_primary_error[i] += (gammapt_star_primary_mean[i]-gammapt_star_primary_sub[i][h])*(gammapt_star_primary_mean[i]-gammapt_star_primary_sub[i][h]);
      gammapt_primary_error[i] += (gammapt_primary_mean[i]-gammapt_primary_sub[i][h])*(gammapt_primary_mean[i]-gammapt_primary_sub[i][h]);
      corr_primary_error[i] += (corr_primary_mean[i]-corr_primary_sub[i][h])*(corr_primary_mean[i]-corr_primary_sub[i][h]);

    }

    mult_sim_error[i] = sqrt(mult_sim_error[i])/dtmp2;
    meanpt_star_sim_error[i] = sqrt(meanpt_star_sim_error[i])/dtmp2;
    meanpt_alice_sim_error[i] = sqrt(meanpt_alice_sim_error[i])/dtmp2;
    delta_sim_error[i] = sqrt(delta_sim_error[i])/dtmp2;
    sigma_sim_error[i] = sqrt(sigma_sim_error[i])/dtmp2;
    dptdpt_sim_error[i] = sqrt(dptdpt_sim_error[i])/dtmp2;
    dptdpt_star_sim_error[i] = sqrt(dptdpt_star_sim_error[i])/dtmp2;
    gammapt_star_sim_error[i] = sqrt(gammapt_star_sim_error[i])/dtmp2;
    gammapt_sim_error[i] = sqrt(gammapt_sim_error[i])/dtmp2;
    corr_sim_error[i] = sqrt(corr_sim_error[i])/dtmp2;




    mult_global_error[i] = sqrt(mult_global_error[i])/dtmp2;
    meanpt_star_global_error[i] = sqrt(meanpt_star_global_error[i])/dtmp2;
    meanpt_alice_global_error[i] = sqrt(meanpt_alice_global_error[i])/dtmp2;
    delta_global_error[i] = sqrt(delta_global_error[i])/dtmp2;
    sigma_global_error[i] = sqrt(sigma_global_error[i])/dtmp2;
    dptdpt_global_error[i] = sqrt(dptdpt_global_error[i])/dtmp2;
    dptdpt_star_global_error[i] = sqrt(dptdpt_star_global_error[i])/dtmp2;
    gammapt_star_global_error[i] = sqrt(gammapt_star_global_error[i])/dtmp2;
    gammapt_global_error[i] = sqrt(gammapt_global_error[i])/dtmp2;
    corr_global_error[i] = sqrt(corr_global_error[i])/dtmp2;


    mult_primary_error[i] = sqrt(mult_primary_error[i])/dtmp2;
    meanpt_star_primary_error[i] = sqrt(meanpt_star_primary_error[i])/dtmp2;
    meanpt_alice_primary_error[i] = sqrt(meanpt_alice_primary_error[i])/dtmp2;
    delta_primary_error[i] = sqrt(delta_primary_error[i])/dtmp2;
    sigma_primary_error[i] = sqrt(sigma_primary_error[i])/dtmp2;
    dptdpt_primary_error[i] = sqrt(dptdpt_primary_error[i])/dtmp2;
    dptdpt_star_primary_error[i] = sqrt(dptdpt_star_primary_error[i])/dtmp2;
    gammapt_star_primary_error[i] = sqrt(gammapt_star_primary_error[i])/dtmp2;
    gammapt_primary_error[i] = sqrt(gammapt_primary_error[i])/dtmp2;
    corr_primary_error[i] = sqrt(corr_primary_error[i])/dtmp2;





    cout << "Delta eta = " << 0.2+i*0.2 << endl;
    cout << "mult_sim_mean[" << i << "] = " << mult_sim_mean[i] << " ; +/- " << mult_sim_error[i] << endl;
    cout << "meanpt_star_sim_mean[" << i << "] = " << meanpt_star_sim_mean[i] << " ; +/- " << meanpt_star_sim_error[i] << endl;
    cout << "meanpt_alice_sim_mean[" << i << "] = " << meanpt_alice_sim_mean[i] << " ; +/- " << meanpt_alice_sim_error[i] << endl;
    cout << "delta_sim_mean[" << i << "] = " << delta_sim_mean[i] << " ; +/- " << delta_sim_error[i] << endl;
    cout << "sigma_sim_mean[" << i << "] = " << sigma_sim_mean[i] << " ; +/- " << sigma_sim_error[i] << endl;
    cout << "dptdpt_star_sim_mean[" << i << "] = " << dptdpt_star_sim_mean[i] << " ; +/- " << dptdpt_star_sim_error[i] << endl;
    cout << "dptdpt_sim_mean[" << i << "] = " << dptdpt_sim_mean[i] << " ; +/- " << dptdpt_sim_error[i] << endl;
    cout << "gammapt_star_sim_mean[" << i << "] = " << gammapt_star_sim_mean[i] << " ; +/- " << gammapt_star_sim_error[i] << endl;
    cout << "gammapt_sim_mean[" << i << "] = " << gammapt_sim_mean[i] << " ; +/- " << gammapt_sim_error[i] << endl;
    cout << "corr_sim_mean[" << i << "] = " << corr_sim_mean[i] << " ; +/- " << corr_sim_error[i] << endl;

    cout << "mult_global_mean[" << i << "] = " << mult_global_mean[i] << " ; +/- " << mult_global_error[i] << endl;
    cout << "meanpt_star_global_mean[" << i << "] = " << meanpt_star_global_mean[i] << " ; +/- " << meanpt_star_global_error[i] << endl;
    cout << "meanpt_alice_global_mean[" << i << "] = " << meanpt_alice_global_mean[i] << " ; +/- " << meanpt_alice_global_error[i] << endl;
    cout << "delta_global_mean[" << i << "] = " << delta_global_mean[i] << " ; +/- " << delta_global_error[i] << endl;
    cout << "sigma_global_mean[" << i << "] = " << sigma_global_mean[i] << " ; +/- " << sigma_global_error[i] << endl;
    cout << "dptdpt_star_global_mean[" << i << "] = " << dptdpt_star_global_mean[i] << " ; +/- " << dptdpt_star_global_error[i] << endl;
    cout << "dptdpt_global_mean[" << i << "] = " << dptdpt_global_mean[i] << " ; +/- " << dptdpt_global_error[i] << endl;
    cout << "gammapt_star_global_mean[" << i << "] = " << gammapt_star_global_mean[i] << " ; +/- " << gammapt_star_global_error[i] << endl;
    cout << "gammapt_global_mean[" << i << "] = " << gammapt_global_mean[i] << " ; +/- " << gammapt_global_error[i] << endl;
    cout << "corr_global_mean[" << i << "] = " << corr_global_mean[i] << " ; +/- " << corr_global_error[i] << endl;

    cout << "mult_primary_mean[" << i << "] = " << mult_primary_mean[i] << " ; +/- " << mult_primary_error[i] << endl;
    cout << "meanpt_star_primary_mean[" << i << "] = " << meanpt_star_primary_mean[i] << " ; +/- " << meanpt_star_primary_error[i] << endl;
    cout << "meanpt_alice_primary_mean[" << i << "] = " << meanpt_alice_primary_mean[i] << " ; +/- " << meanpt_alice_primary_error[i] << endl;
    cout << "delta_primary_mean[" << i << "] = " << delta_primary_mean[i] << " ; +/- " << delta_primary_error[i] << endl;
    cout << "sigma_primary_mean[" << i << "] = " << sigma_primary_mean[i] << " ; +/- " << sigma_primary_error[i] << endl;
    cout << "dptdpt_star_primary_mean[" << i << "] = " << dptdpt_star_primary_mean[i] << " ; +/- " << dptdpt_star_primary_error[i] << endl;
    cout << "dptdpt_primary_mean[" << i << "] = " << dptdpt_primary_mean[i] << " ; +/- " << dptdpt_primary_error[i] << endl;
    cout << "gammapt_star_primary_mean[" << i << "] = " << gammapt_star_primary_mean[i] << " ; +/- " << gammapt_star_primary_error[i] << endl;
    cout << "gammapt_primary_mean[" << i << "] = " << gammapt_primary_mean[i] << " ; +/- " << gammapt_primary_error[i] << endl;
    cout << "corr_primary_mean[" << i << "] = " << corr_primary_mean[i] << " ; +/- " << corr_primary_error[i] << endl;



  }

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();
  pTHist[0][0]->Draw("Hist");
  c1->SaveAs("pt.pdf");

  TCanvas* c2 = new TCanvas("c2","c2",1200,800);
  c2->cd();
  pTMultHist[0][0]->Draw("colz");
  c2->SaveAs("Npt.pdf");

  TCanvas* c3 = new TCanvas("c3","c3",1200,800);
  c3->cd();
  meanpTHist[0][0]->Draw("Hist");
  c3->SaveAs("meanpt.pdf");



}


void analyzeTreeWithSubsamples(TTree* extractedTree, Int_t nSubSamples, Int_t nEventsInSubSample);



void analyzeTreeWithSubsamples(TTree* extractedTree, Int_t nSubSamples, Int_t nEventsInSubSample)
{

  gStyle->SetTitleSize(0.06, "xz");
  gStyle->SetTitleSize(0.07, "y");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetLabelSize(0.06, "xyz");

  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.18);
  gStyle->SetPadTopMargin(0.18);
  gStyle->SetPadBottomMargin(0.18);


  TLegend *leg = new TLegend(0.68,0.92,0.88,0.985);
  TString legEntrySim = "Pure SMASH";
  TString legEntryRec = "Reconstructed SMASH";
  TString legEntryRecPrimary = "Reconstructed Primary SMASH";
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.04);

  TLatex label;
  //label.SetTextFont(132);
  label.SetTextSize(0.045);
  label.SetTextAlign(22);
  label.SetNDC(kTRUE);

  TString strEnergy("_9point46gev");

  char strLowCentrality[35];
  int lowCentrality = 0;
  sprintf(strLowCentrality, "%d", lowCentrality);
  TString strrLowCentrality(strLowCentrality);

  char strHighCentrality[35];
  int highCentrality = 10;
  sprintf(strHighCentrality, "%d", highCentrality);
  TString strrHighCentrality(strHighCentrality);

  char strNSub[35];
  sprintf(strNSub, "%d", nSubSamples);
  TString strrNSub(strNSub);


  std::string beginnning;
  beginnning = std::string("double ")+std::string("b_")+std::string(strrLowCentrality)+std::string("_")+std::string(strrHighCentrality)+std::string(strEnergy)+std::string("_");

  std::string shortbeginnning;
  shortbeginnning = std::string("b_")+std::string(strrLowCentrality)+std::string("_")+std::string(strrHighCentrality)+std::string(strEnergy)+std::string("_");
  TString rootshortbeginning;
  rootshortbeginning = (TString)shortbeginnning;

  std::string titleOutput;
  titleOutput = std::string("output_")+std::string("b_")+std::string(strrLowCentrality)+std::string("_")+std::string(strrHighCentrality)+std::string(strEnergy)+std::string(".txt");



  Double_t nSubSamplesDbl;
  nSubSamplesDbl = (Double_t)nSubSamples;

  Double_t dtmp2 = sqrt(nSubSamplesDbl*(nSubSamplesDbl-1)); //for subsample


  Int_t nWindows;
  nWindows = 16;

  Int_t nOneWindows;
  nOneWindows = 32;

  Int_t forwMult[nWindows];
  Int_t backMult[nWindows];
  Int_t totalMult[nOneWindows];
  Double_t totalPt[nOneWindows];
  Double_t totalPt2[nOneWindows];

  Int_t forwMultPlus[nWindows];
  Int_t backMultPlus[nWindows];
  Int_t plusMult[nOneWindows];
  Double_t plusPt[nOneWindows];
  Double_t plusPt2[nOneWindows];

  Int_t forwMultMinus[nWindows];
  Int_t backMultMinus[nWindows];
  Int_t minusMult[nOneWindows];
  Double_t minusPt[nOneWindows];
  Double_t minusPt2[nOneWindows];


  Int_t forwMultPrimary[nWindows];
  Int_t backMultPrimary[nWindows];
  Int_t totalMultPrimary[nOneWindows];
  Double_t totalPtPrimary[nOneWindows];
  Double_t totalPt2Primary[nOneWindows];

  Int_t forwMultPlusPrimary[nWindows];
  Int_t backMultPlusPrimary[nWindows];
  Int_t plusMultPrimary[nOneWindows];
  Double_t plusPtPrimary[nOneWindows];
  Double_t plusPt2Primary[nOneWindows];

  Int_t forwMultMinusPrimary[nWindows];
  Int_t backMultMinusPrimary[nWindows];
  Int_t minusMultPrimary[nOneWindows];
  Double_t minusPtPrimary[nOneWindows];
  Double_t minusPt2Primary[nOneWindows];

  Int_t forwMultSim[nWindows];
  Int_t backMultSim[nWindows];
  Int_t totalMultSim[nOneWindows];
  Double_t totalPtSim[nOneWindows];
  Double_t totalPt2Sim[nOneWindows];

  Int_t forwMultPlusSim[nWindows];
  Int_t backMultPlusSim[nWindows];
  Int_t plusMultSim[nOneWindows];
  Double_t plusPtSim[nOneWindows];
  Double_t plusPt2Sim[nOneWindows];

  Int_t forwMultMinusSim[nWindows];
  Int_t backMultMinusSim[nWindows];
  Int_t minusMultSim[nOneWindows];
  Double_t minusPtSim[nOneWindows];
  Double_t minusPt2Sim[nOneWindows];


  Int_t eventId;
  Float_t impactPar;
  Float_t vertexX, vertexY, vertexZ;
  Int_t f_SIM_NumberOfTracks; //!
  Float_t f_SIM_TrackPt[NMaxTrack],f_SIM_TrackEta[NMaxTrack];
  Int_t f_SIM_TrackCharge[NMaxTrack];
  Int_t fNumberOfTracks; //!
  Float_t fTrackPt[NMaxTrack],fTrackEta[NMaxTrack];
  Float_t fTrackPtPrimary[NMaxTrack],fTrackEtaPrimary[NMaxTrack];
  Int_t fTrackCharge[NMaxTrack];
  Int_t fnTPChits[NMaxTrack];
  Bool_t fHasPrimaryTrack[NMaxTrack];
  extractedTree->SetBranchAddress("fEventId",&eventId);
  extractedTree->SetBranchAddress("fImpactPar",&impactPar);
  //extractedTree->SetBranchAddress("fVertexX",&vertexX);
  //extractedTree->SetBranchAddress("fVertexY",&vertexY);
  //extractedTree->SetBranchAddress("fVertexZ",&vertexZ);
  extractedTree->SetBranchAddress("SIM_NumberOfTracks",&f_SIM_NumberOfTracks);
  extractedTree->SetBranchAddress("SIM_TrackPt",f_SIM_TrackPt);
  extractedTree->SetBranchAddress("SIM_TrackEta",f_SIM_TrackEta);
  extractedTree->SetBranchAddress("SIM_TrackCharge",f_SIM_TrackCharge);
  extractedTree->SetBranchAddress("nTracks",&fNumberOfTracks);
  extractedTree->SetBranchAddress("trackPt",fTrackPt);
  extractedTree->SetBranchAddress("trackEta",fTrackEta);
  extractedTree->SetBranchAddress("trackPtPrimary",fTrackPtPrimary);
  extractedTree->SetBranchAddress("trackEtaPrimary",fTrackEtaPrimary);
  extractedTree->SetBranchAddress("trackCharge",fTrackCharge);
  //extractedTree->SetBranchAddress("nTPChits",fnTPChits);
  //extractedTree->SetBranchAddress("hasPrimaryTrack",fHasPrimaryTrack);

  TH2D* corrTotalHist[16][nSubSamples];

  TH2D* corrPlusPlusHist[16][nSubSamples];
  TH2D* corrMinusMinusHist[16][nSubSamples];
  TH2D* corrPlusMinusHist[16][nSubSamples];
  TH2D* corrMinusPlusHist[16][nSubSamples];


  TH2D* corrTotalPrimaryHist[16][nSubSamples];

  TH2D* corrPlusPlusPrimaryHist[16][nSubSamples];
  TH2D* corrMinusMinusPrimaryHist[16][nSubSamples];
  TH2D* corrPlusMinusPrimaryHist[16][nSubSamples];
  TH2D* corrMinusPlusPrimaryHist[16][nSubSamples];



  TH2D* corrTotalSimHist[16][nSubSamples];

  TH2D* corrPlusPlusSimHist[16][nSubSamples];
  TH2D* corrMinusMinusSimHist[16][nSubSamples];
  TH2D* corrPlusMinusSimHist[16][nSubSamples];
  TH2D* corrMinusPlusSimHist[16][nSubSamples];

  for (int i=0; i<16; i++) {
    char str[15];
    sprintf(str, "%d", i);
    TString strr(str);

    for (Int_t h=0; h<nSubSamples; h++) {
      char strSub[15];
      sprintf(strSub, "%d", h);
      TString strrSub(strSub);

      corrTotalHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrTotalHist",strrNSub+"x"+strr+"x"+strrSub+"corrTotalHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);

      corrPlusPlusHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrPlusPlusHist",strrNSub+"x"+strr+"x"+strrSub+"corrPlusPlusHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrMinusMinusHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrMinusMinusHist",strrNSub+"x"+strr+"x"+strrSub+"corrMinusMinusHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrPlusMinusHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrPlusMinusHist",strrNSub+"x"+strr+"x"+strrSub+"corrPlusMinusHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrMinusPlusHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrMinusPlusHist",strrNSub+"x"+strr+"x"+strrSub+"corrMinusPlusHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);




      corrTotalPrimaryHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrTotalPrimaryHist",strrNSub+"x"+strr+"x"+strrSub+"corrTotalPrimaryHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);

      corrPlusPlusPrimaryHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrPlusPlusPrimaryHist",strrNSub+"x"+strr+"x"+strrSub+"corrPlusPlusPrimaryHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrMinusMinusPrimaryHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrMinusMinusPrimaryHist",strrNSub+"x"+strr+"x"+strrSub+"corrMinusMinusPrimaryHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrPlusMinusPrimaryHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrPlusMinusPrimaryHist",strrNSub+"x"+strr+"x"+strrSub+"corrPlusMinusPrimaryHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrMinusPlusPrimaryHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrMinusPlusPrimaryHist",strrNSub+"x"+strr+"x"+strrSub+"corrMinusPlusPrimaryHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);




      corrTotalSimHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrTotalSimHist",strrNSub+"x"+strr+"x"+strrSub+"corrTotalSimHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);

      corrPlusPlusSimHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrPlusPlusSimHist",strrNSub+"x"+strr+"x"+strrSub+"corrPlusPlusSimHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrMinusMinusSimHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrMinusMinusSimHist",strrNSub+"x"+strr+"x"+strrSub+"corrMinusMinusSimHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrPlusMinusSimHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrPlusMinusSimHist",strrNSub+"x"+strr+"x"+strrSub+"corrPlusMinusSimHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);
      corrMinusPlusSimHist[i][h] = new TH2D(strrNSub+"x"+strr+"x"+strrSub+"corrMinusPlusSimHist",strrNSub+"x"+strr+"x"+strrSub+"corrMinusPlusSimHist;nF;nB;entries",50,-0.5,49.5,50,-0.5,49.5);


    }

  }




  Int_t nStartingEvent, nEndingEvent;

  nStartingEvent = 0;
  nEndingEvent = nStartingEvent+nEventsInSubSample;


  long double totalMultSub[32][nSubSamples];
  long double totalMultPrimarySub[32][nSubSamples];
  long double totalMultSimSub[32][nSubSamples];


  long double totalMultMean[32]={0};
  long double totalMultPrimaryMean[32]={0};
  long double totalMultSimMean[32]={0};


  long double totalMultError[32]={0};
  long double totalMultPrimaryError[32]={0};
  long double totalMultSimError[32]={0};




  long double totalOmegaSub[32][nSubSamples];
  long double totalOmegaPrimarySub[32][nSubSamples];
  long double totalOmegaSimSub[32][nSubSamples];


  long double totalOmegaMean[32]={0};
  long double totalOmegaPrimaryMean[32]={0};
  long double totalOmegaSimMean[32]={0};


  long double totalOmegaError[32]={0};
  long double totalOmegaPrimaryError[32]={0};
  long double totalOmegaSimError[32]={0};



  long double totalDeltaSub[32][nSubSamples];
  long double totalDeltaPrimarySub[32][nSubSamples];
  long double totalDeltaSimSub[32][nSubSamples];


  long double totalDeltaMean[32]={0};
  long double totalDeltaPrimaryMean[32]={0};
  long double totalDeltaSimMean[32]={0};


  long double totalDeltaError[32]={0};
  long double totalDeltaPrimaryError[32]={0};
  long double totalDeltaSimError[32]={0};



  long double totalSigmaSub[32][nSubSamples];
  long double totalSigmaPrimarySub[32][nSubSamples];
  long double totalSigmaSimSub[32][nSubSamples];


  long double totalSigmaMean[32]={0};
  long double totalSigmaPrimaryMean[32]={0};
  long double totalSigmaSimMean[32]={0};


  long double totalSigmaError[32]={0};
  long double totalSigmaPrimaryError[32]={0};
  long double totalSigmaSimError[32]={0};











  long double plusMultSub[32][nSubSamples];
  long double plusMultPrimarySub[32][nSubSamples];
  long double plusMultSimSub[32][nSubSamples];


  long double plusMultMean[32]={0};
  long double plusMultPrimaryMean[32]={0};
  long double plusMultSimMean[32]={0};


  long double plusMultError[32]={0};
  long double plusMultPrimaryError[32]={0};
  long double plusMultSimError[32]={0};




  long double plusOmegaSub[32][nSubSamples];
  long double plusOmegaPrimarySub[32][nSubSamples];
  long double plusOmegaSimSub[32][nSubSamples];


  long double plusOmegaMean[32]={0};
  long double plusOmegaPrimaryMean[32]={0};
  long double plusOmegaSimMean[32]={0};


  long double plusOmegaError[32]={0};
  long double plusOmegaPrimaryError[32]={0};
  long double plusOmegaSimError[32]={0};



  long double plusDeltaSub[32][nSubSamples];
  long double plusDeltaPrimarySub[32][nSubSamples];
  long double plusDeltaSimSub[32][nSubSamples];


  long double plusDeltaMean[32]={0};
  long double plusDeltaPrimaryMean[32]={0};
  long double plusDeltaSimMean[32]={0};


  long double plusDeltaError[32]={0};
  long double plusDeltaPrimaryError[32]={0};
  long double plusDeltaSimError[32]={0};



  long double plusSigmaSub[32][nSubSamples];
  long double plusSigmaPrimarySub[32][nSubSamples];
  long double plusSigmaSimSub[32][nSubSamples];


  long double plusSigmaMean[32]={0};
  long double plusSigmaPrimaryMean[32]={0};
  long double plusSigmaSimMean[32]={0};


  long double plusSigmaError[32]={0};
  long double plusSigmaPrimaryError[32]={0};
  long double plusSigmaSimError[32]={0};













  long double minusMultSub[32][nSubSamples];
  long double minusMultPrimarySub[32][nSubSamples];
  long double minusMultSimSub[32][nSubSamples];


  long double minusMultMean[32]={0};
  long double minusMultPrimaryMean[32]={0};
  long double minusMultSimMean[32]={0};


  long double minusMultError[32]={0};
  long double minusMultPrimaryError[32]={0};
  long double minusMultSimError[32]={0};




  long double minusOmegaSub[32][nSubSamples];
  long double minusOmegaPrimarySub[32][nSubSamples];
  long double minusOmegaSimSub[32][nSubSamples];


  long double minusOmegaMean[32]={0};
  long double minusOmegaPrimaryMean[32]={0};
  long double minusOmegaSimMean[32]={0};


  long double minusOmegaError[32]={0};
  long double minusOmegaPrimaryError[32]={0};
  long double minusOmegaSimError[32]={0};



  long double minusDeltaSub[32][nSubSamples];
  long double minusDeltaPrimarySub[32][nSubSamples];
  long double minusDeltaSimSub[32][nSubSamples];


  long double minusDeltaMean[32]={0};
  long double minusDeltaPrimaryMean[32]={0};
  long double minusDeltaSimMean[32]={0};


  long double minusDeltaError[32]={0};
  long double minusDeltaPrimaryError[32]={0};
  long double minusDeltaSimError[32]={0};



  long double minusSigmaSub[32][nSubSamples];
  long double minusSigmaPrimarySub[32][nSubSamples];
  long double minusSigmaSimSub[32][nSubSamples];


  long double minusSigmaMean[32]={0};
  long double minusSigmaPrimaryMean[32]={0};
  long double minusSigmaSimMean[32]={0};


  long double minusSigmaError[32]={0};
  long double minusSigmaPrimaryError[32]={0};
  long double minusSigmaSimError[32]={0};

















  for (int h=0; h<nSubSamples; h++) {
    cout << "subsample #" << h << endl;
    long double nEventsInSubSample_D;
    nEventsInSubSample_D=0;

    long double fullSumPt[32] = {0};
    long double fullSumPt2[32] = {0};
    long double fullSumPtpow2[32] = {0};
    long double fullSumMult[32] = {0};
    long double fullSumMultMult[32] = {0};
    long double fullSumMultPt[32] = {0};

    long double fullSumPtPrimary[32] = {0};
    long double fullSumPt2Primary[32] = {0};
    long double fullSumPtpow2Primary[32] = {0};
    long double fullSumMultPrimary[32] = {0};
    long double fullSumMultMultPrimary[32] = {0};
    long double fullSumMultPtPrimary[32] = {0};

    long double fullSumPtSim[32] = {0};
    long double fullSumPt2Sim[32] = {0};
    long double fullSumPtpow2Sim[32] = {0};
    long double fullSumMultSim[32] = {0};
    long double fullSumMultMultSim[32] = {0};
    long double fullSumMultPtSim[32] = {0};



    long double fullPlusSumPt[32] = {0};
    long double fullPlusSumPt2[32] = {0};
    long double fullPlusSumPtpow2[32] = {0};
    long double fullPlusSumMult[32] = {0};
    long double fullPlusSumMultMult[32] = {0};
    long double fullPlusSumMultPt[32] = {0};

    long double fullPlusSumPtPrimary[32] = {0};
    long double fullPlusSumPt2Primary[32] = {0};
    long double fullPlusSumPtpow2Primary[32] = {0};
    long double fullPlusSumMultPrimary[32] = {0};
    long double fullPlusSumMultMultPrimary[32] = {0};
    long double fullPlusSumMultPtPrimary[32] = {0};

    long double fullPlusSumPtSim[32] = {0};
    long double fullPlusSumPt2Sim[32] = {0};
    long double fullPlusSumPtpow2Sim[32] = {0};
    long double fullPlusSumMultSim[32] = {0};
    long double fullPlusSumMultMultSim[32] = {0};
    long double fullPlusSumMultPtSim[32] = {0};



    long double fullMinusSumPt[32] = {0};
    long double fullMinusSumPt2[32] = {0};
    long double fullMinusSumPtpow2[32] = {0};
    long double fullMinusSumMult[32] = {0};
    long double fullMinusSumMultMult[32] = {0};
    long double fullMinusSumMultPt[32] = {0};

    long double fullMinusSumPtPrimary[32] = {0};
    long double fullMinusSumPt2Primary[32] = {0};
    long double fullMinusSumPtpow2Primary[32] = {0};
    long double fullMinusSumMultPrimary[32] = {0};
    long double fullMinusSumMultMultPrimary[32] = {0};
    long double fullMinusSumMultPtPrimary[32] = {0};

    long double fullMinusSumPtSim[32] = {0};
    long double fullMinusSumPt2Sim[32] = {0};
    long double fullMinusSumPtpow2Sim[32] = {0};
    long double fullMinusSumMultSim[32] = {0};
    long double fullMinusSumMultMultSim[32] = {0};
    long double fullMinusSumMultPtSim[32] = {0};



    long double meanSumPt[32]={0};
    long double meanSumPt2[32]={0};
    long double meanPt[32]={0};
    long double meanPt2[32]={0};
    long double meanMultPt[32]={0};
    long double meanX2[32]={0};
    long double meanX[32]={0};
    long double varX[32]={0};
    long double ptOmega[32]={0};
    long double sumPtOmega[32]={0};

    long double meanSumPtPrimary[32]={0};
    long double meanSumPt2Primary[32]={0};
    long double meanPtPrimary[32]={0};
    long double meanPt2Primary[32]={0};
    long double meanMultPtPrimary[32]={0};
    long double meanX2Primary[32]={0};
    long double meanXPrimary[32]={0};
    long double varXPrimary[32]={0};
    long double ptOmegaPrimary[32]={0};
    long double sumPtOmegaPrimary[32]={0};


    long double meanSumPtSim[32]={0};
    long double meanSumPt2Sim[32]={0};
    long double meanPtSim[32]={0};
    long double meanPt2Sim[32]={0};
    long double meanMultPtSim[32]={0};
    long double meanX2Sim[32]={0};
    long double meanXSim[32]={0};
    long double varXSim[32]={0};
    long double ptOmegaSim[32]={0};
    long double sumPtOmegaSim[32]={0};




    long double plusmeanSumPt[32]={0};
    long double plusmeanSumPt2[32]={0};
    long double plusmeanPt[32]={0};
    long double plusmeanPt2[32]={0};
    long double plusmeanMultPt[32]={0};
    long double plusmeanX2[32]={0};
    long double plusmeanX[32]={0};
    long double plusvarX[32]={0};
    long double plusptOmega[32]={0};
    long double plussumPtOmega[32]={0};

    long double plusmeanSumPtPrimary[32]={0};
    long double plusmeanSumPt2Primary[32]={0};
    long double plusmeanPtPrimary[32]={0};
    long double plusmeanPt2Primary[32]={0};
    long double plusmeanMultPtPrimary[32]={0};
    long double plusmeanX2Primary[32]={0};
    long double plusmeanXPrimary[32]={0};
    long double plusvarXPrimary[32]={0};
    long double plusptOmegaPrimary[32]={0};
    long double plussumPtOmegaPrimary[32]={0};


    long double plusmeanSumPtSim[32]={0};
    long double plusmeanSumPt2Sim[32]={0};
    long double plusmeanPtSim[32]={0};
    long double plusmeanPt2Sim[32]={0};
    long double plusmeanMultPtSim[32]={0};
    long double plusmeanX2Sim[32]={0};
    long double plusmeanXSim[32]={0};
    long double plusvarXSim[32]={0};
    long double plusptOmegaSim[32]={0};
    long double plussumPtOmegaSim[32]={0};





    long double minusmeanSumPt[32]={0};
    long double minusmeanSumPt2[32]={0};
    long double minusmeanPt[32]={0};
    long double minusmeanPt2[32]={0};
    long double minusmeanMultPt[32]={0};
    long double minusmeanX2[32]={0};
    long double minusmeanX[32]={0};
    long double minusvarX[32]={0};
    long double minusptOmega[32]={0};
    long double minussumPtOmega[32]={0};


    long double minusmeanSumPtPrimary[32]={0};
    long double minusmeanSumPt2Primary[32]={0};
    long double minusmeanPtPrimary[32]={0};
    long double minusmeanPt2Primary[32]={0};
    long double minusmeanMultPtPrimary[32]={0};
    long double minusmeanX2Primary[32]={0};
    long double minusmeanXPrimary[32]={0};
    long double minusvarXPrimary[32]={0};
    long double minusptOmegaPrimary[32]={0};
    long double minussumPtOmegaPrimary[32]={0};


    long double minusmeanSumPtSim[32]={0};
    long double minusmeanSumPt2Sim[32]={0};
    long double minusmeanPtSim[32]={0};
    long double minusmeanPt2Sim[32]={0};
    long double minusmeanMultPtSim[32]={0};
    long double minusmeanX2Sim[32]={0};
    long double minusmeanXSim[32]={0};
    long double minusvarXSim[32]={0};
    long double minusptOmegaSim[32]={0};
    long double minussumPtOmegaSim[32]={0};


    for (int j=nStartingEvent; j<nEndingEvent; j++) {
      extractedTree->GetEntry(j);

      for (Int_t iW=0; iW<nWindows; iW++) {
        forwMultSim[iW] = 0;
        backMultSim[iW] = 0;
        forwMultPlusSim[iW] = 0;
        backMultPlusSim[iW] = 0;
        forwMultMinusSim[iW] = 0;
        backMultMinusSim[iW] = 0;
      }

      for (Int_t jW=0; jW<nOneWindows; jW++) {
        totalMultSim[jW] = 0;
        totalPtSim[jW] = 0;
        totalPt2Sim[jW] = 0;

        plusMultSim[jW] = 0;
        plusPtSim[jW] = 0;
        plusPt2Sim[jW] = 0;

        minusMultSim[jW] = 0;
        minusPtSim[jW] = 0;
        minusPt2Sim[jW] = 0;
      }

      for (int i=0; i<f_SIM_NumberOfTracks; i++) {
        for (Int_t jW=0; jW<nOneWindows; jW++) {
          if (f_SIM_TrackEta[i]<0.05+jW*0.05 && f_SIM_TrackEta[i]>-0.05-jW*0.05) {
            totalMultSim[jW] = totalMultSim[jW] + 1;
            totalPtSim[jW] = totalPtSim[jW] + f_SIM_TrackPt[i];
            totalPt2Sim[jW] = totalPt2Sim[jW] + f_SIM_TrackPt[i]*f_SIM_TrackPt[i];
          }
        }

        for (Int_t iW=0; iW<nWindows; iW++) {
          if (f_SIM_TrackEta[i]<0.1+iW*0.1 && f_SIM_TrackEta[i]>0+iW*0.1) {
            forwMultSim[iW] = forwMultSim[iW] + 1;
          }
          if (f_SIM_TrackEta[i]<0-iW*0.1 && f_SIM_TrackEta[i]>-0.1-iW*0.1) {
            backMultSim[iW] = backMultSim[iW] + 1;
          }
        }

        if (f_SIM_TrackCharge[i]>0) {
          //                    plusMult++;
          //                    plusPt = plusPt + pt;
          //                    plusPt2 = plusPt2 + pt*pt;

          for (Int_t jW=0; jW<nOneWindows; jW++) {
            if (f_SIM_TrackEta[i]<0.05+jW*0.05 && f_SIM_TrackEta[i]>-0.05-jW*0.05) {
              plusMultSim[jW] = plusMultSim[jW] + 1;
              plusPtSim[jW] = plusPtSim[jW] + f_SIM_TrackPt[i];
              plusPt2Sim[jW] = plusPt2Sim[jW] + f_SIM_TrackPt[i]*f_SIM_TrackPt[i];
            }
          }

          for (Int_t iW=0; iW<nWindows; iW++) {
            if (f_SIM_TrackEta[i]<0.1+iW*0.1 && f_SIM_TrackEta[i]>0+iW*0.1) {
              forwMultPlusSim[iW] = forwMultPlusSim[iW] + 1;
            }
            if (f_SIM_TrackEta[i]<0-iW*0.1 && f_SIM_TrackEta[i]>-0.1-iW*0.1) {
              backMultPlusSim[iW] = backMultPlusSim[iW] + 1;
            }
          }
        } else {
          //                    minusMult++;
          //                    minusPt = minusPt + pt;
          //                    minusPt2 = minusPt2 + pt*pt;

          for (Int_t jW=0; jW<nOneWindows; jW++) {
            if (f_SIM_TrackEta[i]<0.05+jW*0.05 && f_SIM_TrackEta[i]>-0.05-jW*0.05) {
              minusMultSim[jW] = minusMultSim[jW] + 1;
              minusPtSim[jW] = minusPtSim[jW] + f_SIM_TrackPt[i];
              minusPt2Sim[jW] = minusPt2Sim[jW] + f_SIM_TrackPt[i]*f_SIM_TrackPt[i];
            }
          }

          for (Int_t iW=0; iW<nWindows; iW++) {
            if (f_SIM_TrackEta[i]<0.1+iW*0.1 && f_SIM_TrackEta[i]>0+iW*0.1) {
              forwMultMinusSim[iW] = forwMultMinusSim[iW] + 1;
            }
            if (f_SIM_TrackEta[i]<0-iW*0.1 && f_SIM_TrackEta[i]>-0.1-iW*0.1) {
              backMultMinusSim[iW] = backMultMinusSim[iW] + 1;
            }
          }
        }
      }



      for (Int_t iW=0; iW<nWindows; iW++) {
        forwMult[iW] = 0;
        backMult[iW] = 0;
        forwMultPlus[iW] = 0;
        backMultPlus[iW] = 0;
        forwMultMinus[iW] = 0;
        backMultMinus[iW] = 0;
      }

      for (Int_t jW=0; jW<nOneWindows; jW++) {
        totalMult[jW] = 0;
        totalPt[jW] = 0;
        totalPt2[jW] = 0;

        plusMult[jW] = 0;
        plusPt[jW] = 0;
        plusPt2[jW] = 0;

        minusMult[jW] = 0;
        minusPt[jW] = 0;
        minusPt2[jW] = 0;
      }



      for (Int_t iW=0; iW<nWindows; iW++) {
        forwMultPrimary[iW] = 0;
        backMultPrimary[iW] = 0;
        forwMultPlusPrimary[iW] = 0;
        backMultPlusPrimary[iW] = 0;
        forwMultMinusPrimary[iW] = 0;
        backMultMinusPrimary[iW] = 0;
      }

      for (Int_t jW=0; jW<nOneWindows; jW++) {
        totalMultPrimary[jW] = 0;
        totalPtPrimary[jW] = 0;
        totalPt2Primary[jW] = 0;

        plusMultPrimary[jW] = 0;
        plusPtPrimary[jW] = 0;
        plusPt2Primary[jW] = 0;

        minusMultPrimary[jW] = 0;
        minusPtPrimary[jW] = 0;
        minusPt2Primary[jW] = 0;
      }

      for (int i=0; i<fNumberOfTracks; i++) {
        for (Int_t jW=0; jW<nOneWindows; jW++) {
          if (fTrackEta[i]<0.05+jW*0.05 && fTrackEta[i]>-0.05-jW*0.05) {
            totalMult[jW] = totalMult[jW] + 1;
            totalPt[jW] = totalPt[jW] + fTrackPt[i];
            totalPt2[jW] = totalPt2[jW] + fTrackPt[i]*fTrackPt[i];
          }
        }

        for (Int_t iW=0; iW<nWindows; iW++) {
          if ((fTrackEta[i]<0.1+iW*0.1) && (fTrackEta[i]>0+iW*0.1)) {
            forwMult[iW] = forwMult[iW] + 1;
          }
          if ((fTrackEta[i]<0-iW*0.1) && (fTrackEta[i]>-0.1-iW*0.1)) {
            backMult[iW] = backMult[iW] + 1;
          }
        }




        for (Int_t jW=0; jW<nOneWindows; jW++) {
          if (fTrackEtaPrimary[i]<0.05+jW*0.05 && fTrackEtaPrimary[i]>-0.05-jW*0.05) {
            totalMultPrimary[jW] = totalMultPrimary[jW] + 1;
            totalPtPrimary[jW] = totalPtPrimary[jW] + fTrackPtPrimary[i];
            totalPt2Primary[jW] = totalPt2Primary[jW] + fTrackPtPrimary[i]*fTrackPtPrimary[i];
          }
        }

        for (Int_t iW=0; iW<nWindows; iW++) {
          if ((fTrackEtaPrimary[i]<0.1+iW*0.1) && (fTrackEtaPrimary[i]>0+iW*0.1)) {
            forwMultPrimary[iW] = forwMultPrimary[iW] + 1;
          }
          if ((fTrackEtaPrimary[i]<0-iW*0.1) && (fTrackEtaPrimary[i]>-0.1-iW*0.1)) {
            backMultPrimary[iW] = backMultPrimary[iW] + 1;
          }
        }






        if (fTrackCharge[i]>0) {

          for (Int_t jW=0; jW<nOneWindows; jW++) {
            if (fTrackEta[i]<0.05+jW*0.05 && fTrackEta[i]>-0.05-jW*0.05) {
              plusMult[jW] = plusMult[jW] + 1;
              plusPt[jW] = plusPt[jW] + fTrackPt[i];
              plusPt2[jW] = plusPt2[jW] + fTrackPt[i]*fTrackPt[i];
            }
          }

          for (Int_t iW=0; iW<nWindows; iW++) {
            if (fTrackEta[i]<0.1+iW*0.1 && fTrackEta[i]>0+iW*0.1) {
              forwMultPlus[iW] = forwMultPlus[iW] + 1;
            }
            if (fTrackEta[i]<0-iW*0.1 && fTrackEta[i]>-0.1-iW*0.1) {
              backMultPlus[iW] = backMultPlus[iW] + 1;
            }
          }



          for (Int_t jW=0; jW<nOneWindows; jW++) {
            if (fTrackEtaPrimary[i]<0.05+jW*0.05 && fTrackEtaPrimary[i]>-0.05-jW*0.05) {
              plusMultPrimary[jW] = plusMultPrimary[jW] + 1;
              plusPtPrimary[jW] = plusPtPrimary[jW] + fTrackPtPrimary[i];
              plusPt2Primary[jW] = plusPt2Primary[jW] + fTrackPtPrimary[i]*fTrackPtPrimary[i];
            }
          }

          for (Int_t iW=0; iW<nWindows; iW++) {
            if (fTrackEtaPrimary[i]<0.1+iW*0.1 && fTrackEtaPrimary[i]>0+iW*0.1) {
              forwMultPlusPrimary[iW] = forwMultPlusPrimary[iW] + 1;
            }
            if (fTrackEtaPrimary[i]<0-iW*0.1 && fTrackEtaPrimary[i]>-0.1-iW*0.1) {
              backMultPlusPrimary[iW] = backMultPlusPrimary[iW] + 1;
            }
          }





        } else {

          for (Int_t jW=0; jW<nOneWindows; jW++) {
            if (fTrackEta[i]<0.05+jW*0.05 && fTrackEta[i]>-0.05-jW*0.05) {
              minusMult[jW] = minusMult[jW] + 1;
              minusPt[jW] = minusPt[jW] + fTrackPt[i];
              minusPt2[jW] = minusPt2[jW] + fTrackPt[i]*fTrackPt[i];
            }
          }

          for (Int_t iW=0; iW<nWindows; iW++) {
            if (fTrackEta[i]<0.1+iW*0.1 && fTrackEta[i]>0+iW*0.1) {
              forwMultMinus[iW] = forwMultMinus[iW] + 1;
            }
            if (fTrackEta[i]<0-iW*0.1 && fTrackEta[i]>-0.1-iW*0.1) {
              backMultMinus[iW] = backMultMinus[iW] + 1;
            }
          }



          for (Int_t jW=0; jW<nOneWindows; jW++) {
            if (fTrackEtaPrimary[i]<0.05+jW*0.05 && fTrackEtaPrimary[i]>-0.05-jW*0.05) {
              minusMultPrimary[jW] = minusMultPrimary[jW] + 1;
              minusPtPrimary[jW] = minusPtPrimary[jW] + fTrackPtPrimary[i];
              minusPt2Primary[jW] = minusPt2Primary[jW] + fTrackPtPrimary[i]*fTrackPtPrimary[i];
            }
          }

          for (Int_t iW=0; iW<nWindows; iW++) {
            if (fTrackEtaPrimary[i]<0.1+iW*0.1 && fTrackEtaPrimary[i]>0+iW*0.1) {
              forwMultMinusPrimary[iW] = forwMultMinusPrimary[iW] + 1;
            }
            if (fTrackEtaPrimary[i]<0-iW*0.1 && fTrackEtaPrimary[i]>-0.1-iW*0.1) {
              backMultMinusPrimary[iW] = backMultMinusPrimary[iW] + 1;
            }
          }


        }

      }


      //implement event cuts below!!!

      double bForTreeLow, bForTreeHigh;
      bForTreeLow = -1;
      bForTreeHigh = 100;


      switch (lowCentrality) {
        case 0:
        bForTreeLow = -1;
        break;
        case 5:
        bForTreeLow = 3.23611;
        break;
        case 10:
        bForTreeLow = 4.59567;
        break;
        case 15:
        bForTreeLow = 5.63829;
        break;
        case 20:
        bForTreeLow = 6.51336;
        break;
        case 25:
        bForTreeLow = 7.28364;
        break;
        case 30:
        bForTreeLow = 7.97662;
        break;
        case 35:
        bForTreeLow = 8.61874;
        break;
        case 40:
        bForTreeLow = 9.21147;
        break;
        case 45:
        bForTreeLow = 9.76649;
        break;
        case 50:
        bForTreeLow = 10.2912;
        break;
        case 55:
        bForTreeLow = 10.7946;
        break;
        case 60:
        bForTreeLow = 11.2812;
        break;
        case 65:
        bForTreeLow = 11.742;
        break;
        case 70:
        bForTreeLow = 12.1956;
        break;
        case 75:
        bForTreeLow = 12.624;
        break;
        case 80:
        bForTreeLow = 13.0411;
        break;
        case 85:
        bForTreeLow = 13.466;
        break;
        case 90:
        bForTreeLow = 13.9349;
        break;
        case 95:
        bForTreeLow = 14.5656;
        break;
        default:
        bForTreeLow = -1;
        break;
      }




      switch (highCentrality) {
        case 0:
        cout << "don't do that, idiot";
        break;
        case 5:
        bForTreeHigh = 3.23611;
        break;
        case 10:
        bForTreeHigh = 4.59567;
        break;
        case 15:
        bForTreeHigh = 5.63829;
        break;
        case 20:
        bForTreeHigh = 6.51336;
        break;
        case 25:
        bForTreeHigh = 7.28364;
        break;
        case 30:
        bForTreeHigh = 7.97662;
        break;
        case 35:
        bForTreeHigh = 8.61874;
        break;
        case 40:
        bForTreeHigh = 9.21147;
        break;
        case 45:
        bForTreeHigh = 9.76649;
        break;
        case 50:
        bForTreeHigh = 10.2912;
        break;
        case 55:
        bForTreeHigh = 10.7946;
        break;
        case 60:
        bForTreeHigh = 11.2812;
        break;
        case 65:
        bForTreeHigh = 11.742;
        break;
        case 70:
        bForTreeHigh = 12.1956;
        break;
        case 75:
        bForTreeHigh = 12.624;
        break;
        case 80:
        bForTreeHigh = 13.0411;
        break;
        case 85:
        bForTreeHigh = 13.466;
        break;
        case 90:
        bForTreeHigh = 13.9349;
        break;
        case 95:
        bForTreeHigh = 14.5656;
        break;
        default:
        bForTreeHigh = 100;
        break;
      }










      if (impactPar>bForTreeHigh || impactPar<bForTreeLow) {
        continue;
      }







      nEventsInSubSample_D++;
      //event was analyzed, now starting to fill counters

      for (int i=0; i<32; i++) {


        fullSumPt[i] += totalPt[i];
        fullSumPt2[i] += totalPt2[i];
        fullSumPtpow2[i] += totalPt[i]*totalPt[i];
        fullSumMult[i] += totalMult[i];
        fullSumMultMult[i] += totalMult[i]*totalMult[i];
        fullSumMultPt[i] += totalMult[i]*totalPt[i];

        fullSumPtPrimary[i] += totalPtPrimary[i];
        fullSumPt2Primary[i] += totalPt2Primary[i];
        fullSumPtpow2Primary[i] += totalPtPrimary[i]*totalPtPrimary[i];
        fullSumMultPrimary[i] += totalMultPrimary[i];
        fullSumMultMultPrimary[i] += totalMultPrimary[i]*totalMultPrimary[i];
        fullSumMultPtPrimary[i] += totalMultPrimary[i]*totalPtPrimary[i];

        fullSumPtSim[i] += totalPtSim[i];
        fullSumPt2Sim[i] += totalPt2Sim[i];
        fullSumPtpow2Sim[i] += totalPtSim[i]*totalPtSim[i];
        fullSumMultSim[i] += totalMultSim[i];
        fullSumMultMultSim[i] += totalMultSim[i]*totalMultSim[i];
        fullSumMultPtSim[i] += totalMultSim[i]*totalPtSim[i];






        fullPlusSumPt[i] += plusPt[i];
        fullPlusSumPt2[i] += plusPt2[i];
        fullPlusSumPtpow2[i] += plusPt[i]*plusPt[i];
        fullPlusSumMult[i] += plusMult[i];
        fullPlusSumMultMult[i] += plusMult[i]*plusMult[i];
        fullPlusSumMultPt[i] += plusMult[i]*plusPt[i];

        fullPlusSumPtPrimary[i] += plusPtPrimary[i];
        fullPlusSumPt2Primary[i] += plusPt2Primary[i];
        fullPlusSumPtpow2Primary[i] += plusPtPrimary[i]*plusPtPrimary[i];
        fullPlusSumMultPrimary[i] += plusMultPrimary[i];
        fullPlusSumMultMultPrimary[i] += plusMultPrimary[i]*plusMultPrimary[i];
        fullPlusSumMultPtPrimary[i] += plusMultPrimary[i]*plusPtPrimary[i];

        fullPlusSumPtSim[i] += plusPtSim[i];
        fullPlusSumPt2Sim[i] += plusPt2Sim[i];
        fullPlusSumPtpow2Sim[i] += plusPtSim[i]*plusPtSim[i];
        fullPlusSumMultSim[i] += plusMultSim[i];
        fullPlusSumMultMultSim[i] += plusMultSim[i]*plusMultSim[i];
        fullPlusSumMultPtSim[i] += plusMultSim[i]*plusPtSim[i];





        fullMinusSumPt[i] += minusPt[i];
        fullMinusSumPt2[i] += minusPt2[i];
        fullMinusSumPtpow2[i] += minusPt[i]*minusPt[i];
        fullMinusSumMult[i] += minusMult[i];
        fullMinusSumMultMult[i] += minusMult[i]*minusMult[i];
        fullMinusSumMultPt[i] += minusMult[i]*minusPt[i];

        fullMinusSumPtPrimary[i] += minusPtPrimary[i];
        fullMinusSumPt2Primary[i] += minusPt2Primary[i];
        fullMinusSumPtpow2Primary[i] += minusPtPrimary[i]*minusPtPrimary[i];
        fullMinusSumMultPrimary[i] += minusMultPrimary[i];
        fullMinusSumMultMultPrimary[i] += minusMultPrimary[i]*minusMultPrimary[i];
        fullMinusSumMultPtPrimary[i] += minusMultPrimary[i]*minusPtPrimary[i];

        fullMinusSumPtSim[i] += minusPtSim[i];
        fullMinusSumPt2Sim[i] += minusPt2Sim[i];
        fullMinusSumPtpow2Sim[i] += minusPtSim[i]*minusPtSim[i];
        fullMinusSumMultSim[i] += minusMultSim[i];
        fullMinusSumMultMultSim[i] += minusMultSim[i]*minusMultSim[i];
        fullMinusSumMultPtSim[i] += minusMultSim[i]*minusPtSim[i];


      }


      for (int i=0; i<16; i++) {
        corrTotalHist[i][h]->Fill(forwMult[i],backMult[i]);

        corrPlusPlusHist[i][h]->Fill(forwMultPlus[i],backMultPlus[i]);
        corrMinusMinusHist[i][h]->Fill(forwMultMinus[i],backMultMinus[i]);
        corrPlusMinusHist[i][h]->Fill(forwMultPlus[i],backMultMinus[i]);
        corrMinusPlusHist[i][h]->Fill(forwMultMinus[i],backMultPlus[i]);


        corrTotalPrimaryHist[i][h]->Fill(forwMultPrimary[i],backMultPrimary[i]);

        corrPlusPlusPrimaryHist[i][h]->Fill(forwMultPlusPrimary[i],backMultPlusPrimary[i]);
        corrMinusMinusPrimaryHist[i][h]->Fill(forwMultMinusPrimary[i],backMultMinusPrimary[i]);
        corrPlusMinusPrimaryHist[i][h]->Fill(forwMultPlusPrimary[i],backMultMinusPrimary[i]);
        corrMinusPlusPrimaryHist[i][h]->Fill(forwMultMinusPrimary[i],backMultPlusPrimary[i]);


        corrTotalSimHist[i][h]->Fill(forwMultSim[i],backMultSim[i]);

        corrPlusPlusSimHist[i][h]->Fill(forwMultPlusSim[i],backMultPlusSim[i]);
        corrMinusMinusSimHist[i][h]->Fill(forwMultMinusSim[i],backMultMinusSim[i]);
        corrPlusMinusSimHist[i][h]->Fill(forwMultPlusSim[i],backMultMinusSim[i]);
        corrMinusPlusSimHist[i][h]->Fill(forwMultMinusSim[i],backMultPlusSim[i]);
      }



    }

    nStartingEvent = nEndingEvent;
    nEndingEvent = nEndingEvent + nEventsInSubSample;


    for (int i=0; i<32; i++) {
      meanSumPt[i] = fullSumPt[i] / nEventsInSubSample_D;
      meanSumPt2[i] = fullSumPtpow2[i] / nEventsInSubSample_D;
      meanPt[i] = fullSumPt[i] / fullSumMult[i];
      meanPt2[i] = fullSumPt2[i] / fullSumMult[i];
      meanMultPt[i] = fullSumMultPt[i] / nEventsInSubSample_D;

      meanX2[i] = fullSumMultMult[i] / nEventsInSubSample_D;
      meanX[i] = fullSumMult[i] / nEventsInSubSample_D;

      varX[i] = meanX2[i] - meanX[i]*meanX[i];

      ptOmega[i] = (meanPt2[i] - meanPt[i]*meanPt[i]) / meanPt[i];
      //cout << "ptOmega " << ptOmega << endl;
      sumPtOmega[i] = (meanSumPt2[i] - meanSumPt[i]*meanSumPt[i]) / meanSumPt[i];



      totalMultSub[i][h] = meanX[i];
      totalOmegaSub[i][h] = varX[i]/meanX[i];
      totalDeltaSub[i][h] = (meanX[i]*sumPtOmega[i] - meanSumPt[i]*totalOmegaSub[i][h]) / (meanX[i]*ptOmega[i]);
      totalSigmaSub[i][h] = (meanX[i]*sumPtOmega[i] + meanSumPt[i]*totalOmegaSub[i][h] - 2*(meanMultPt[i]-meanX[i]*meanSumPt[i])) / (meanX[i]*ptOmega[i]);






      meanSumPtPrimary[i] = fullSumPtPrimary[i] / nEventsInSubSample_D;
      meanSumPt2Primary[i] = fullSumPtpow2Primary[i] / nEventsInSubSample_D;
      meanPtPrimary[i] = fullSumPtPrimary[i] / fullSumMultPrimary[i];
      meanPt2Primary[i] = fullSumPt2Primary[i] / fullSumMultPrimary[i];
      meanMultPtPrimary[i] = fullSumMultPtPrimary[i] / nEventsInSubSample_D;

      meanX2Primary[i] = fullSumMultMultPrimary[i] / nEventsInSubSample_D;
      meanXPrimary[i] = fullSumMultPrimary[i] / nEventsInSubSample_D;

      varXPrimary[i] = meanX2Primary[i] - meanXPrimary[i]*meanXPrimary[i];

      ptOmegaPrimary[i] = (meanPt2Primary[i] - meanPtPrimary[i]*meanPtPrimary[i]) / meanPtPrimary[i];
      //cout << "ptOmega " << ptOmega << endl;
      sumPtOmegaPrimary[i] = (meanSumPt2Primary[i] - meanSumPtPrimary[i]*meanSumPtPrimary[i]) / meanSumPtPrimary[i];



      totalMultPrimarySub[i][h] = meanXPrimary[i];
      totalOmegaPrimarySub[i][h] = varXPrimary[i]/meanXPrimary[i];
      totalDeltaPrimarySub[i][h] = (meanXPrimary[i]*sumPtOmegaPrimary[i] - meanSumPtPrimary[i]*totalOmegaPrimarySub[i][h]) / (meanXPrimary[i]*ptOmegaPrimary[i]);
      totalSigmaPrimarySub[i][h] = (meanXPrimary[i]*sumPtOmegaPrimary[i] + meanSumPtPrimary[i]*totalOmegaPrimarySub[i][h] - 2*(meanMultPtPrimary[i]-meanXPrimary[i]*meanSumPtPrimary[i])) / (meanXPrimary[i]*ptOmegaPrimary[i]);






      meanSumPtSim[i] = fullSumPtSim[i] / nEventsInSubSample_D;
      meanSumPt2Sim[i] = fullSumPtpow2Sim[i] / nEventsInSubSample_D;
      meanPtSim[i] = fullSumPtSim[i] / fullSumMultSim[i];
      meanPt2Sim[i] = fullSumPt2Sim[i] / fullSumMultSim[i];
      meanMultPtSim[i] = fullSumMultPtSim[i] / nEventsInSubSample_D;

      meanX2Sim[i] = fullSumMultMultSim[i] / nEventsInSubSample_D;
      meanXSim[i] = fullSumMultSim[i] / nEventsInSubSample_D;

      varXSim[i] = meanX2Sim[i] - meanXSim[i]*meanXSim[i];

      ptOmegaSim[i] = (meanPt2Sim[i] - meanPtSim[i]*meanPtSim[i]) / meanPtSim[i];
      //cout << "ptOmega " << ptOmega << endl;
      sumPtOmegaSim[i] = (meanSumPt2Sim[i] - meanSumPtSim[i]*meanSumPtSim[i]) / meanSumPtSim[i];



      totalMultSimSub[i][h] = meanXSim[i];
      totalOmegaSimSub[i][h] = varXSim[i]/meanXSim[i];
      totalDeltaSimSub[i][h] = (meanXSim[i]*sumPtOmegaSim[i] - meanSumPtSim[i]*totalOmegaSimSub[i][h]) / (meanXSim[i]*ptOmegaSim[i]);
      totalSigmaSimSub[i][h] = (meanXSim[i]*sumPtOmegaSim[i] + meanSumPtSim[i]*totalOmegaSimSub[i][h] - 2*(meanMultPtSim[i]-meanXSim[i]*meanSumPtSim[i])) / (meanXSim[i]*ptOmegaSim[i]);

      //cout << "i = " << i << " ; h = " << h << " totalSigmaSub = " << totalSigmaSub[i][h] << " ; totalSigmaSimSub = " << totalSigmaSimSub[i][h] << endl;



      plusmeanSumPt[i] = fullPlusSumPt[i] / nEventsInSubSample_D;
      plusmeanSumPt2[i] = fullPlusSumPtpow2[i] / nEventsInSubSample_D;
      plusmeanPt[i] = fullPlusSumPt[i] / fullPlusSumMult[i];
      plusmeanPt2[i] = fullPlusSumPt2[i] / fullPlusSumMult[i];
      plusmeanMultPt[i] = fullPlusSumMultPt[i] / nEventsInSubSample_D;

      plusmeanX2[i] = fullPlusSumMultMult[i] / nEventsInSubSample_D;
      plusmeanX[i] = fullPlusSumMult[i] / nEventsInSubSample_D;

      plusvarX[i] = plusmeanX2[i] - plusmeanX[i]*plusmeanX[i];

      plusptOmega[i] = (plusmeanPt2[i] - plusmeanPt[i]*plusmeanPt[i]) / plusmeanPt[i];
      //cout << "ptOmega " << ptOmega << endl;
      plussumPtOmega[i] = (plusmeanSumPt2[i] - plusmeanSumPt[i]*plusmeanSumPt[i]) / plusmeanSumPt[i];



      plusMultSub[i][h] = plusmeanX[i];
      plusOmegaSub[i][h] = plusvarX[i]/plusmeanX[i];
      plusDeltaSub[i][h] = (plusmeanX[i]*plussumPtOmega[i] - plusmeanSumPt[i]*plusOmegaSub[i][h]) / (plusmeanX[i]*plusptOmega[i]);
      plusSigmaSub[i][h] = (plusmeanX[i]*plussumPtOmega[i] + plusmeanSumPt[i]*plusOmegaSub[i][h] - 2*(plusmeanMultPt[i]-plusmeanX[i]*plusmeanSumPt[i])) / (plusmeanX[i]*plusptOmega[i]);






      plusmeanSumPtPrimary[i] = fullPlusSumPtPrimary[i] / nEventsInSubSample_D;
      plusmeanSumPt2Primary[i] = fullPlusSumPtpow2Primary[i] / nEventsInSubSample_D;
      plusmeanPtPrimary[i] = fullPlusSumPtPrimary[i] / fullPlusSumMultPrimary[i];
      plusmeanPt2Primary[i] = fullPlusSumPt2Primary[i] / fullPlusSumMultPrimary[i];
      plusmeanMultPtPrimary[i] = fullPlusSumMultPtPrimary[i] / nEventsInSubSample_D;

      plusmeanX2Primary[i] = fullPlusSumMultMultPrimary[i] / nEventsInSubSample_D;
      plusmeanXPrimary[i] = fullPlusSumMultPrimary[i] / nEventsInSubSample_D;

      plusvarXPrimary[i] = plusmeanX2Primary[i] - plusmeanXPrimary[i]*plusmeanXPrimary[i];

      plusptOmegaPrimary[i] = (plusmeanPt2Primary[i] - plusmeanPtPrimary[i]*plusmeanPtPrimary[i]) / plusmeanPtPrimary[i];
      //cout << "ptOmega " << ptOmega << endl;
      plussumPtOmegaPrimary[i] = (plusmeanSumPt2Primary[i] - plusmeanSumPtPrimary[i]*plusmeanSumPtPrimary[i]) / plusmeanSumPtPrimary[i];



      plusMultPrimarySub[i][h] = plusmeanXPrimary[i];
      plusOmegaPrimarySub[i][h] = plusvarXPrimary[i]/plusmeanXPrimary[i];
      plusDeltaPrimarySub[i][h] = (plusmeanXPrimary[i]*plussumPtOmegaPrimary[i] - plusmeanSumPtPrimary[i]*plusOmegaPrimarySub[i][h]) / (plusmeanXPrimary[i]*plusptOmegaPrimary[i]);
      plusSigmaPrimarySub[i][h] = (plusmeanXPrimary[i]*plussumPtOmegaPrimary[i] + plusmeanSumPtPrimary[i]*plusOmegaPrimarySub[i][h] - 2*(plusmeanMultPtPrimary[i]-plusmeanXPrimary[i]*plusmeanSumPtPrimary[i])) / (plusmeanXPrimary[i]*plusptOmegaPrimary[i]);









      plusmeanSumPtSim[i] = fullPlusSumPtSim[i] / nEventsInSubSample_D;
      plusmeanSumPt2Sim[i] = fullPlusSumPtpow2Sim[i] / nEventsInSubSample_D;
      plusmeanPtSim[i] = fullPlusSumPtSim[i] / fullPlusSumMultSim[i];
      plusmeanPt2Sim[i] = fullPlusSumPt2Sim[i] / fullPlusSumMultSim[i];
      plusmeanMultPtSim[i] = fullPlusSumMultPtSim[i] / nEventsInSubSample_D;

      plusmeanX2Sim[i] = fullPlusSumMultMultSim[i] / nEventsInSubSample_D;
      plusmeanXSim[i] = fullPlusSumMultSim[i] / nEventsInSubSample_D;

      plusvarXSim[i] = plusmeanX2Sim[i] - plusmeanXSim[i]*plusmeanXSim[i];

      plusptOmegaSim[i] = (plusmeanPt2Sim[i] - plusmeanPtSim[i]*plusmeanPtSim[i]) / plusmeanPtSim[i];
      //cout << "ptOmega " << ptOmega << endl;
      plussumPtOmegaSim[i] = (plusmeanSumPt2Sim[i] - plusmeanSumPtSim[i]*plusmeanSumPtSim[i]) / plusmeanSumPtSim[i];



      plusMultSimSub[i][h] = plusmeanXSim[i];
      plusOmegaSimSub[i][h] = plusvarXSim[i]/plusmeanXSim[i];
      plusDeltaSimSub[i][h] = (plusmeanXSim[i]*plussumPtOmegaSim[i] - plusmeanSumPtSim[i]*plusOmegaSimSub[i][h]) / (plusmeanXSim[i]*plusptOmegaSim[i]);
      plusSigmaSimSub[i][h] = (plusmeanXSim[i]*plussumPtOmegaSim[i] + plusmeanSumPtSim[i]*plusOmegaSimSub[i][h] - 2*(plusmeanMultPtSim[i]-plusmeanXSim[i]*plusmeanSumPtSim[i])) / (plusmeanXSim[i]*plusptOmegaSim[i]);










      minusmeanSumPt[i] = fullMinusSumPt[i] / nEventsInSubSample_D;
      minusmeanSumPt2[i] = fullMinusSumPtpow2[i] / nEventsInSubSample_D;
      minusmeanPt[i] = fullMinusSumPt[i] / fullMinusSumMult[i];
      minusmeanPt2[i] = fullMinusSumPt2[i] / fullMinusSumMult[i];
      minusmeanMultPt[i] = fullMinusSumMultPt[i] / nEventsInSubSample_D;

      minusmeanX2[i] = fullMinusSumMultMult[i] / nEventsInSubSample_D;
      minusmeanX[i] = fullMinusSumMult[i] / nEventsInSubSample_D;

      minusvarX[i] = minusmeanX2[i] - minusmeanX[i]*minusmeanX[i];

      minusptOmega[i] = (minusmeanPt2[i] - minusmeanPt[i]*minusmeanPt[i]) / minusmeanPt[i];
      //cout << "ptOmega " << ptOmega << endl;
      minussumPtOmega[i] = (minusmeanSumPt2[i] - minusmeanSumPt[i]*minusmeanSumPt[i]) / minusmeanSumPt[i];



      minusMultSub[i][h] = minusmeanX[i];
      minusOmegaSub[i][h] = minusvarX[i]/minusmeanX[i];
      minusDeltaSub[i][h] = (minusmeanX[i]*minussumPtOmega[i] - minusmeanSumPt[i]*minusOmegaSub[i][h]) / (minusmeanX[i]*minusptOmega[i]);
      minusSigmaSub[i][h] = (minusmeanX[i]*minussumPtOmega[i] + minusmeanSumPt[i]*minusOmegaSub[i][h] - 2*(minusmeanMultPt[i]-minusmeanX[i]*minusmeanSumPt[i])) / (minusmeanX[i]*minusptOmega[i]);







      minusmeanSumPtPrimary[i] = fullMinusSumPtPrimary[i] / nEventsInSubSample_D;
      minusmeanSumPt2Primary[i] = fullMinusSumPtpow2Primary[i] / nEventsInSubSample_D;
      minusmeanPtPrimary[i] = fullMinusSumPtPrimary[i] / fullMinusSumMultPrimary[i];
      minusmeanPt2Primary[i] = fullMinusSumPt2Primary[i] / fullMinusSumMultPrimary[i];
      minusmeanMultPtPrimary[i] = fullMinusSumMultPtPrimary[i] / nEventsInSubSample_D;

      minusmeanX2Primary[i] = fullMinusSumMultMultPrimary[i] / nEventsInSubSample_D;
      minusmeanXPrimary[i] = fullMinusSumMultPrimary[i] / nEventsInSubSample_D;

      minusvarXPrimary[i] = minusmeanX2Primary[i] - minusmeanXPrimary[i]*minusmeanXPrimary[i];

      minusptOmegaPrimary[i] = (minusmeanPt2Primary[i] - minusmeanPtPrimary[i]*minusmeanPtPrimary[i]) / minusmeanPtPrimary[i];
      //cout << "ptOmega " << ptOmega << endl;
      minussumPtOmegaPrimary[i] = (minusmeanSumPt2Primary[i] - minusmeanSumPtPrimary[i]*minusmeanSumPtPrimary[i]) / minusmeanSumPtPrimary[i];



      minusMultPrimarySub[i][h] = minusmeanXPrimary[i];
      minusOmegaPrimarySub[i][h] = minusvarXPrimary[i]/minusmeanXPrimary[i];
      minusDeltaPrimarySub[i][h] = (minusmeanXPrimary[i]*minussumPtOmegaPrimary[i] - minusmeanSumPtPrimary[i]*minusOmegaPrimarySub[i][h]) / (minusmeanXPrimary[i]*minusptOmegaPrimary[i]);
      minusSigmaPrimarySub[i][h] = (minusmeanXPrimary[i]*minussumPtOmegaPrimary[i] + minusmeanSumPtPrimary[i]*minusOmegaPrimarySub[i][h] - 2*(minusmeanMultPtPrimary[i]-minusmeanXPrimary[i]*minusmeanSumPtPrimary[i])) / (minusmeanXPrimary[i]*minusptOmegaPrimary[i]);









      minusmeanSumPtSim[i] = fullMinusSumPtSim[i] / nEventsInSubSample_D;
      minusmeanSumPt2Sim[i] = fullMinusSumPtpow2Sim[i] / nEventsInSubSample_D;
      minusmeanPtSim[i] = fullMinusSumPtSim[i] / fullMinusSumMultSim[i];
      minusmeanPt2Sim[i] = fullMinusSumPt2Sim[i] / fullMinusSumMultSim[i];
      minusmeanMultPtSim[i] = fullMinusSumMultPtSim[i] / nEventsInSubSample_D;

      minusmeanX2Sim[i] = fullMinusSumMultMultSim[i] / nEventsInSubSample_D;
      minusmeanXSim[i] = fullMinusSumMultSim[i] / nEventsInSubSample_D;

      minusvarXSim[i] = minusmeanX2Sim[i] - minusmeanXSim[i]*minusmeanXSim[i];

      minusptOmegaSim[i] = (minusmeanPt2Sim[i] - minusmeanPtSim[i]*minusmeanPtSim[i]) / minusmeanPtSim[i];
      //cout << "ptOmega " << ptOmega << endl;
      minussumPtOmegaSim[i] = (minusmeanSumPt2Sim[i] - minusmeanSumPtSim[i]*minusmeanSumPtSim[i]) / minusmeanSumPtSim[i];



      minusMultSimSub[i][h] = minusmeanXSim[i];
      minusOmegaSimSub[i][h] = minusvarXSim[i]/minusmeanXSim[i];
      minusDeltaSimSub[i][h] = (minusmeanXSim[i]*minussumPtOmegaSim[i] - minusmeanSumPtSim[i]*minusOmegaSimSub[i][h]) / (minusmeanXSim[i]*minusptOmegaSim[i]);
      minusSigmaSimSub[i][h] = (minusmeanXSim[i]*minussumPtOmegaSim[i] + minusmeanSumPtSim[i]*minusOmegaSimSub[i][h] - 2*(minusmeanMultPtSim[i]-minusmeanXSim[i]*minusmeanSumPtSim[i])) / (minusmeanXSim[i]*minusptOmegaSim[i]);






    }

  }





  for (int i=0; i<32; i++) {
    for (Int_t h=0; h<nSubSamples; h++) {
      totalMultMean[i] = totalMultMean[i] + totalMultSub[i][h];
      totalOmegaMean[i] = totalOmegaMean[i] + totalOmegaSub[i][h];
      totalDeltaMean[i] = totalDeltaMean[i] + totalDeltaSub[i][h];
      totalSigmaMean[i] = totalSigmaMean[i] + totalSigmaSub[i][h];

      totalMultPrimaryMean[i] = totalMultPrimaryMean[i] + totalMultPrimarySub[i][h];
      totalOmegaPrimaryMean[i] = totalOmegaPrimaryMean[i] + totalOmegaPrimarySub[i][h];
      totalDeltaPrimaryMean[i] = totalDeltaPrimaryMean[i] + totalDeltaPrimarySub[i][h];
      totalSigmaPrimaryMean[i] = totalSigmaPrimaryMean[i] + totalSigmaPrimarySub[i][h];

      totalMultSimMean[i] = totalMultSimMean[i] + totalMultSimSub[i][h];
      totalOmegaSimMean[i] = totalOmegaSimMean[i] + totalOmegaSimSub[i][h];
      totalDeltaSimMean[i] = totalDeltaSimMean[i] + totalDeltaSimSub[i][h];
      totalSigmaSimMean[i] = totalSigmaSimMean[i] + totalSigmaSimSub[i][h];


      plusMultMean[i] = plusMultMean[i] + plusMultSub[i][h];
      plusOmegaMean[i] = plusOmegaMean[i] + plusOmegaSub[i][h];
      plusDeltaMean[i] = plusDeltaMean[i] + plusDeltaSub[i][h];
      plusSigmaMean[i] = plusSigmaMean[i] + plusSigmaSub[i][h];

      plusMultPrimaryMean[i] = plusMultPrimaryMean[i] + plusMultPrimarySub[i][h];
      plusOmegaPrimaryMean[i] = plusOmegaPrimaryMean[i] + plusOmegaPrimarySub[i][h];
      plusDeltaPrimaryMean[i] = plusDeltaPrimaryMean[i] + plusDeltaPrimarySub[i][h];
      plusSigmaPrimaryMean[i] = plusSigmaPrimaryMean[i] + plusSigmaPrimarySub[i][h];

      plusMultSimMean[i] = plusMultSimMean[i] + plusMultSimSub[i][h];
      plusOmegaSimMean[i] = plusOmegaSimMean[i] + plusOmegaSimSub[i][h];
      plusDeltaSimMean[i] = plusDeltaSimMean[i] + plusDeltaSimSub[i][h];
      plusSigmaSimMean[i] = plusSigmaSimMean[i] + plusSigmaSimSub[i][h];


      minusMultMean[i] = minusMultMean[i] + minusMultSub[i][h];
      minusOmegaMean[i] = minusOmegaMean[i] + minusOmegaSub[i][h];
      minusDeltaMean[i] = minusDeltaMean[i] + minusDeltaSub[i][h];
      minusSigmaMean[i] = minusSigmaMean[i] + minusSigmaSub[i][h];

      minusMultPrimaryMean[i] = minusMultPrimaryMean[i] + minusMultPrimarySub[i][h];
      minusOmegaPrimaryMean[i] = minusOmegaPrimaryMean[i] + minusOmegaPrimarySub[i][h];
      minusDeltaPrimaryMean[i] = minusDeltaPrimaryMean[i] + minusDeltaPrimarySub[i][h];
      minusSigmaPrimaryMean[i] = minusSigmaPrimaryMean[i] + minusSigmaPrimarySub[i][h];

      minusMultSimMean[i] = minusMultSimMean[i] + minusMultSimSub[i][h];
      minusOmegaSimMean[i] = minusOmegaSimMean[i] + minusOmegaSimSub[i][h];
      minusDeltaSimMean[i] = minusDeltaSimMean[i] + minusDeltaSimSub[i][h];
      minusSigmaSimMean[i] = minusSigmaSimMean[i] + minusSigmaSimSub[i][h];
    }

    totalMultMean[i] = totalMultMean[i]/nSubSamplesDbl;
    totalOmegaMean[i] = totalOmegaMean[i]/nSubSamplesDbl;
    totalDeltaMean[i] = totalDeltaMean[i]/nSubSamplesDbl;
    totalSigmaMean[i] = totalSigmaMean[i]/nSubSamplesDbl;

    totalMultPrimaryMean[i] = totalMultPrimaryMean[i]/nSubSamplesDbl;
    totalOmegaPrimaryMean[i] = totalOmegaPrimaryMean[i]/nSubSamplesDbl;
    totalDeltaPrimaryMean[i] = totalDeltaPrimaryMean[i]/nSubSamplesDbl;
    totalSigmaPrimaryMean[i] = totalSigmaPrimaryMean[i]/nSubSamplesDbl;

    totalMultSimMean[i] = totalMultSimMean[i]/nSubSamplesDbl;
    totalOmegaSimMean[i] = totalOmegaSimMean[i]/nSubSamplesDbl;
    totalDeltaSimMean[i] = totalDeltaSimMean[i]/nSubSamplesDbl;
    totalSigmaSimMean[i] = totalSigmaSimMean[i]/nSubSamplesDbl;


    plusMultMean[i] = plusMultMean[i]/nSubSamplesDbl;
    plusOmegaMean[i] = plusOmegaMean[i]/nSubSamplesDbl;
    plusDeltaMean[i] = plusDeltaMean[i]/nSubSamplesDbl;
    plusSigmaMean[i] = plusSigmaMean[i]/nSubSamplesDbl;

    plusMultPrimaryMean[i] = plusMultPrimaryMean[i]/nSubSamplesDbl;
    plusOmegaPrimaryMean[i] = plusOmegaPrimaryMean[i]/nSubSamplesDbl;
    plusDeltaPrimaryMean[i] = plusDeltaPrimaryMean[i]/nSubSamplesDbl;
    plusSigmaPrimaryMean[i] = plusSigmaPrimaryMean[i]/nSubSamplesDbl;

    plusMultSimMean[i] = plusMultSimMean[i]/nSubSamplesDbl;
    plusOmegaSimMean[i] = plusOmegaSimMean[i]/nSubSamplesDbl;
    plusDeltaSimMean[i] = plusDeltaSimMean[i]/nSubSamplesDbl;
    plusSigmaSimMean[i] = plusSigmaSimMean[i]/nSubSamplesDbl;



    minusMultMean[i] = minusMultMean[i]/nSubSamplesDbl;
    minusOmegaMean[i] = minusOmegaMean[i]/nSubSamplesDbl;
    minusDeltaMean[i] = minusDeltaMean[i]/nSubSamplesDbl;
    minusSigmaMean[i] = minusSigmaMean[i]/nSubSamplesDbl;

    minusMultPrimaryMean[i] = minusMultPrimaryMean[i]/nSubSamplesDbl;
    minusOmegaPrimaryMean[i] = minusOmegaPrimaryMean[i]/nSubSamplesDbl;
    minusDeltaPrimaryMean[i] = minusDeltaPrimaryMean[i]/nSubSamplesDbl;
    minusSigmaPrimaryMean[i] = minusSigmaPrimaryMean[i]/nSubSamplesDbl;

    minusMultSimMean[i] = minusMultSimMean[i]/nSubSamplesDbl;
    minusOmegaSimMean[i] = minusOmegaSimMean[i]/nSubSamplesDbl;
    minusDeltaSimMean[i] = minusDeltaSimMean[i]/nSubSamplesDbl;
    minusSigmaSimMean[i] = minusSigmaSimMean[i]/nSubSamplesDbl;



    for (Int_t h=0; h<nSubSamples; h++) {
      totalMultError[i] += (totalMultMean[i]-totalMultSub[i][h])*(totalMultMean[i]-totalMultSub[i][h]);
      totalOmegaError[i] += (totalOmegaMean[i]-totalOmegaSub[i][h])*(totalOmegaMean[i]-totalOmegaSub[i][h]);
      totalDeltaError[i] += (totalDeltaMean[i]-totalDeltaSub[i][h])*(totalDeltaMean[i]-totalDeltaSub[i][h]);
      totalSigmaError[i] += (totalSigmaMean[i]-totalSigmaSub[i][h])*(totalSigmaMean[i]-totalSigmaSub[i][h]);

      totalMultPrimaryError[i] += (totalMultPrimaryMean[i]-totalMultPrimarySub[i][h])*(totalMultPrimaryMean[i]-totalMultPrimarySub[i][h]);
      totalOmegaPrimaryError[i] += (totalOmegaPrimaryMean[i]-totalOmegaPrimarySub[i][h])*(totalOmegaPrimaryMean[i]-totalOmegaPrimarySub[i][h]);
      totalDeltaPrimaryError[i] += (totalDeltaPrimaryMean[i]-totalDeltaPrimarySub[i][h])*(totalDeltaPrimaryMean[i]-totalDeltaPrimarySub[i][h]);
      totalSigmaPrimaryError[i] += (totalSigmaPrimaryMean[i]-totalSigmaPrimarySub[i][h])*(totalSigmaPrimaryMean[i]-totalSigmaPrimarySub[i][h]);

      totalMultSimError[i] += (totalMultSimMean[i]-totalMultSimSub[i][h])*(totalMultSimMean[i]-totalMultSimSub[i][h]);
      totalOmegaSimError[i] += (totalOmegaSimMean[i]-totalOmegaSimSub[i][h])*(totalOmegaSimMean[i]-totalOmegaSimSub[i][h]);
      totalDeltaSimError[i] += (totalDeltaSimMean[i]-totalDeltaSimSub[i][h])*(totalDeltaSimMean[i]-totalDeltaSimSub[i][h]);
      totalSigmaSimError[i] += (totalSigmaSimMean[i]-totalSigmaSimSub[i][h])*(totalSigmaSimMean[i]-totalSigmaSimSub[i][h]);



      plusMultError[i] += (plusMultMean[i]-plusMultSub[i][h])*(plusMultMean[i]-plusMultSub[i][h]);
      plusOmegaError[i] += (plusOmegaMean[i]-plusOmegaSub[i][h])*(plusOmegaMean[i]-plusOmegaSub[i][h]);
      plusDeltaError[i] += (plusDeltaMean[i]-plusDeltaSub[i][h])*(plusDeltaMean[i]-plusDeltaSub[i][h]);
      plusSigmaError[i] += (plusSigmaMean[i]-plusSigmaSub[i][h])*(plusSigmaMean[i]-plusSigmaSub[i][h]);

      plusMultPrimaryError[i] += (plusMultPrimaryMean[i]-plusMultPrimarySub[i][h])*(plusMultPrimaryMean[i]-plusMultPrimarySub[i][h]);
      plusOmegaPrimaryError[i] += (plusOmegaPrimaryMean[i]-plusOmegaPrimarySub[i][h])*(plusOmegaPrimaryMean[i]-plusOmegaPrimarySub[i][h]);
      plusDeltaPrimaryError[i] += (plusDeltaPrimaryMean[i]-plusDeltaPrimarySub[i][h])*(plusDeltaPrimaryMean[i]-plusDeltaPrimarySub[i][h]);
      plusSigmaPrimaryError[i] += (plusSigmaPrimaryMean[i]-plusSigmaPrimarySub[i][h])*(plusSigmaPrimaryMean[i]-plusSigmaPrimarySub[i][h]);

      plusMultSimError[i] += (plusMultSimMean[i]-plusMultSimSub[i][h])*(plusMultSimMean[i]-plusMultSimSub[i][h]);
      plusOmegaSimError[i] += (plusOmegaSimMean[i]-plusOmegaSimSub[i][h])*(plusOmegaSimMean[i]-plusOmegaSimSub[i][h]);
      plusDeltaSimError[i] += (plusDeltaSimMean[i]-plusDeltaSimSub[i][h])*(plusDeltaSimMean[i]-plusDeltaSimSub[i][h]);
      plusSigmaSimError[i] += (plusSigmaSimMean[i]-plusSigmaSimSub[i][h])*(plusSigmaSimMean[i]-plusSigmaSimSub[i][h]);



      minusMultError[i] += (minusMultMean[i]-minusMultSub[i][h])*(minusMultMean[i]-minusMultSub[i][h]);
      minusOmegaError[i] += (minusOmegaMean[i]-minusOmegaSub[i][h])*(minusOmegaMean[i]-minusOmegaSub[i][h]);
      minusDeltaError[i] += (minusDeltaMean[i]-minusDeltaSub[i][h])*(minusDeltaMean[i]-minusDeltaSub[i][h]);
      minusSigmaError[i] += (minusSigmaMean[i]-minusSigmaSub[i][h])*(minusSigmaMean[i]-minusSigmaSub[i][h]);

      minusMultPrimaryError[i] += (minusMultPrimaryMean[i]-minusMultPrimarySub[i][h])*(minusMultPrimaryMean[i]-minusMultPrimarySub[i][h]);
      minusOmegaPrimaryError[i] += (minusOmegaPrimaryMean[i]-minusOmegaPrimarySub[i][h])*(minusOmegaPrimaryMean[i]-minusOmegaPrimarySub[i][h]);
      minusDeltaPrimaryError[i] += (minusDeltaPrimaryMean[i]-minusDeltaPrimarySub[i][h])*(minusDeltaPrimaryMean[i]-minusDeltaPrimarySub[i][h]);
      minusSigmaPrimaryError[i] += (minusSigmaPrimaryMean[i]-minusSigmaPrimarySub[i][h])*(minusSigmaPrimaryMean[i]-minusSigmaPrimarySub[i][h]);

      minusMultSimError[i] += (minusMultSimMean[i]-minusMultSimSub[i][h])*(minusMultSimMean[i]-minusMultSimSub[i][h]);
      minusOmegaSimError[i] += (minusOmegaSimMean[i]-minusOmegaSimSub[i][h])*(minusOmegaSimMean[i]-minusOmegaSimSub[i][h]);
      minusDeltaSimError[i] += (minusDeltaSimMean[i]-minusDeltaSimSub[i][h])*(minusDeltaSimMean[i]-minusDeltaSimSub[i][h]);
      minusSigmaSimError[i] += (minusSigmaSimMean[i]-minusSigmaSimSub[i][h])*(minusSigmaSimMean[i]-minusSigmaSimSub[i][h]);


    }

    totalMultError[i] = sqrt(totalMultError[i])/dtmp2;
    totalOmegaError[i] = sqrt(totalOmegaError[i])/dtmp2;
    totalDeltaError[i] = sqrt(totalDeltaError[i])/dtmp2;
    totalSigmaError[i] = sqrt(totalSigmaError[i])/dtmp2;

    totalMultPrimaryError[i] = sqrt(totalMultPrimaryError[i])/dtmp2;
    totalOmegaPrimaryError[i] = sqrt(totalOmegaPrimaryError[i])/dtmp2;
    totalDeltaPrimaryError[i] = sqrt(totalDeltaPrimaryError[i])/dtmp2;
    totalSigmaPrimaryError[i] = sqrt(totalSigmaPrimaryError[i])/dtmp2;

    totalMultSimError[i] = sqrt(totalMultSimError[i])/dtmp2;
    totalOmegaSimError[i] = sqrt(totalOmegaSimError[i])/dtmp2;
    totalDeltaSimError[i] = sqrt(totalDeltaSimError[i])/dtmp2;
    totalSigmaSimError[i] = sqrt(totalSigmaSimError[i])/dtmp2;



    plusMultError[i] = sqrt(plusMultError[i])/dtmp2;
    plusOmegaError[i] = sqrt(plusOmegaError[i])/dtmp2;
    plusDeltaError[i] = sqrt(plusDeltaError[i])/dtmp2;
    plusSigmaError[i] = sqrt(plusSigmaError[i])/dtmp2;

    plusMultPrimaryError[i] = sqrt(plusMultPrimaryError[i])/dtmp2;
    plusOmegaPrimaryError[i] = sqrt(plusOmegaPrimaryError[i])/dtmp2;
    plusDeltaPrimaryError[i] = sqrt(plusDeltaPrimaryError[i])/dtmp2;
    plusSigmaPrimaryError[i] = sqrt(plusSigmaPrimaryError[i])/dtmp2;

    plusMultSimError[i] = sqrt(plusMultSimError[i])/dtmp2;
    plusOmegaSimError[i] = sqrt(plusOmegaSimError[i])/dtmp2;
    plusDeltaSimError[i] = sqrt(plusDeltaSimError[i])/dtmp2;
    plusSigmaSimError[i] = sqrt(plusSigmaSimError[i])/dtmp2;




    minusMultError[i] = sqrt(minusMultError[i])/dtmp2;
    minusOmegaError[i] = sqrt(minusOmegaError[i])/dtmp2;
    minusDeltaError[i] = sqrt(minusDeltaError[i])/dtmp2;
    minusSigmaError[i] = sqrt(minusSigmaError[i])/dtmp2;

    minusMultPrimaryError[i] = sqrt(minusMultPrimaryError[i])/dtmp2;
    minusOmegaPrimaryError[i] = sqrt(minusOmegaPrimaryError[i])/dtmp2;
    minusDeltaPrimaryError[i] = sqrt(minusDeltaPrimaryError[i])/dtmp2;
    minusSigmaPrimaryError[i] = sqrt(minusSigmaPrimaryError[i])/dtmp2;

    minusMultSimError[i] = sqrt(minusMultSimError[i])/dtmp2;
    minusOmegaSimError[i] = sqrt(minusOmegaSimError[i])/dtmp2;
    minusDeltaSimError[i] = sqrt(minusDeltaSimError[i])/dtmp2;
    minusSigmaSimError[i] = sqrt(minusSigmaSimError[i])/dtmp2;




    cout << "Delta eta = " << 0.1+i*0.1 << endl;
    cout << "totalMultMean[" << i << "] = " << totalMultMean[i] << " ; +/- " << totalMultError[i] << endl;
    cout << "totalOmegaMean[" << i << "] = " << totalOmegaMean[i] << " ; +/- " << totalOmegaError[i] << endl;
    cout << "totalDeltaMean[" << i << "] = " << totalDeltaMean[i] << " ; +/- " << totalDeltaError[i] << endl;
    cout << "totalSigmaMean[" << i << "] = " << totalSigmaMean[i] << " ; +/- " << totalSigmaError[i] << endl;

    cout << "totalMultPrimaryMean[" << i << "] = " << totalMultPrimaryMean[i] << " ; +/- " << totalMultPrimaryError[i] << endl;
    cout << "totalOmegaPrimaryMean[" << i << "] = " << totalOmegaPrimaryMean[i] << " ; +/- " << totalOmegaPrimaryError[i] << endl;
    cout << "totalDeltaPrimaryMean[" << i << "] = " << totalDeltaPrimaryMean[i] << " ; +/- " << totalDeltaPrimaryError[i] << endl;
    cout << "totalSigmaPrimaryMean[" << i << "] = " << totalSigmaPrimaryMean[i] << " ; +/- " << totalSigmaPrimaryError[i] << endl;

    cout << "totalMultSimMean[" << i << "] = " << totalMultSimMean[i] << " ; +/- " << totalMultSimError[i] << endl;
    cout << "totalOmegaSimMean[" << i << "] = " << totalOmegaSimMean[i] << " ; +/- " << totalOmegaSimError[i] << endl;
    cout << "totalDeltaSimMean[" << i << "] = " << totalDeltaSimMean[i] << " ; +/- " << totalDeltaSimError[i] << endl;
    cout << "totalSigmaSimMean[" << i << "] = " << totalSigmaSimMean[i] << " ; +/- " << totalSigmaSimError[i] << endl;


    cout << "plusMultMean[" << i << "] = " << plusMultMean[i] << " ; +/- " << plusMultError[i] << endl;
    cout << "plusOmegaMean[" << i << "] = " << plusOmegaMean[i] << " ; +/- " << plusOmegaError[i] << endl;
    cout << "plusDeltaMean[" << i << "] = " << plusDeltaMean[i] << " ; +/- " << plusDeltaError[i] << endl;
    cout << "plusSigmaMean[" << i << "] = " << plusSigmaMean[i] << " ; +/- " << plusSigmaError[i] << endl;

    cout << "plusMultPrimaryMean[" << i << "] = " << plusMultPrimaryMean[i] << " ; +/- " << plusMultPrimaryError[i] << endl;
    cout << "plusOmegaPrimaryMean[" << i << "] = " << plusOmegaPrimaryMean[i] << " ; +/- " << plusOmegaPrimaryError[i] << endl;
    cout << "plusDeltaPrimaryMean[" << i << "] = " << plusDeltaPrimaryMean[i] << " ; +/- " << plusDeltaPrimaryError[i] << endl;
    cout << "plusSigmaPrimaryMean[" << i << "] = " << plusSigmaPrimaryMean[i] << " ; +/- " << plusSigmaPrimaryError[i] << endl;

    cout << "plusMultSimMean[" << i << "] = " << plusMultSimMean[i] << " ; +/- " << plusMultSimError[i] << endl;
    cout << "plusOmegaSimMean[" << i << "] = " << plusOmegaSimMean[i] << " ; +/- " << plusOmegaSimError[i] << endl;
    cout << "plusDeltaSimMean[" << i << "] = " << plusDeltaSimMean[i] << " ; +/- " << plusDeltaSimError[i] << endl;
    cout << "plusSigmaSimMean[" << i << "] = " << plusSigmaSimMean[i] << " ; +/- " << plusSigmaSimError[i] << endl;



    cout << "minusMultMean[" << i << "] = " << minusMultMean[i] << " ; +/- " << minusMultError[i] << endl;
    cout << "minusOmegaMean[" << i << "] = " << minusOmegaMean[i] << " ; +/- " << minusOmegaError[i] << endl;
    cout << "minusDeltaMean[" << i << "] = " << minusDeltaMean[i] << " ; +/- " << minusDeltaError[i] << endl;
    cout << "minusSigmaMean[" << i << "] = " << minusSigmaMean[i] << " ; +/- " << minusSigmaError[i] << endl;

    cout << "minusMultPrimaryMean[" << i << "] = " << minusMultPrimaryMean[i] << " ; +/- " << minusMultPrimaryError[i] << endl;
    cout << "minusOmegaPrimaryMean[" << i << "] = " << minusOmegaPrimaryMean[i] << " ; +/- " << minusOmegaPrimaryError[i] << endl;
    cout << "minusDeltaPrimaryMean[" << i << "] = " << minusDeltaPrimaryMean[i] << " ; +/- " << minusDeltaPrimaryError[i] << endl;
    cout << "minusSigmaPrimaryMean[" << i << "] = " << minusSigmaPrimaryMean[i] << " ; +/- " << minusSigmaPrimaryError[i] << endl;

    cout << "minusMultSimMean[" << i << "] = " << minusMultSimMean[i] << " ; +/- " << minusMultSimError[i] << endl;
    cout << "minusOmegaSimMean[" << i << "] = " << minusOmegaSimMean[i] << " ; +/- " << minusOmegaSimError[i] << endl;
    cout << "minusDeltaSimMean[" << i << "] = " << minusDeltaSimMean[i] << " ; +/- " << minusDeltaSimError[i] << endl;
    cout << "minusSigmaSimMean[" << i << "] = " << minusSigmaSimMean[i] << " ; +/- " << minusSigmaSimError[i] << endl;


  }

  TGraphErrors* graphMeanTotal = new TGraphErrors(24);
  TGraphErrors* graphMeanTotalPrimary = new TGraphErrors(24);
  TGraphErrors* graphMeanTotalSim = new TGraphErrors(24);
  TGraphErrors* graphOmegaTotal = new TGraphErrors(24);
  TGraphErrors* graphOmegaTotalPrimary = new TGraphErrors(24);
  TGraphErrors* graphOmegaTotalSim = new TGraphErrors(24);
  TGraphErrors* graphDeltaTotal = new TGraphErrors(24);
  TGraphErrors* graphDeltaTotalPrimary = new TGraphErrors(24);
  TGraphErrors* graphDeltaTotalSim = new TGraphErrors(24);
  TGraphErrors* graphSigmaTotal = new TGraphErrors(24);
  TGraphErrors* graphSigmaTotalPrimary = new TGraphErrors(24);
  TGraphErrors* graphSigmaTotalSim = new TGraphErrors(24);


  TGraphErrors* graphMeanPlus = new TGraphErrors(24);
  TGraphErrors* graphMeanPlusPrimary = new TGraphErrors(24);
  TGraphErrors* graphMeanPlusSim = new TGraphErrors(24);
  TGraphErrors* graphOmegaPlus = new TGraphErrors(24);
  TGraphErrors* graphOmegaPlusPrimary = new TGraphErrors(24);
  TGraphErrors* graphOmegaPlusSim = new TGraphErrors(24);
  TGraphErrors* graphDeltaPlus = new TGraphErrors(24);
  TGraphErrors* graphDeltaPlusPrimary = new TGraphErrors(24);
  TGraphErrors* graphDeltaPlusSim = new TGraphErrors(24);
  TGraphErrors* graphSigmaPlus = new TGraphErrors(24);
  TGraphErrors* graphSigmaPlusPrimary = new TGraphErrors(24);
  TGraphErrors* graphSigmaPlusSim = new TGraphErrors(24);


  TGraphErrors* graphMeanMinus = new TGraphErrors(24);
  TGraphErrors* graphMeanMinusPrimary = new TGraphErrors(24);
  TGraphErrors* graphMeanMinusSim = new TGraphErrors(24);
  TGraphErrors* graphOmegaMinus = new TGraphErrors(24);
  TGraphErrors* graphOmegaMinusPrimary = new TGraphErrors(24);
  TGraphErrors* graphOmegaMinusSim = new TGraphErrors(24);
  TGraphErrors* graphDeltaMinus = new TGraphErrors(24);
  TGraphErrors* graphDeltaMinusPrimary = new TGraphErrors(24);
  TGraphErrors* graphDeltaMinusSim = new TGraphErrors(24);
  TGraphErrors* graphSigmaMinus = new TGraphErrors(24);
  TGraphErrors* graphSigmaMinusPrimary = new TGraphErrors(24);
  TGraphErrors* graphSigmaMinusSim = new TGraphErrors(24);


  for (int i=0; i<24; i++) {
    graphMeanTotal->SetPoint(i,0.1+i*0.1,totalMultMean[i]);
    graphMeanTotal->SetPointError(i,0,totalMultError[i]);
    graphOmegaTotal->SetPoint(i,0.1+i*0.1,totalOmegaMean[i]);
    graphOmegaTotal->SetPointError(i,0,totalOmegaError[i]);
    graphDeltaTotal->SetPoint(i,0.1+i*0.1,totalDeltaMean[i]);
    graphDeltaTotal->SetPointError(i,0,totalDeltaError[i]);
    graphSigmaTotal->SetPoint(i,0.1+i*0.1,totalSigmaMean[i]);
    graphSigmaTotal->SetPointError(i,0,totalSigmaError[i]);

    graphMeanTotalPrimary->SetPoint(i,0.1+i*0.1,totalMultPrimaryMean[i]);
    graphMeanTotalPrimary->SetPointError(i,0,totalMultPrimaryError[i]);
    graphOmegaTotalPrimary->SetPoint(i,0.1+i*0.1,totalOmegaPrimaryMean[i]);
    graphOmegaTotalPrimary->SetPointError(i,0,totalOmegaPrimaryError[i]);
    graphDeltaTotalPrimary->SetPoint(i,0.1+i*0.1,totalDeltaPrimaryMean[i]);
    graphDeltaTotalPrimary->SetPointError(i,0,totalDeltaPrimaryError[i]);
    graphSigmaTotalPrimary->SetPoint(i,0.1+i*0.1,totalSigmaPrimaryMean[i]);
    graphSigmaTotalPrimary->SetPointError(i,0,totalSigmaPrimaryError[i]);

    graphMeanTotalSim->SetPoint(i,0.1+i*0.1,totalMultSimMean[i]);
    graphMeanTotalSim->SetPointError(i,0,totalMultSimError[i]);
    graphOmegaTotalSim->SetPoint(i,0.1+i*0.1,totalOmegaSimMean[i]);
    graphOmegaTotalSim->SetPointError(i,0,totalOmegaSimError[i]);
    graphDeltaTotalSim->SetPoint(i,0.1+i*0.1,totalDeltaSimMean[i]);
    graphDeltaTotalSim->SetPointError(i,0,totalDeltaSimError[i]);
    graphSigmaTotalSim->SetPoint(i,0.1+i*0.1,totalSigmaSimMean[i]);
    graphSigmaTotalSim->SetPointError(i,0,totalSigmaSimError[i]);

    graphMeanPlus->SetPoint(i,0.1+i*0.1,plusMultMean[i]);
    graphMeanPlus->SetPointError(i,0,plusMultError[i]);
    graphOmegaPlus->SetPoint(i,0.1+i*0.1,plusOmegaMean[i]);
    graphOmegaPlus->SetPointError(i,0,plusOmegaError[i]);
    graphDeltaPlus->SetPoint(i,0.1+i*0.1,plusDeltaMean[i]);
    graphDeltaPlus->SetPointError(i,0,plusDeltaError[i]);
    graphSigmaPlus->SetPoint(i,0.1+i*0.1,plusSigmaMean[i]);
    graphSigmaPlus->SetPointError(i,0,plusSigmaError[i]);

    graphMeanPlusPrimary->SetPoint(i,0.1+i*0.1,plusMultPrimaryMean[i]);
    graphMeanPlusPrimary->SetPointError(i,0,plusMultPrimaryError[i]);
    graphOmegaPlusPrimary->SetPoint(i,0.1+i*0.1,plusOmegaPrimaryMean[i]);
    graphOmegaPlusPrimary->SetPointError(i,0,plusOmegaPrimaryError[i]);
    graphDeltaPlusPrimary->SetPoint(i,0.1+i*0.1,plusDeltaPrimaryMean[i]);
    graphDeltaPlusPrimary->SetPointError(i,0,plusDeltaPrimaryError[i]);
    graphSigmaPlusPrimary->SetPoint(i,0.1+i*0.1,plusSigmaPrimaryMean[i]);
    graphSigmaPlusPrimary->SetPointError(i,0,plusSigmaPrimaryError[i]);

    graphMeanPlusSim->SetPoint(i,0.1+i*0.1,plusMultSimMean[i]);
    graphMeanPlusSim->SetPointError(i,0,plusMultSimError[i]);
    graphOmegaPlusSim->SetPoint(i,0.1+i*0.1,plusOmegaSimMean[i]);
    graphOmegaPlusSim->SetPointError(i,0,plusOmegaSimError[i]);
    graphDeltaPlusSim->SetPoint(i,0.1+i*0.1,plusDeltaSimMean[i]);
    graphDeltaPlusSim->SetPointError(i,0,plusDeltaSimError[i]);
    graphSigmaPlusSim->SetPoint(i,0.1+i*0.1,plusSigmaSimMean[i]);
    graphSigmaPlusSim->SetPointError(i,0,plusSigmaSimError[i]);


    graphMeanMinus->SetPoint(i,0.1+i*0.1,minusMultMean[i]);
    graphMeanMinus->SetPointError(i,0,minusMultError[i]);
    graphOmegaMinus->SetPoint(i,0.1+i*0.1,minusOmegaMean[i]);
    graphOmegaMinus->SetPointError(i,0,minusOmegaError[i]);
    graphDeltaMinus->SetPoint(i,0.1+i*0.1,minusDeltaMean[i]);
    graphDeltaMinus->SetPointError(i,0,minusDeltaError[i]);
    graphSigmaMinus->SetPoint(i,0.1+i*0.1,minusSigmaMean[i]);
    graphSigmaMinus->SetPointError(i,0,minusSigmaError[i]);

    graphMeanMinusPrimary->SetPoint(i,0.1+i*0.1,minusMultPrimaryMean[i]);
    graphMeanMinusPrimary->SetPointError(i,0,minusMultPrimaryError[i]);
    graphOmegaMinusPrimary->SetPoint(i,0.1+i*0.1,minusOmegaPrimaryMean[i]);
    graphOmegaMinusPrimary->SetPointError(i,0,minusOmegaPrimaryError[i]);
    graphDeltaMinusPrimary->SetPoint(i,0.1+i*0.1,minusDeltaPrimaryMean[i]);
    graphDeltaMinusPrimary->SetPointError(i,0,minusDeltaPrimaryError[i]);
    graphSigmaMinusPrimary->SetPoint(i,0.1+i*0.1,minusSigmaPrimaryMean[i]);
    graphSigmaMinusPrimary->SetPointError(i,0,minusSigmaPrimaryError[i]);

    graphMeanMinusSim->SetPoint(i,0.1+i*0.1,minusMultSimMean[i]);
    graphMeanMinusSim->SetPointError(i,0,minusMultSimError[i]);
    graphOmegaMinusSim->SetPoint(i,0.1+i*0.1,minusOmegaSimMean[i]);
    graphOmegaMinusSim->SetPointError(i,0,minusOmegaSimError[i]);
    graphDeltaMinusSim->SetPoint(i,0.1+i*0.1,minusDeltaSimMean[i]);
    graphDeltaMinusSim->SetPointError(i,0,minusDeltaSimError[i]);
    graphSigmaMinusSim->SetPoint(i,0.1+i*0.1,minusSigmaSimMean[i]);
    graphSigmaMinusSim->SetPointError(i,0,minusSigmaSimError[i]);

  }

  TLegendEntry *ler = leg->AddEntry(graphMeanTotal,legEntryRec,"lp");
  ler->SetTextColor(kBlack);
  TLegendEntry *les = leg->AddEntry(graphMeanTotalSim,legEntrySim,"lp");
  les->SetTextColor(kBlue);

  graphMeanTotal->GetXaxis()->SetTitle("#delta #eta");
  graphOmegaTotal->GetXaxis()->SetTitle("#delta #eta");
  graphDeltaTotal->GetXaxis()->SetTitle("#delta #eta");
  graphSigmaTotal->GetXaxis()->SetTitle("#delta #eta");

  graphMeanTotal->GetYaxis()->SetTitle("#LT N#GT");
  graphOmegaTotal->GetYaxis()->SetTitle("#omega[N]");
  graphDeltaTotal->GetYaxis()->SetTitle("#Delta[N,P_{T}]");
  graphSigmaTotal->GetYaxis()->SetTitle("#Sigma[N,P_{T}]");

  graphMeanTotal->SetTitle("#LT N#GT");
  graphOmegaTotal->SetTitle("#omega[N]");
  graphDeltaTotal->SetTitle("#Delta[N,P_{T}]");
  graphSigmaTotal->SetTitle("#Sigma[N,P_{T}]");

  TCanvas* cSigmaTotal = new TCanvas("cSigmaTotal","cSigmaTotal",1000,1000);
  cSigmaTotal->Divide(2,2);
  cSigmaTotal->cd(1);
  graphMeanTotal->Draw();
  graphMeanTotalSim->SetLineColor(kBlue);
  graphMeanTotalSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaTotal->cd(2);
  graphOmegaTotal->Draw();
  graphOmegaTotalSim->SetLineColor(kBlue);
  graphOmegaTotalSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaTotal->cd(3);

  graphDeltaTotal->SetMinimum(0.7);
  graphDeltaTotal->SetMaximum(1.5);

  graphDeltaTotal->Draw();
  graphDeltaTotalSim->SetLineColor(kBlue);
  graphDeltaTotalSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaTotal->cd(4);

  graphSigmaTotal->SetMinimum(0.96);
  graphSigmaTotal->SetMaximum(1.11);

  graphSigmaTotal->Draw();
  graphSigmaTotalSim->SetLineColor(kBlue);
  graphSigmaTotalSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaTotal->SaveAs(rootshortbeginning+"cSigmaTotal"+".pdf");



  graphMeanPlus->GetXaxis()->SetTitle("#delta #eta");
  graphOmegaPlus->GetXaxis()->SetTitle("#delta #eta");
  graphDeltaPlus->GetXaxis()->SetTitle("#delta #eta");
  graphSigmaPlus->GetXaxis()->SetTitle("#delta #eta");

  graphMeanPlus->GetYaxis()->SetTitle("#LT N^{+}#GT");
  graphOmegaPlus->GetYaxis()->SetTitle("#omega[N^{+}]");
  graphDeltaPlus->GetYaxis()->SetTitle("#Delta[N^{+},P_{T}^{+}]");
  graphSigmaPlus->GetYaxis()->SetTitle("#Sigma[N^{+},P_{T}^{+}]");

  graphMeanPlus->SetTitle("#LT N^{+}#GT");
  graphOmegaPlus->SetTitle("#omega[N^{+}]");
  graphDeltaPlus->SetTitle("#Delta[N^{+},P_{T}^{+}]");
  graphSigmaPlus->SetTitle("#Sigma[N^{+},P_{T}^{+}]");






  TCanvas* cSigmaPlus = new TCanvas("cSigmaPlus","cSigmaPlus",1000,1000);
  cSigmaPlus->Divide(2,2);
  cSigmaPlus->cd(1);
  graphMeanPlus->Draw();
  graphMeanPlusSim->SetLineColor(kBlue);
  graphMeanPlusSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaPlus->cd(2);
  graphOmegaPlus->Draw();
  graphOmegaPlusSim->SetLineColor(kBlue);
  graphOmegaPlusSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaPlus->cd(3);
  graphDeltaPlus->SetMinimum(0.7);
  graphDeltaPlus->SetMaximum(1.5);

  graphDeltaPlus->Draw();
  graphDeltaPlusSim->SetLineColor(kBlue);
  graphDeltaPlusSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaPlus->cd(4);

  graphSigmaPlus->SetMinimum(0.96);
  graphSigmaPlus->SetMaximum(1.11);

  graphSigmaPlus->Draw();
  graphSigmaPlusSim->SetLineColor(kBlue);
  graphSigmaPlusSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaPlus->SaveAs(rootshortbeginning+"cSigmaPlus.pdf");




  graphMeanMinus->GetXaxis()->SetTitle("#delta #eta");
  graphOmegaMinus->GetXaxis()->SetTitle("#delta #eta");
  graphDeltaMinus->GetXaxis()->SetTitle("#delta #eta");
  graphSigmaMinus->GetXaxis()->SetTitle("#delta #eta");

  graphMeanMinus->GetYaxis()->SetTitle("#LT N^{-}#GT");
  graphOmegaMinus->GetYaxis()->SetTitle("#omega[N^{-}]");
  graphDeltaMinus->GetYaxis()->SetTitle("#Delta[N^{-},P_{T}^{-}]");
  graphSigmaMinus->GetYaxis()->SetTitle("#Sigma[N^{-},P_{T}^{-}]");

  graphMeanMinus->SetTitle("#LT N^{-}#GT");
  graphOmegaMinus->SetTitle("#omega[N^{-}]");
  graphDeltaMinus->SetTitle("#Delta[N^{-},P_{T}^{-}]");
  graphSigmaMinus->SetTitle("#Sigma[N^{-},P_{T}^{-}]");

  TCanvas* cSigmaMinus = new TCanvas("cSigmaMinus","cSigmaMinus",1000,1000);
  cSigmaMinus->Divide(2,2);
  cSigmaMinus->cd(1);
  graphMeanMinus->Draw();
  graphMeanMinusSim->SetLineColor(kBlue);
  graphMeanMinusSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaMinus->cd(2);
  graphOmegaMinus->Draw();
  graphOmegaMinusSim->SetLineColor(kBlue);
  graphOmegaMinusSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaMinus->cd(3);
  graphDeltaMinus->SetMinimum(0.7);
  graphDeltaMinus->SetMaximum(1.5);

  graphDeltaMinus->Draw();
  graphDeltaMinusSim->SetLineColor(kBlue);
  graphDeltaMinusSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaMinus->cd(4);

  graphSigmaMinus->SetMinimum(0.96);
  graphSigmaMinus->SetMaximum(1.11);

  graphSigmaMinus->Draw();
  graphSigmaMinusSim->SetLineColor(kBlue);
  graphSigmaMinusSim->Draw("Same");
  leg->Draw();
  label.DrawLatex(0.27,0.93, "0.2<p_{T}<2 GeV/c");
  label.DrawLatex(0.5,0.87, "Bi+Bi@#sqrt{s_{NN}}=9.46 GeV, 0-10% centrality");
  cSigmaMinus->SaveAs(rootshortbeginning+"cSigmaMinus.pdf");





}





void readTree(){


  TFile *f = new TFile("smash_mpdroot_bibi_9point46.root");
  TTree *t1 = (TTree*)f->Get("tree_sim_rec");

  Int_t nEvents = (Int_t)t1->GetEntries();
  std::cout << nEvents << " events in root file" << std::endl;

  Int_t nSubSamples;
  Int_t nEventsInSubSample;

  nSubSamples = 30;
  nEventsInSubSample = nEvents/nSubSamples;
  std::cout << "now analyzing with " << nSubSamples << " nSubSamples of " << nEventsInSubSample << " events" << std::endl;
  //analyzeTreeWithSubsamples(t1,nSubSamples,nEventsInSubSample);
  analyzeTree(t1,nSubSamples,nEventsInSubSample);

  return;



}
