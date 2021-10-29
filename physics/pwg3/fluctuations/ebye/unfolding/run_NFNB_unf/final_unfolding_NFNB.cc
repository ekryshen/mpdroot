// ==================================================================================================================================================================================
// This code performes 2D-unfolding of NF-NB distribution (the correlation between multiplicities calculated in separated forward and backward pseudorapidity intervals) using the pre-prepared tree that contains sim and rec info about the events to build detector Response Matrix (RM).
// In this code the pseudo pure and rec data is sampled from the same tree (15% for pseudo data, the rest to build RM). In general, for this purpose another Event generator
// can be used. Or this code can be in future applied to the data itself to perform the corrections of the detector bias.

// Author: Daria Prokhorova (SPbSU), 10.2020
// ==================================================================================================================================================================================

#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3.h>
#include <TCutG.h>
#include <TProfile.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>
//#include <utl/Vector.h>
#include <TCanvas.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TKey.h>
#include <TROOT.h>
#include <TChain.h>
#include "TTree.h"
#include "TBranch.h"
#include <list>
#include "TString.h"
#include "TRandom3.h"

//***************************************************************************************************************************
// to run locally do
// root -l
// gSystem->Load("/*add_your_path*/Unfolding/RooUnfold/libRooUnfold");
// .L final_unfolding_NFNB.cc+
// .x final_unfolding_NFNB( 1, 20 )
//***************************************************************************************************************************

// before doing unfolding the RooUnfold package should be installed from https://gitlab.cern.ch/RooUnfold/RooUnfold

#include "/lhep/users/dprokhor/fairSoft/Unfolding/RooUnfold/src/RooUnfoldResponse.h"
#include "/lhep/users/dprokhor/fairSoft/Unfolding/RooUnfold/src/RooUnfoldBayes.h"

using namespace std;
using namespace TMath;

// n_subsamples - defines the number of bootstrap samples. n_subsamples = 1 will be used to calculate the central value without bootstrap, other values of n_subsamples will be used to calculate statistical uncertainty. n_subsamples is used only for "data". No bootstrap is applied to the RM creation: all bootstrap "data" histo are unfolded by the same RM - as in the first approximation

// n_unf_iterations - defines the number of unfolding iterations in the iterative bayesian unfolding. Be aware that too big number of iterations might start to increase the error of unfolding from some point (when unfolding starts to recognize some structures in statistical fluctuations)

int final_unfolding_NFNB(int n_subsamples, int n_unf_iterations)
{
    int Numb_winds = 10; // number of pseudorapidity intervals pairs
    int charge_combo_FB = 4; // number of charge combinations: FallBall, FposBpos, FnegBneg, FposBneg
    double fraction = 0.15; // set fraction of all MC to be pseudo-data (should be less then 0.3)
    
    // the file to extract pseudo-data and MC for RM
     TFile *file_to_an = new TFile("out.root"); // put here the path to the root file with tree of rec/sim info about each event
    TKey *key = 0;
    const char *name = 0;
    TIter next(file_to_an->GetListOfKeys());
    key = (TKey*)next();
    name = key->GetName();
    cout << "name " << name << endl;
    
    // create variables to study
    int multF_rec[Numb_winds],multB_rec[Numb_winds], multPlusF_rec[Numb_winds],multPlusB_rec[Numb_winds], multMinusF_rec[Numb_winds], multMinusB_rec[Numb_winds],multF_sim[Numb_winds],multB_sim[Numb_winds], multPlusF_sim[Numb_winds],multPlusB_sim[Numb_winds], multMinusF_sim[Numb_winds], multMinusB_sim[Numb_winds];
    
    TRandom3 rand;
    rand.SetSeed(0);
    
    TString filenameFB;
    TString end = ".root";
    TString outpath = "unf_NFNB_"; // put here the path where to save output root file
    TString Sboot;
    TString Sunf;
    Sboot = to_string(n_subsamples);
    Sunf = to_string(n_unf_iterations);
    filenameFB = outpath+Sboot+"_"+Sunf+end;
    
    cout<<"file writing"<<endl;
    TFile *fFileFB = new TFile(filenameFB, "RECREATE");
    
    // access tree with sim and rec data
    TTree *T = (TTree*)file_to_an->Get(name);
    T->SetBranchAddress("multF_rec", multF_rec);
    T->SetBranchAddress("multB_rec", multB_rec);
    T->SetBranchAddress("multPlusF_rec", multPlusF_rec);
    T->SetBranchAddress("multMinusF_rec", multMinusF_rec);
    T->SetBranchAddress("multPlusB_rec", multPlusB_rec);
    T->SetBranchAddress("multMinusB_rec", multMinusB_rec);
    
    T->SetBranchAddress("multF_sim", multF_sim);
    T->SetBranchAddress("multB_sim", multB_sim);
    T->SetBranchAddress("multPlusF_sim", multPlusF_sim);
    T->SetBranchAddress("multMinusF_sim", multMinusF_sim);
    T->SetBranchAddress("multPlusB_sim", multPlusB_sim);
    T->SetBranchAddress("multMinusB_sim", multMinusB_sim);
    
    // create histograms to calculate value for RM and value+error from subsamples for data
    TH2D*                   fHistMultFB_sim_MC[Numb_winds][charge_combo_FB];
    TH2D*                   fHistMultFB_sim_data[Numb_winds][charge_combo_FB][n_subsamples];
    
    TH2D*                   fHistMultFB_rec_MC[Numb_winds][charge_combo_FB];
    TH2D*                   fHistMultFB_rec_data[Numb_winds][charge_combo_FB][n_subsamples];
    
    TString FBname_hist_sim_MC = "fHistMultFB_sim_MC";
    TString FBname_hist_sim_data = "fHistMultFB_sim_data";
    TString FBname_hist_rec_MC = "fHistMultFB_rec_MC";
    TString FBname_hist_rec_data = "fHistMultFB_rec_data";
    
    char wind_numb[50];
    TString FBname_sim_MC[Numb_winds][charge_combo_FB],FBname_rec_MC[Numb_winds][charge_combo_FB],FBname_sim_data[Numb_winds][charge_combo_FB][n_subsamples],FBname_rec_data[Numb_winds][charge_combo_FB][n_subsamples];
    
    TString name_charge[charge_combo_FB];
    name_charge[0] = "_all";
    name_charge[1] = "_pos";
    name_charge[2] = "_neg";
    name_charge[3] = "_FposBneg";
    char sub_numb[50];
    
    int n_bins_all = 100;
    int n_bins_charges = 55;
    
    double bin_all_max = 99.5;
    double bin_charge_max = 54.5;
    
    double bin_all_min = -0.5;
    double bin_charge_min =  -0.5;
    
    for(int i = 0; i<Numb_winds;i++)
    {
        sprintf(wind_numb, "_%i", i);
        for(int k = 0;k<charge_combo_FB;k++)
        {
            FBname_sim_MC[i][k] = FBname_hist_sim_MC + wind_numb+name_charge[k];
            FBname_rec_MC[i][k] = FBname_hist_rec_MC + wind_numb+name_charge[k];
            
            if(k==0)
            {
            fHistMultFB_sim_MC[i][k] = new TH2D(FBname_sim_MC[i][k], FBname_sim_MC[i][k], n_bins_all, bin_all_min, bin_all_max, n_bins_all, bin_all_min, bin_all_max);
            fHistMultFB_rec_MC[i][k] = new TH2D(FBname_rec_MC[i][k], FBname_rec_MC[i][k],n_bins_all, bin_all_min, bin_all_max, n_bins_all, bin_all_min, bin_all_max);
            }
            else
            {
                fHistMultFB_sim_MC[i][k] = new TH2D(FBname_sim_MC[i][k], FBname_sim_MC[i][k], n_bins_charges, bin_charge_min, bin_charge_max, n_bins_charges, bin_charge_min, bin_charge_max);
                fHistMultFB_rec_MC[i][k] = new TH2D(FBname_rec_MC[i][k], FBname_rec_MC[i][k],n_bins_charges, bin_charge_min, bin_charge_max, n_bins_charges, bin_charge_min, bin_charge_max);
            }
            
            
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                sprintf(sub_numb, "_%i", sub_sample);
                FBname_sim_data[i][k][sub_sample] = FBname_hist_sim_data + wind_numb+name_charge[k]+sub_numb;
                FBname_rec_data[i][k][sub_sample] = FBname_hist_rec_data + wind_numb+name_charge[k]+sub_numb;
            
                if(k==0)
                {
                fHistMultFB_sim_data[i][k][sub_sample] = new TH2D(FBname_sim_data[i][k][sub_sample], FBname_sim_data[i][k][sub_sample], n_bins_all, bin_all_min, bin_all_max, n_bins_all, bin_all_min, bin_all_max);
                fHistMultFB_rec_data[i][k][sub_sample] = new TH2D(FBname_rec_data[i][k][sub_sample], FBname_rec_data[i][k][sub_sample], n_bins_all, bin_all_min, bin_all_max, n_bins_all, bin_all_min, bin_all_max);
                }
                else
                {
                    fHistMultFB_sim_data[i][k][sub_sample] = new TH2D(FBname_sim_data[i][k][sub_sample], FBname_sim_data[i][k][sub_sample], n_bins_charges, bin_charge_min, bin_charge_max, n_bins_charges, bin_charge_min, bin_charge_max);
                    fHistMultFB_rec_data[i][k][sub_sample] = new TH2D(FBname_rec_data[i][k][sub_sample], FBname_rec_data[i][k][sub_sample], n_bins_charges, bin_charge_min, bin_charge_max, n_bins_charges, bin_charge_min, bin_charge_max);
                }
            }
        }
    }
    
    // setup data and MC events
    int entries;
    entries = T->GetEntries();
    int nEv_data = entries*fraction;
    int nEv_MC = entries-nEv_data;
    
    cout << "entries all = " << entries << endl;
    cout << "entries data = " << nEv_data << endl;
    cout << "entries MC = " << nEv_MC << endl;
    cout<<" "<<endl;
    cout << "Numb_winds = "<<Numb_winds<<endl;
    cout<<" "<<endl;
    cout << "n_subsamples = "<<n_subsamples<<endl;
    cout << "n_unf_iterations = "<<n_unf_iterations<<endl;
    
    // create response matrix based on rec and sim histograms
    RooUnfoldResponse* responseForBack[Numb_winds][charge_combo_FB];
    for(int k = 0; k<Numb_winds;k++)
    {
        for(int charge = 0; charge<charge_combo_FB;charge++)
        {
            responseForBack[k][charge] = new RooUnfoldResponse(fHistMultFB_rec_MC[k][charge], fHistMultFB_sim_MC[k][charge]);
        }
    }
    
    // start filling data histograms (bootstrap method to calculate errors, n_subsamples = 1 corresponds to the full statistics to obtain the value, other n_subsamples to calculate errors)
    int n_ev_in_subsample;
    cout<<"start data"<<endl;
    
    for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
    {
        cout<<"subsample #"<<sub_sample<<endl;
        n_ev_in_subsample=0;
        
        while(n_ev_in_subsample<nEv_data)
        {
            for(int i = 0; i<Numb_winds;i++)
            {
                multB_rec[i]=0;multPlusB_rec[i]=0;multMinusB_rec[i]=0;
                multF_rec[i]=0;multPlusF_rec[i]=0;multMinusF_rec[i]=0;
                multB_sim[i]=0;multPlusB_sim[i]=0;multMinusB_sim[i]=0;
                multF_sim[i]=0;multPlusF_sim[i]=0;multMinusF_sim[i]=0;
            }
            
            if(n_subsamples==1)
            {
                T->GetEntry(n_ev_in_subsample);
            }
            else
            {
                T->GetEntry(rand.Integer(nEv_data));
            }
            
            n_ev_in_subsample++;
            // filling pure (sim) and biased (rec) data histograms
            for(int i = 0; i<Numb_winds;i++)
            {
                fHistMultFB_rec_data[i][0][sub_sample]->Fill(multF_rec[i], multB_rec[i]);
                fHistMultFB_sim_data[i][0][sub_sample]->Fill(multF_sim[i], multB_sim[i]);
                
                fHistMultFB_rec_data[i][1][sub_sample]->Fill(multPlusF_rec[i], multPlusB_rec[i]);
                fHistMultFB_sim_data[i][1][sub_sample]->Fill(multPlusF_sim[i], multPlusB_sim[i]);
                
                fHistMultFB_rec_data[i][2][sub_sample]->Fill(multMinusF_rec[i], multMinusB_rec[i]);
                fHistMultFB_sim_data[i][2][sub_sample]->Fill(multMinusF_sim[i], multMinusB_sim[i]);
                
                fHistMultFB_rec_data[i][3][sub_sample]->Fill(multPlusF_rec[i], multMinusB_rec[i]);
                fHistMultFB_sim_data[i][3][sub_sample]->Fill(multPlusF_sim[i], multMinusB_sim[i]);
            }
        }
    }
    
    // filling RM created previously with rec and sim MC
    cout<<"start feeling RM"<<endl;
    for (int j = nEv_data; j < entries; j++)
    {
        for(int i = 0; i<Numb_winds;i++)
        {
            multB_rec[i]=0;multPlusB_rec[i]=0;multMinusB_rec[i]=0;
            multF_rec[i]=0;multPlusF_rec[i]=0;multMinusF_rec[i]=0;
            multB_sim[i]=0;multPlusB_sim[i]=0;multMinusB_sim[i]=0;
            multF_sim[i]=0;multPlusF_sim[i]=0;multMinusF_sim[i]=0;
        }
        
        T->GetEntry(j);
        
        for(int i = 0; i<Numb_winds;i++)
        {
            fHistMultFB_rec_MC[i][0]->Fill(multF_rec[i], multB_rec[i]);
            fHistMultFB_sim_MC[i][0]->Fill(multF_sim[i], multB_sim[i]);
            
            fHistMultFB_rec_MC[i][1]->Fill(multPlusF_rec[i], multPlusB_rec[i]);
            fHistMultFB_sim_MC[i][1]->Fill(multPlusF_sim[i], multPlusB_sim[i]);
            
            fHistMultFB_rec_MC[i][2]->Fill(multMinusF_rec[i], multMinusB_rec[i]);
            fHistMultFB_sim_MC[i][2]->Fill(multMinusF_sim[i], multMinusB_sim[i]);
            
            fHistMultFB_rec_MC[i][3]->Fill(multPlusF_rec[i], multMinusB_rec[i]);
            fHistMultFB_sim_MC[i][3]->Fill(multPlusF_sim[i], multMinusB_sim[i]);
            
             // if value of the variable is equal to -1 in REC that means the event was lost by rec event cuts -> we should put it in "miss" to normalize histogram
            if((multF_rec[i]==-1)&&(multB_rec[i]==-1))
               {
                   responseForBack[i][0]->Miss(multF_sim[i], multB_sim[i]);
               }
               else
               {
                   responseForBack[i][0]->Fill(multF_rec[i], multB_rec[i], multF_sim[i], multB_sim[i]);
               }
               
               if((multPlusF_rec[i]==-1)&&(multPlusB_rec[i]==-1))
               {
                   responseForBack[i][1]->Miss(multPlusF_sim[i], multPlusB_sim[i]);
               }
               else
               {
                   responseForBack[i][1]->Fill(multPlusF_rec[i], multPlusB_rec[i], multPlusF_sim[i], multPlusB_sim[i]);
               }
               
               if((multMinusF_rec[i]==-1)&&(multMinusB_rec[i]==-1))
               {
                   responseForBack[i][2]->Miss(multMinusF_sim[i], multMinusB_sim[i]);
               }
               else
               {
                   responseForBack[i][2]->Fill(multMinusF_rec[i], multMinusB_rec[i], multMinusF_sim[i], multMinusB_sim[i]);
               }
               
               if((multPlusF_rec[i]==-1)&&(multMinusB_rec[i]==-1))
               {
                   responseForBack[i][3]->Miss(multPlusF_sim[i], multMinusB_sim[i]);
               }
               else
               {
                   responseForBack[i][3]->Fill(multPlusF_rec[i], multMinusB_rec[i], multPlusF_sim[i], multMinusB_sim[i]);
               }
               
        }// end of windows
   }//end of data events
               
               // create unfold function (here you can change RooUnfold method: bin-by-bin, RooUnfoldBayes, RooUnfoldSvd) using RM, biased rec data to unfold and the number of iterations in the case of IBU
               RooUnfoldBayes*   unfoldForBack[Numb_winds][charge_combo_FB][n_subsamples];
               
               for(int m = 0; m<Numb_winds;m++)
               {
                   for(int charge = 0; charge<charge_combo_FB;charge++)
                   {
                       for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
                       {
                           unfoldForBack[m][charge][sub_sample] = new RooUnfoldBayes(responseForBack[m][charge], fHistMultFB_rec_data[m][charge][sub_sample], n_unf_iterations);
                       }
                   }
               }
               
               cout<<"unfold objects created"<<endl;
               
               // create histogram for unfolded data, perform unfolding and rename unfolded histo
               TH2D* hRecoForBack[Numb_winds][charge_combo_FB][n_subsamples];
               
               TString hist_FB_name_unf = "unfolded_FB";
               TString FB_name_unf[Numb_winds][charge_combo_FB][n_subsamples];
               cout<<"FB unfolding"<<endl;
               
               for(int j = 0; j<Numb_winds;j++)
               {
                   sprintf(wind_numb, "_%i", j);
                   for(int charge = 0; charge<charge_combo_FB;charge++)
                   {
                       for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
                       {
                           hRecoForBack[j][charge][sub_sample]=(TH2D*) unfoldForBack[j][charge][sub_sample]->Hreco();
                           sprintf(sub_numb, "_%i", sub_sample);
                           FB_name_unf[j][charge][sub_sample] = hist_FB_name_unf+wind_numb+name_charge[charge]+sub_numb;
                           hRecoForBack[j][charge][sub_sample]-> SetName(FB_name_unf[j][charge][sub_sample]);
                       }
                   }
               }
               
               // create output file and fill it with pure, biased and unfolded data
    
               for(int j = 0; j<Numb_winds;j++)
               {
                   for(int charge = 0; charge<charge_combo_FB;charge++)
                   {
                       for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
                       {
                           fHistMultFB_sim_data[j][charge][sub_sample]->Write();
                           fHistMultFB_rec_data[j][charge][sub_sample]->Write();
                           hRecoForBack[j][charge][sub_sample]->Write();
                           
                        delete fHistMultFB_sim_data[j][charge][sub_sample];
                        delete fHistMultFB_rec_data[j][charge][sub_sample];
                        delete hRecoForBack[j][charge][sub_sample];
                           
                       }
                   }
               }
               fFileFB->Close();
               cout<<"job is done!"<<endl;
               
               return 0;
}
               
               
