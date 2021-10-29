// ==================================================================================================================================================================================
// This code performes 1D-unfolding of Multiplicity distribution using the pre-prepared tree that contains sim and rec info about the events to build detector Response Matrix (RM).
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
#include "TSystem.h"

//***************************************************************************************************************************
// to run locally do
// root -l
// gSystem->Load("/*add_your_path*/RooUnfold/libRooUnfold");
// .L final_unfolding_N.cc+
// final_unfolding_N( 1, 50 )
//***************************************************************************************************************************

// before doing unfolding the RooUnfold package should be installed from https://gitlab.cern.ch/RooUnfold/RooUnfold

#include "/*add_your_path*/RooUnfold/src/RooUnfoldResponse.h"
#include "/*add_your_path*/RooUnfold/src/RooUnfoldBayes.h"

using namespace std;
using namespace TMath;

// n_subsamples - defines the number of bootstrap samples. n_subsamples = 1 will be used to calculate the central value without bootstrap, other values of n_subsamples will be used to calculate statistical uncertainty. n_subsamples is used only for data. No bootstrap is applied to the RM creation: all bootstrap data histo are unfolded by the same RM.

//n_unf_iterations - defines the number of unfolding iterations in the iterative bayesian unfolding. Be aware that too big number of iterations might start to increase the error of unfolding from some point.

int final_unfolding_N(int n_subsamples, int n_unf_iterations)
{
    int NumbOne_winds = 10; // number of pseudorapidity intervals
    int charge_combo = 3; // number of charge combinations: all, pos, neg
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
    int multOne_rec[NumbOne_winds],multOnePlus_rec[NumbOne_winds], multOneMinus_rec[NumbOne_winds], multOne_sim[NumbOne_winds],multOnePlus_sim[NumbOne_winds], multOneMinus_sim[NumbOne_winds];
    
    TRandom3 rand;
    rand.SetSeed(0);
    
    // access tree with sim and rec data
    TTree *T = (TTree*)file_to_an->Get(name);
    T->SetBranchAddress("multOne_rec", multOne_rec);
    T->SetBranchAddress("multOnePlus_rec", multOnePlus_rec);
    T->SetBranchAddress("multOneMinus_rec", multOneMinus_rec);
    
    T->SetBranchAddress("multOne_sim", multOne_sim);
    T->SetBranchAddress("multOnePlus_sim", multOnePlus_sim);
    T->SetBranchAddress("multOneMinus_sim", multOneMinus_sim);
    
    // create histograms to calculate value for RM and value+error from subsamples for data
    TH1D*                   fHistMultONE_sim_MC[NumbOne_winds][charge_combo];
    TH1D*                   fHistMultONE_sim_data[NumbOne_winds][charge_combo][n_subsamples];
    TH1D*                   fHistMultONE_rec_MC[NumbOne_winds][charge_combo];
    TH1D*                   fHistMultONE_rec_data[NumbOne_winds][charge_combo][n_subsamples];
    TH2D *hRM = new TH2D("hRM","Response Matrix", 1000, -0.5, 999.5, 1000, -0.5, 999.5);
    
    TString ONEname_hist_sim_MC = "fHistMultONE_sim_MC";
    TString ONEname_hist_sim_data = "fHistMultONE_sim_data";
    TString ONEname_hist_rec_MC = "fHistMultONE_rec_MC";
    TString ONEname_hist_rec_data = "fHistMultONE_rec_data";
    
    char windONE_numb[50];
    TString ONEname_sim_MC[NumbOne_winds][charge_combo],ONEname_rec_MC[NumbOne_winds][charge_combo],ONEname_sim_data[NumbOne_winds][charge_combo][n_subsamples],ONEname_rec_data[NumbOne_winds][charge_combo][n_subsamples];
    
    TString name_charge[charge_combo];
    name_charge[0] = "_all";
    name_charge[1] = "_pos";
    name_charge[2] = "_neg";
    
    char sub_numb[50];
    char wind_numb[50];

    for(int i = 0; i<NumbOne_winds;i++)
    {
        sprintf(windONE_numb, "_%i", i);
        for(int k = 0;k<charge_combo;k++)
        {
            ONEname_sim_MC[i][k] = ONEname_hist_sim_MC + windONE_numb+name_charge[k];
            ONEname_rec_MC[i][k] = ONEname_hist_rec_MC + windONE_numb+name_charge[k];
            fHistMultONE_sim_MC[i][k] = new TH1D(ONEname_sim_MC[i][k], ONEname_sim_MC[i][k], 1000, -0.5, 999.5);
            fHistMultONE_rec_MC[i][k] = new TH1D(ONEname_rec_MC[i][k], ONEname_rec_MC[i][k],1000, -0.5, 999.5);
            
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                sprintf(sub_numb, "_%i", sub_sample);
                ONEname_sim_data[i][k][sub_sample] = ONEname_hist_sim_data + windONE_numb+name_charge[k]+sub_numb;
                ONEname_rec_data[i][k][sub_sample] = ONEname_hist_rec_data + windONE_numb+name_charge[k]+sub_numb;
                
                fHistMultONE_sim_data[i][k][sub_sample] = new TH1D(ONEname_sim_data[i][k][sub_sample], ONEname_sim_data[i][k][sub_sample], 1000, -0.5, 999.5);
                fHistMultONE_rec_data[i][k][sub_sample] = new TH1D(ONEname_rec_data[i][k][sub_sample], ONEname_rec_data[i][k][sub_sample], 1000, -0.5, 999.5);
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
    cout << "NumbOne_winds = "<<NumbOne_winds<<endl;
    cout<<" "<<endl;
    cout << "n_subsamples = "<<n_subsamples<<endl;
    cout << "n_unf_iterations = "<<n_unf_iterations<<endl;
    
    // create response matrix based on rec and sim histograms
    RooUnfoldResponse* responseOneWind[NumbOne_winds][charge_combo];
    for(int k = 0; k<NumbOne_winds;k++)
    {
        for(int charge = 0; charge<charge_combo;charge++)
        {
            responseOneWind[k][charge] = new RooUnfoldResponse(fHistMultONE_rec_MC[k][charge], fHistMultONE_sim_MC[k][charge]);
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
            for(int i = 0; i<NumbOne_winds;i++)
            {
                multOne_rec[i]=0;multOnePlus_rec[i]=0;multOneMinus_rec[i]=0;
                multOne_sim[i]=0;multOnePlus_sim[i]=0;multOneMinus_sim[i]=0;
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
            for(int i = 0; i<NumbOne_winds;i++)
            {
                fHistMultONE_rec_data[i][0][sub_sample]->Fill(multOne_rec[i]);
                fHistMultONE_sim_data[i][0][sub_sample]->Fill(multOne_sim[i]);

                fHistMultONE_rec_data[i][1][sub_sample]->Fill(multOnePlus_rec[i]);
                fHistMultONE_sim_data[i][1][sub_sample]->Fill(multOnePlus_sim[i]);

                fHistMultONE_rec_data[i][2][sub_sample]->Fill(multOneMinus_rec[i]);
                fHistMultONE_sim_data[i][2][sub_sample]->Fill(multOneMinus_sim[i]);
            }
        }
    }
    
    // filling RM created previously with rec and sim MC
    cout<<"start feeling RM"<<endl;
    for (int j = nEv_data; j < entries; j++)
    {
        for(int i = 0; i<NumbOne_winds;i++)
        {
            multOne_rec[i]=0;multOnePlus_rec[i]=0;multOneMinus_rec[i]=0;
            multOne_sim[i]=0;multOnePlus_sim[i]=0;multOneMinus_sim[i]=0;
        }

        T->GetEntry(j);

        for(int i = 0; i<NumbOne_winds;i++)
        {
            fHistMultONE_rec_MC[i][0]->Fill(multOne_rec[i]);
            fHistMultONE_sim_MC[i][0]->Fill(multOne_sim[i]);
            
            fHistMultONE_rec_MC[i][1]->Fill(multOnePlus_rec[i]);
            fHistMultONE_sim_MC[i][1]->Fill(multOnePlus_sim[i]);
            
            fHistMultONE_rec_MC[i][2]->Fill(multOneMinus_rec[i]);
            fHistMultONE_sim_MC[i][2]->Fill(multOneMinus_sim[i]);
            
            // if value of the variable is equal to -1 in REC that means the event was lost by rec event cuts -> we should put it in "miss" to normalize response matrix
            if(multOne_rec[i]==-1)
            {
                responseOneWind[i][0]->Miss(multOne_sim[i]);
            }
            else
            {
                responseOneWind[i][0]->Fill(multOne_rec[i], multOne_sim[i]);
                if(i==0)
                {hRM->Fill(multOne_rec[0],multOne_sim[0]);}
            }
            
            if(multOnePlus_rec[i]==-1)
            {
                responseOneWind[i][1]->Miss(multOnePlus_sim[i]);
            }
            else
            {
                responseOneWind[i][1]->Fill(multOnePlus_rec[i], multOnePlus_sim[i]);
            }
            
            if(multOneMinus_rec[i]==-1)
            {
                responseOneWind[i][2]->Miss(multOneMinus_sim[i]);
            }
            else
            {
                responseOneWind[i][2]->Fill(multOneMinus_rec[i], multOneMinus_sim[i]);
            }
        }
    } //end of data events
    
    // create unfold function (here you can change RooUnfold method: bin-by-bin, RooUnfoldBayes (IBU), RooUnfoldSvd) using RM, biased rec data to unfold and the number of iterations in the case of IBU
    RooUnfoldBayes*   unfoldOneWind[NumbOne_winds][charge_combo][n_subsamples];
    for(int m = 0; m<NumbOne_winds;m++)
    {
        for(int charge = 0; charge<charge_combo;charge++)
        {
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                unfoldOneWind[m][charge][sub_sample] = new RooUnfoldBayes(responseOneWind[m][charge], fHistMultONE_rec_data[m][charge][sub_sample], n_unf_iterations);
            }
        }
    }
    
    // create histogram for unfolded data, perform unfolding and rename unfolded histo
    TH1D* hRecoOneWind[NumbOne_winds][charge_combo][n_subsamples];
    TString hist_ONE_name_unf = "unfolded_ONE";
    TString ONE_name_unf[NumbOne_winds][charge_combo][n_subsamples];
    
    cout<<"ONE wind unfolding"<<endl;
    for(int j = 0; j<NumbOne_winds;j++)
    {
        sprintf(wind_numb, "_%i", j);
        for(int charge = 0; charge<charge_combo;charge++)
        {
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                sprintf(sub_numb, "_%i", sub_sample);
                ONE_name_unf[j][charge][sub_sample] = hist_ONE_name_unf+wind_numb+name_charge[charge]+sub_numb;
                hRecoOneWind[j][charge][sub_sample] = (TH1D*) unfoldOneWind[j][charge][sub_sample]->Hreco();
                hRecoOneWind[j][charge][sub_sample]-> SetName(ONE_name_unf[j][charge][sub_sample]);
            }
        }
    }
    
    // create output file and fill it with pure, biased and unfolded data
    TString filenameONE;
    TString end = ".root";
    TString outpath = "unf_N_"; // put here the path where to save output root file
    TString Sboot;
    TString Sunf;
    
    Sboot = to_string(n_subsamples);
    Sunf = to_string(n_unf_iterations);
    filenameONE = outpath+"ONEwindN_"+"_"+Sboot+"_"+Sunf+end;
    
    TFile *fFileONE = new TFile(filenameONE, "RECREATE");
    for(int j = 0; j<NumbOne_winds;j++)
    {
        for(int charge = 0; charge<charge_combo;charge++)
        {
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                fHistMultONE_sim_data[j][charge][sub_sample]->Write();
                fHistMultONE_rec_data[j][charge][sub_sample]->Write();
                hRecoOneWind[j][charge][sub_sample]->Write();
            }
        }
    }
    fFileONE->Close();
    
    cout<<"job is done!"<<endl;
    
    return 0;
}
