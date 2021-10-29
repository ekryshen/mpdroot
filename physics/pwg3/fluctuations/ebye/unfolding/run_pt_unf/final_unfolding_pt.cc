// ==================================================================================================================================================================================
// This code performes 1D-unfolding of pT spectrum using the pre-prepared tree that contains sim and rec info about the tracks to build detector Response Matrix (RM).
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
#include <TSystem.h>
#include <TGraph.h>

//***************************************************************************************************************************
// to run do
// root -l
// gSystem->Load("/*put_your_path*/RooUnfold/libRooUnfold");
// .L final_unfolding_pt.cc+
// final_unfolding_pt( 100, 50 )
//***************************************************************************************************************************

// before doing unfolding the RooUnfold package should be installed from https://gitlab.cern.ch/RooUnfold/RooUnfold

#include "/*put_your_path*/RooUnfold/src/RooUnfoldResponse.h"
#include "/*put_your_path*/RooUnfold/src/RooUnfoldBayes.h"

using namespace std;
using namespace TMath;

// n_subsamples - defines the number of bootstrap samples. n_subsamples = 1 will be used to calculate the central value without bootstrap, other values of n_subsamples will be used to calculate statistical uncertainty. n_subsamples is used only for data. No bootstrap is applied to the RM creation: all bootstrap data histo are unfolded by the same RM.

//n_unf_iterations - defines the number of unfolding iterations in the iterative bayesian unfolding. Be aware that too big number of iterations might start to increase the error of unfolding from some point.

int final_unfolding_pt(int n_subsamples, int n_unf_iterations_1d)
{
    int NumbOne_winds = 10; // number of pseudorapidity intervals
    int charge_combo_pt = 3; // number of charge combinations: all, pos, neg
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
    vector<vector<Double_t> > * particles_pt_in_event_sim =0;
    vector<vector<Double_t> > * particles_ptPlus_in_event_sim =0;
    vector<vector<Double_t> > * particles_ptMinus_in_event_sim =0;
    
    vector<vector<Double_t> > * particles_pt_in_event_rec =0;
    vector<vector<Double_t> > * particles_ptPlus_in_event_rec =0;
    vector<vector<Double_t> > * particles_ptMinus_in_event_rec =0;
    
    TRandom3 rand;
    rand.SetSeed(0);
    
    // access tree with sim and rec data
    TTree *T = (TTree*)file_to_an->Get(name);
    T->SetBranchAddress("ptOne_sim", &particles_pt_in_event_sim);
    T->SetBranchAddress("ptOnePlus_sim", &particles_ptPlus_in_event_sim);
    T->SetBranchAddress("ptOneMinus_sim", &particles_ptMinus_in_event_sim);
    
    T->SetBranchAddress("ptOne_rec", &particles_pt_in_event_rec);
    T->SetBranchAddress("ptOnePlus_rec", &particles_ptPlus_in_event_rec);
    T->SetBranchAddress("ptOneMinus_rec", &particles_ptMinus_in_event_rec);
    
    // create histograms to calculate value for RM and value+error from subsamples for data
    TH1D*                   fHistpt_sim_MC[NumbOne_winds][charge_combo_pt];
    TH1D*                   fHistpt_sim_data[NumbOne_winds][charge_combo_pt][n_subsamples];
    
    TH1D*                   fHistpt_rec_MC[NumbOne_winds][charge_combo_pt];
    TH1D*                   fHistpt_rec_data[NumbOne_winds][charge_combo_pt][n_subsamples];
    TH2D *hRM = new TH2D("hRM","Response Matrix", 200,0, 2, 200,0, 2);

    TString ptname_hist_sim_MC = "fHistpt_sim_MC";
    TString ptname_hist_sim_data = "fHistpt_sim_data";
    TString ptname_hist_rec_MC = "fHistpt_rec_MC";
    TString ptname_hist_rec_data = "fHistpt_rec_data";
    
    char wind_numb[50];
    TString ptname_sim_MC[NumbOne_winds][charge_combo_pt], ptname_rec_MC[NumbOne_winds][charge_combo_pt], ptname_sim_data[NumbOne_winds][charge_combo_pt][n_subsamples], ptname_rec_data[NumbOne_winds][charge_combo_pt][n_subsamples];
    
    TString name_charge[charge_combo_pt];
    name_charge[0] = "_all";
    name_charge[1] = "_pos";
    name_charge[2] = "_neg";
   
    char sub_numb[50];
    
    for(int i = 0; i<NumbOne_winds;i++)
    {
        sprintf(wind_numb, "_%i", i);
        
        for(int k = 0;k<charge_combo_pt;k++)
        {
            ptname_sim_MC[i][k] =ptname_hist_sim_MC + wind_numb+name_charge[k];
            ptname_rec_MC[i][k] =ptname_hist_rec_MC + wind_numb+name_charge[k];
            
            fHistpt_sim_MC[i][k] = new TH1D(ptname_sim_MC[i][k], ptname_sim_MC[i][k], 200,0, 2);
            fHistpt_rec_MC[i][k] = new TH1D(ptname_rec_MC[i][k], ptname_rec_MC[i][k],200,0, 2);
            
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                sprintf(sub_numb, "_%i", sub_sample);
                
                ptname_sim_data[i][k][sub_sample] =ptname_hist_sim_data + wind_numb+name_charge[k]+sub_numb;
                ptname_rec_data[i][k][sub_sample] =ptname_hist_rec_data+ wind_numb+name_charge[k]+sub_numb;
                fHistpt_sim_data[i][k][sub_sample] = new TH1D(ptname_sim_data[i][k][sub_sample], ptname_sim_data[i][k][sub_sample], 200,0, 2);
                fHistpt_rec_data[i][k][sub_sample] = new TH1D(ptname_rec_data[i][k][sub_sample], ptname_rec_data[i][k][sub_sample], 200,0, 2);
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
    cout << "n_unf_iterations_1d = "<<n_unf_iterations_1d<<endl;
    
    int n_ev_in_subsample;
    // start filling data histograms (bootstrap method to calculate errors, n_subsamples = 1 corresponds to the full statistics to obtain the value, other n_subsamples to calculate errors)
    cout<<"start data"<<endl;
    
    for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
    {
        cout<<"subsample #"<<sub_sample<<endl;
        n_ev_in_subsample=0;
    
        while(n_ev_in_subsample<nEv_data)
        {
            (*particles_pt_in_event_sim).clear();
            (*particles_ptPlus_in_event_sim).clear();
            (*particles_ptMinus_in_event_sim).clear();
            
            if(n_subsamples==1)
            {
                T->GetEntry(n_ev_in_subsample);
            }
            else
            {
                T->GetEntry(rand.Integer(nEv_data));
            }
            
            // filling pure (sim) and biased (rec) data histograms
            for(int i = 0; i<NumbOne_winds;i++)
            {
                for(long unsigned int j=0;j<(*particles_pt_in_event_sim)[i].size();j++)
                {
                    fHistpt_sim_data[i][0][sub_sample]->Fill((*particles_pt_in_event_sim)[i].at(j));
                }
                
                for(long unsigned int j=0;j<(*particles_pt_in_event_rec)[i].size();j++)
                {
                    fHistpt_rec_data[i][0][sub_sample]->Fill((*particles_pt_in_event_rec)[i].at(j));
                }
                
                for(long unsigned int j=0;j<(*particles_ptPlus_in_event_sim)[i].size();j++)
                {
                    fHistpt_sim_data[i][1][sub_sample]->Fill((*particles_ptPlus_in_event_sim)[i].at(j));
                }
                
                for(long unsigned int j=0;j<(*particles_ptPlus_in_event_rec)[i].size();j++)
                {
                    fHistpt_rec_data[i][1][sub_sample]->Fill((*particles_ptPlus_in_event_rec)[i].at(j));
                }
                
                for(long unsigned int j=0;j<(*particles_ptMinus_in_event_sim)[i].size();j++)
                {
                    fHistpt_sim_data[i][2][sub_sample]->Fill((*particles_ptMinus_in_event_sim)[i].at(j));
                }
                
                for(long unsigned int j=0;j<(*particles_ptMinus_in_event_rec)[i].size();j++)
                {
                    fHistpt_rec_data[i][2][sub_sample]->Fill((*particles_ptMinus_in_event_rec)[i].at(j));
                }
            }
            
            (*particles_pt_in_event_sim).clear();
            (*particles_ptPlus_in_event_sim).clear();
            (*particles_ptMinus_in_event_sim).clear();
            
            n_ev_in_subsample++;
        }// n_ev_in_subsample is enaugh
    }// end of subsamples
    
    // create response matrix based on rec and sim histograms
    RooUnfoldResponse* responsept[NumbOne_winds][charge_combo_pt];
    for(int k = 0; k<NumbOne_winds;k++)
    {
        for(int charge = 0; charge<charge_combo_pt;charge++)
        {
            responsept[k][charge] = new RooUnfoldResponse(fHistpt_rec_MC[k][charge], fHistpt_sim_MC[k][charge]);
        }
    }

    cout<<"start feeling RM"<<endl;
    // filling RM created previously with rec and sim MC
    for (int iEvent = nEv_data; iEvent < entries; iEvent++)
    {
        T->GetEntry(iEvent);
        
        for(int i = 0; i<NumbOne_winds;i++)
        {
            for(long unsigned int j=0;j<(*particles_pt_in_event_sim)[i].size();j++)
            {
                fHistpt_sim_MC[i][0]->Fill((*particles_pt_in_event_sim)[i].at(j));
            }
            for(long unsigned int j=0;j<(*particles_pt_in_event_rec)[i].size();j++)
            {
                fHistpt_rec_MC[i][0]->Fill((*particles_pt_in_event_rec)[i].at(j));
            }
            for(long unsigned int j=0;j<(*particles_ptPlus_in_event_sim)[i].size();j++)
            {
                fHistpt_sim_MC[i][1]->Fill((*particles_ptPlus_in_event_sim)[i].at(j));
            }
            for(long unsigned int j=0;j<(*particles_ptPlus_in_event_rec)[i].size();j++)
            {
                fHistpt_rec_MC[i][1]->Fill((*particles_ptPlus_in_event_rec)[i].at(j));
            }
            for(long unsigned int j=0;j<(*particles_ptMinus_in_event_sim)[i].size();j++)
            {
                fHistpt_sim_MC[i][2]->Fill((*particles_ptMinus_in_event_sim)[i].at(j));
            }
            for(long unsigned int j=0;j<(*particles_ptMinus_in_event_rec)[i].size();j++)
            {
                fHistpt_rec_MC[i][2]->Fill((*particles_ptMinus_in_event_rec)[i].at(j));
            }
     
            // if value of the variable is equal to -1 in REC that means the particle was lost by rec track cuts -> we should put it in "miss" to normalize histogram, if value of the variable is equal to -1 in SIM that means the REC particle is fake -> we should put it in "fake" to normalize histogram
            for(int j=0;j<(*particles_pt_in_event_sim)[i].size();j++)
            {
                if((*particles_pt_in_event_rec)[i].at(j) < 0)
                {
                    responsept[i][0]->Miss((*particles_pt_in_event_sim)[i].at(j));
                }
                else if((*particles_pt_in_event_sim)[i].at(j) <0)
                {
                    responsept[i][0]->Fake((*particles_pt_in_event_rec)[i].at(j));
                }
                else
                {
                    responsept[i][0]->Fill((*particles_pt_in_event_rec)[i].at(j),(*particles_pt_in_event_sim)[i].at(j));
                    if(i==0)
                    {
                    hRM->Fill((*particles_pt_in_event_rec)[0].at(j),(*particles_pt_in_event_sim)[0].at(j));
                    }
                }
            }
            
            
            for(int j=0;j<(*particles_ptPlus_in_event_sim)[i].size();j++)
            {
                if((*particles_ptPlus_in_event_rec)[i].at(j) < 0)
                {
                    responsept[i][1]->Miss((*particles_ptPlus_in_event_sim)[i].at(j));
                }
                else if((*particles_ptPlus_in_event_sim)[i].at(j) <0)
                {
                    responsept[i][1]->Fake((*particles_ptPlus_in_event_rec)[i].at(j));
                }
                else
                {
                    responsept[i][1]->Fill((*particles_ptPlus_in_event_rec)[i].at(j),(*particles_ptPlus_in_event_sim)[i].at(j));
                }
            }
            
            for(int j=0;j<(*particles_ptMinus_in_event_sim)[i].size();j++)
            {
                if((*particles_ptMinus_in_event_rec)[i].at(j) < 0)
                {
                    responsept[i][2]->Miss((*particles_ptMinus_in_event_sim)[i].at(j));
                }
                else if((*particles_ptMinus_in_event_sim)[i].at(j) <0)
                {
                    responsept[i][2]->Fake((*particles_ptMinus_in_event_rec)[i].at(j));
                }
                else
                {
                    responsept[i][2]->Fill((*particles_ptMinus_in_event_rec)[i].at(j),(*particles_ptMinus_in_event_sim)[i].at(j));
                }
            }
        } // end of window loop
        
        (*particles_pt_in_event_sim).clear();
        (*particles_ptPlus_in_event_sim).clear();
        (*particles_ptMinus_in_event_sim).clear();
        
        (*particles_pt_in_event_rec).clear();
        (*particles_ptPlus_in_event_rec).clear();
        (*particles_ptMinus_in_event_rec).clear();
    } // end of MC events
    
    cout<<"creating RM"<<endl;
    // create unfold function (here you can change RooUnfold method: bin-by-bin, RooUnfoldBayes, RooUnfoldSvd) using RM, biased rec data to unfold and the number of iterations in the case of IBU
    RooUnfoldBayes*   unfoldpt[NumbOne_winds][charge_combo_pt][n_subsamples];
    
    for(int m = 0; m<NumbOne_winds;m++)
    {
        for(int charge = 0; charge<charge_combo_pt;charge++)
        {
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                unfoldpt[m][charge][sub_sample] = new RooUnfoldBayes(responsept[m][charge], fHistpt_rec_data[m][charge][sub_sample], n_unf_iterations_1d);
            }
        }
    }
    
     // create histogram for unfolded data, perform unfolding and rename unfolded histo
    TH1D* hRecopt[NumbOne_winds][charge_combo_pt][n_subsamples];
    TString hist_pt_name_unf = "unfolded_pt";
    TString pt_name_unf[NumbOne_winds][charge_combo_pt][n_subsamples];
    cout<<"pt unfolding"<<endl;
    
    for(int j = 0; j<NumbOne_winds;j++)
    {
        sprintf(wind_numb, "_%i", j);
        for(int charge = 0; charge<charge_combo_pt;charge++)
        {
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                hRecopt[j][charge][sub_sample]=(TH1D*) unfoldpt[j][charge][sub_sample]->Hreco();
                sprintf(sub_numb, "_%i", sub_sample);
                pt_name_unf[j][charge][sub_sample] = hist_pt_name_unf+wind_numb+name_charge[charge]+sub_numb;
                hRecopt[j][charge][sub_sample]-> SetName(pt_name_unf[j][charge][sub_sample]);
            }
        }
    }
    
    // create output file and fill it with pure, biased and unfolded data
    TString end = ".root";
    TString outpath = "Unf_pt_"; // put here the path where to save output root file
    TString Sboot;
    TString Sunf_1d;
    
    Sboot = to_string(n_subsamples);
    Sunf_1d = to_string(n_unf_iterations_1d);
    cout<<"unfold objects created"<<endl;
    TString filenamept;
    filenamept = outpath+Sboot+"_"+Sunf_1d+end;
    
    cout<<"file writing "<<filenamept<<endl;

    TFile *fFilept = new TFile(filenamept, "RECREATE");
    for(int j = 0; j<NumbOne_winds;j++)
    {
        for(int charge = 0; charge<charge_combo_pt;charge++)
        {
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                fHistpt_sim_data[j][charge][sub_sample]->Write();
                fHistpt_rec_data[j][charge][sub_sample]->Write();
                hRecopt[j][charge][sub_sample]->Write();
            }
        }
    }
    fFilept->Close();

    cout<<"job is done"<<endl;
    
    return 0;
}
