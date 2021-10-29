// ==================================================================================================================================================================================
// This code performes 2D-unfolding of PT-N distribution (the correlation between the scalar sum of event transverse momentum and event multiplicity) using the pre-prepared tree that contains sim and rec info about the events to build detector Response Matrix (RM).
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

//***************************************************************************************************************************
// to run locally do
// root -l
// gSystem->Load("/*put_your_path*/RooUnfold/libRooUnfold");
// .L final_unfolding_PTN.cc+
// .x final_unfolding_PTN( 100, 50 , 500)
//***************************************************************************************************************************

// before doing unfolding the RooUnfold package should be installed from https://gitlab.cern.ch/RooUnfold/RooUnfold

#include "/*put_your_path*/RooUnfold/src/RooUnfoldResponse.h"
#include "/*put_your_path*/RooUnfold/src/RooUnfoldBayes.h"

using namespace std;
using namespace TMath;

// n_subsamples - defines the number of bootstrap samples. n_subsamples = 1 will be used to calculate the central value without bootstrap, other values of n_subsamples will be used to calculate statistical uncertainty. n_subsamples is used only for data. No bootstrap is applied to the RM creation: all bootstrap data histo are unfolded by the same RM.

//n_unf_iterations - defines the number of unfolding iterations in the iterative bayesian unfolding. Be aware that too big number of iterations might start to increase the error of unfolding from some point.

// n_bins - defines the number of bins for PT histogram

int final_unfolding_PTN(int n_subsamples, int n_unf_iterations_2d, int n_bins)
{
    int NumbOne_winds = 10; // number of pseudorapidity intervals
    int charge_combo = 5;  // number of charge combinations: all, pos, neg, posneg, negpos
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
    int multOne_rec[NumbOne_winds], multOnePlus_rec[NumbOne_winds], multOneMinus_rec[NumbOne_winds], multOne_sim[NumbOne_winds],  multOnePlus_sim[NumbOne_winds], multOneMinus_sim[NumbOne_winds];
    
    Double_t PTOne_rec[NumbOne_winds], PT2One_rec[NumbOne_winds], PTOnePlus_rec[NumbOne_winds],PT2OnePlus_rec[NumbOne_winds],PTOneMinus_rec[NumbOne_winds], PT2OneMinus_rec[NumbOne_winds], PTOne_sim[NumbOne_winds], PT2One_sim[NumbOne_winds],PTOnePlus_sim[NumbOne_winds], PT2OnePlus_sim[NumbOne_winds], PTOneMinus_sim[NumbOne_winds], PT2OneMinus_sim[NumbOne_winds];
    
    TRandom3 rand;
    rand.SetSeed(0);
    
    // access tree with sim and rec data
    TTree *T = (TTree*)file_to_an->Get(name);
    T->SetBranchAddress("multOne_rec", multOne_rec);
    T->SetBranchAddress("PTOne_rec", PTOne_rec);
    T->SetBranchAddress("multOnePlus_rec", multOnePlus_rec);
    T->SetBranchAddress("multOneMinus_rec", multOneMinus_rec);
    T->SetBranchAddress("PTOnePlus_rec", PTOnePlus_rec);
    T->SetBranchAddress("PTOneMinus_rec", PTOneMinus_rec);
    
    T->SetBranchAddress("multOne_sim", multOne_sim);
    T->SetBranchAddress("PTOne_sim", PTOne_sim);
    T->SetBranchAddress("multOnePlus_sim", multOnePlus_sim);
    T->SetBranchAddress("multOneMinus_sim", multOneMinus_sim);
    T->SetBranchAddress("PTOnePlus_sim", PTOnePlus_sim);
    T->SetBranchAddress("PTOneMinus_sim", PTOneMinus_sim);
    
    // create histograms to calculate value for RM and value+error from subsamples for data
    TH2D*                   fHistPTN_sim_MC[NumbOne_winds][charge_combo];
    TH2D*                   fHistPTN_sim_data[NumbOne_winds][charge_combo][n_subsamples];
    
    TH2D*                   fHistPTN_rec_MC[NumbOne_winds][charge_combo];
    TH2D*                   fHistPTN_rec_data[NumbOne_winds][charge_combo][n_subsamples];
    
    TString PTNname_hist_sim_MC = "fHistPTN_sim_MC";
    TString PTNname_hist_sim_data = "fHistPTN_sim_data";
    TString PTNname_hist_rec_MC = "fHistPTN_rec_MC";
    TString PTNname_hist_rec_data = "fHistPTN_rec_data";
    
    char wind_numb[50];
    TString PTNname_sim_MC[NumbOne_winds][charge_combo],PTNname_rec_MC[NumbOne_winds][charge_combo],PTNname_sim_data[NumbOne_winds][charge_combo][n_subsamples],PTNname_rec_data[NumbOne_winds][charge_combo][n_subsamples];
    
    TString name_charge[charge_combo];
    name_charge[0] = "_all";
    name_charge[1] = "_pos";
    name_charge[2] = "_neg";
    name_charge[3] = "_PTposNneg";
    name_charge[4] = "_NposPTneg";
    char sub_numb[50];
    
//    int n_bins_all = 502;
//    int pt_bin_all = 1230;
//
//    int n_bins_charges_plus = 292;
//    int pt_bins_charges_plus = 780;
//
//    int n_bins_charges_minus = 242;
//    int pt_bins_charges_minus = 540;
//
//    double bin_all_max = 501.5;
//    double bin_pt_all_max = 246;
//
//    double bin_charge_max_plus = 291.5;
//    double bin_pt_charge_max_plus = 156;
//
//    double bin_charge_max_minus = 241.5;
//    double bin_pt_charge_max_minus = 108;
    
    int n_bins_all = 502;
    int pt_bin_all = n_bins;
    
    int n_bins_charges_plus = 292;
    int pt_bins_charges_plus = n_bins;
    
    int n_bins_charges_minus = 242;
    int pt_bins_charges_minus = n_bins;
    
    double bin_all_max = 501.5;
    double bin_pt_all_max = 246;
    
    double bin_charge_max_plus = 291.5;
    double bin_pt_charge_max_plus = 156;
    
    double bin_charge_max_minus = 241.5;
    double bin_pt_charge_max_minus = 108;
    
    double bin_min = -0.5;
    
    for(int i = 0; i<NumbOne_winds;i++)
    {
        sprintf(wind_numb, "_%i", i);
        
        for(int k = 0;k<charge_combo;k++)
        {
            PTNname_sim_MC[i][k] = PTNname_hist_sim_MC + wind_numb+name_charge[k];
            PTNname_rec_MC[i][k] = PTNname_hist_rec_MC + wind_numb+name_charge[k];
            
            if(k==0)
            {
            fHistPTN_sim_MC[i][0] = new TH2D(PTNname_sim_MC[i][0], PTNname_sim_MC[i][0], n_bins_all, bin_min, bin_all_max, pt_bin_all, 0, bin_pt_all_max);//240, 0, 8
            fHistPTN_rec_MC[i][0] = new TH2D(PTNname_rec_MC[i][0], PTNname_rec_MC[i][0], n_bins_all, bin_min, bin_all_max, pt_bin_all, 0, bin_pt_all_max);
            }
            if(k==1)
            {
                fHistPTN_sim_MC[i][1] = new TH2D(PTNname_sim_MC[i][1], PTNname_sim_MC[i][1], n_bins_charges_plus, bin_min, bin_charge_max_plus, pt_bins_charges_plus, 0, bin_pt_charge_max_plus);//240, 0, 8
                fHistPTN_rec_MC[i][1] = new TH2D(PTNname_rec_MC[i][1], PTNname_rec_MC[i][1], n_bins_charges_plus, bin_min, bin_charge_max_plus, pt_bins_charges_plus, 0, bin_pt_charge_max_plus);
            }
            if(k==2)
            {
                fHistPTN_sim_MC[i][2] = new TH2D(PTNname_sim_MC[i][2], PTNname_sim_MC[i][2], n_bins_charges_minus, bin_min, bin_charge_max_minus, pt_bins_charges_minus, 0, bin_pt_charge_max_minus);//240, 0, 8
                fHistPTN_rec_MC[i][2] = new TH2D(PTNname_rec_MC[i][2], PTNname_rec_MC[i][2], n_bins_charges_minus, bin_min, bin_charge_max_minus, pt_bins_charges_minus, 0, bin_pt_charge_max_minus);
            }
            if(k==3)
            {
                fHistPTN_sim_MC[i][3] = new TH2D(PTNname_sim_MC[i][3], PTNname_sim_MC[i][3], n_bins_charges_minus, bin_min, bin_charge_max_minus, pt_bins_charges_plus, 0, bin_pt_charge_max_plus);//240, 0, 8
                fHistPTN_rec_MC[i][3] = new TH2D(PTNname_rec_MC[i][3], PTNname_rec_MC[i][3], n_bins_charges_minus, bin_min, bin_charge_max_minus, pt_bins_charges_plus, 0, bin_pt_charge_max_plus);
            }
            if(k==4)
            {
                fHistPTN_sim_MC[i][4] = new TH2D(PTNname_sim_MC[i][4], PTNname_sim_MC[i][4], n_bins_charges_plus, bin_min, bin_charge_max_plus, pt_bins_charges_minus, 0, bin_pt_charge_max_minus);//240, 0, 8
                fHistPTN_rec_MC[i][4] = new TH2D(PTNname_rec_MC[i][4], PTNname_rec_MC[i][4], n_bins_charges_plus, bin_min, bin_charge_max_plus, pt_bins_charges_minus, 0, bin_pt_charge_max_minus);
            }
            
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                sprintf(sub_numb, "_%i", sub_sample);
                PTNname_sim_data[i][k][sub_sample] = PTNname_hist_sim_data + wind_numb+name_charge[k]+sub_numb;
                PTNname_rec_data[i][k][sub_sample] = PTNname_hist_rec_data + wind_numb+name_charge[k]+sub_numb;
                
                if(k==0)
                {
                    fHistPTN_sim_data[i][0][sub_sample] = new TH2D(PTNname_sim_data[i][0][sub_sample], PTNname_sim_data[i][0][sub_sample], n_bins_all, bin_min, bin_all_max, pt_bin_all, 0, bin_pt_all_max);//240, 0, 8
                    fHistPTN_rec_data[i][0][sub_sample] = new TH2D(PTNname_rec_data[i][0][sub_sample], PTNname_rec_data[i][0][sub_sample], n_bins_all, bin_min, bin_all_max, pt_bin_all, 0, bin_pt_all_max);
                }
                if(k==1)
                {
                    fHistPTN_sim_data[i][1][sub_sample] = new TH2D(PTNname_sim_data[i][1][sub_sample], PTNname_sim_data[i][1][sub_sample], n_bins_charges_plus, bin_min, bin_charge_max_plus, pt_bins_charges_plus, 0, bin_pt_charge_max_plus);//240, 0, 8
                    fHistPTN_rec_data[i][1][sub_sample] = new TH2D(PTNname_rec_data[i][1][sub_sample], PTNname_rec_data[i][1][sub_sample], n_bins_charges_plus, bin_min, bin_charge_max_plus, pt_bins_charges_plus, 0, bin_pt_charge_max_plus);
                }
                if(k==2)
                {
                    fHistPTN_sim_data[i][2][sub_sample] = new TH2D(PTNname_sim_data[i][2][sub_sample], PTNname_sim_data[i][2][sub_sample], n_bins_charges_minus, bin_min, bin_charge_max_minus, pt_bins_charges_minus, 0, bin_pt_charge_max_minus);//240, 0, 8
                    fHistPTN_rec_data[i][2][sub_sample] = new TH2D(PTNname_rec_data[i][2][sub_sample], PTNname_rec_data[i][2][sub_sample], n_bins_charges_minus, bin_min, bin_charge_max_minus, pt_bins_charges_minus, 0, bin_pt_charge_max_minus);
                }
                if(k==3)
                {
                    fHistPTN_sim_data[i][3][sub_sample] = new TH2D(PTNname_sim_data[i][3][sub_sample], PTNname_sim_data[i][3][sub_sample], n_bins_charges_minus, bin_min, bin_charge_max_minus, pt_bins_charges_plus, 0, bin_pt_charge_max_plus);//240, 0, 8
                    fHistPTN_rec_data[i][3][sub_sample] = new TH2D(PTNname_rec_data[i][3][sub_sample], PTNname_rec_data[i][3][sub_sample], n_bins_charges_minus, bin_min, bin_charge_max_minus, pt_bins_charges_plus, 0, bin_pt_charge_max_plus);
                }
                if(k==4)
                {
                    fHistPTN_sim_data[i][4][sub_sample] = new TH2D(PTNname_sim_data[i][4][sub_sample], PTNname_sim_data[i][4][sub_sample], n_bins_charges_plus, bin_min, bin_charge_max_plus, pt_bins_charges_minus, 0, bin_pt_charge_max_minus);//240, 0, 8
                    fHistPTN_rec_data[i][4][sub_sample] = new TH2D(PTNname_rec_data[i][4][sub_sample], PTNname_rec_data[i][4][sub_sample], n_bins_charges_plus, bin_min, bin_charge_max_plus, pt_bins_charges_minus, 0, bin_pt_charge_max_minus);
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
    cout << "NumbOne_winds = "<<NumbOne_winds<<endl;
    cout<<" "<<endl;
    cout << "n_subsamples = "<<n_subsamples<<endl;
    cout << "n_unf_iterations_2d = "<<n_unf_iterations_2d<<endl;

    int n_ev_in_subsample;
    // start filling data histograms (bootstrap method to calculate errors, n_subsamples = 1 corresponds to the full statistics to obtain the value, other n_subsamples to calculate errors)
    cout<<"start data"<<endl;
    
    for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
    {
        cout<<"subsample #"<<sub_sample<<endl;
        n_ev_in_subsample=0;
        
        while(n_ev_in_subsample<nEv_data)
        {
            for(int i = 0; i<NumbOne_winds;i++)
            {
                multOne_rec[i]=0;PTOne_rec[i]=0; multOnePlus_rec[i]=0; multOneMinus_rec[i]=0; PTOnePlus_rec[i]=0;PTOneMinus_rec[i]=0; multOne_sim[i]=0; PTOne_sim[i]=0; multOnePlus_sim[i]=0; multOneMinus_sim[i]=0; PTOnePlus_sim[i]=0;  PTOneMinus_sim[i]=0;
            }
            
            if(n_subsamples==1)
            {
                T->GetEntry(n_ev_in_subsample);
            }
            else
            {
                T->GetEntry(rand.Integer(nEv_data));
            }
            
            for(int i = 0; i<NumbOne_winds;i++)
            {
                // filling pure (sim) and biased (rec) data histograms
                fHistPTN_rec_data[i][0][sub_sample]->Fill(multOne_rec[i], PTOne_rec[i]);
                fHistPTN_sim_data[i][0][sub_sample]->Fill(multOne_sim[i], PTOne_sim[i]);
                
                fHistPTN_rec_data[i][1][sub_sample]->Fill(multOnePlus_rec[i], PTOnePlus_rec[i]);
                fHistPTN_sim_data[i][1][sub_sample]->Fill(multOnePlus_sim[i], PTOnePlus_sim[i]);
                
                fHistPTN_rec_data[i][2][sub_sample]->Fill(multOneMinus_rec[i], PTOneMinus_rec[i]);
                fHistPTN_sim_data[i][2][sub_sample]->Fill(multOneMinus_sim[i], PTOneMinus_sim[i]);
                
                fHistPTN_rec_data[i][3][sub_sample]->Fill(multOneMinus_rec[i], PTOnePlus_rec[i]);
                fHistPTN_sim_data[i][3][sub_sample]->Fill(multOneMinus_sim[i], PTOnePlus_sim[i]);
                
                fHistPTN_rec_data[i][4][sub_sample]->Fill(multOnePlus_rec[i], PTOneMinus_rec[i]);
                fHistPTN_sim_data[i][4][sub_sample]->Fill(multOnePlus_sim[i], PTOneMinus_sim[i]);
            }
            n_ev_in_subsample++;
        }// n_ev_in_subsample is enaugh
    }// end of subsamples
    
    // create response matrix based on rec and sim histograms
    RooUnfoldResponse* responsePTN[NumbOne_winds][charge_combo];
    for(int k = 0; k<NumbOne_winds;k++)
    {
        for(int charge = 0; charge<charge_combo;charge++)
        {
            responsePTN[k][charge] = new RooUnfoldResponse(fHistPTN_rec_MC[k][charge], fHistPTN_sim_MC[k][charge]);
        }
    }
        
    cout<<"start feeling RM"<<endl;
    // filling RM created previously with rec and sim MC

    for (int iEvent = nEv_data; iEvent < entries; iEvent++)
    {
        for(int i = 0; i<NumbOne_winds;i++)
        {
            multOne_rec[i]=0;PTOne_rec[i]=0; multOnePlus_rec[i]=0; multOneMinus_rec[i]=0; PTOnePlus_rec[i]=0;PTOneMinus_rec[i]=0;
            multOne_sim[i]=0; PTOne_sim[i]=0; multOnePlus_sim[NumbOne_winds]=0; multOneMinus_sim[i]=0; PTOnePlus_sim[i]=0;  PTOneMinus_sim[i]=0;
        }
        
        T->GetEntry(iEvent);
        
        for(int i = 0; i<NumbOne_winds;i++)
        {
            fHistPTN_rec_MC[i][0]->Fill(multOne_rec[i], PTOne_rec[i]);
            fHistPTN_sim_MC[i][0]->Fill(multOne_sim[i], PTOne_sim[i]);
            
            fHistPTN_rec_MC[i][1]->Fill(multOnePlus_rec[i], PTOnePlus_rec[i]);
            fHistPTN_sim_MC[i][1]->Fill(multOnePlus_sim[i], PTOnePlus_sim[i]);
            
            fHistPTN_rec_MC[i][2]->Fill(multOneMinus_rec[i], PTOneMinus_rec[i]);
            fHistPTN_sim_MC[i][2]->Fill(multOneMinus_sim[i], PTOneMinus_sim[i]);
            
            fHistPTN_rec_MC[i][3]->Fill(multOneMinus_rec[i], PTOnePlus_rec[i]);
            fHistPTN_sim_MC[i][3]->Fill(multOneMinus_sim[i], PTOnePlus_sim[i]);
            
            fHistPTN_rec_MC[i][4]->Fill(multOnePlus_rec[i], PTOneMinus_rec[i]);
            fHistPTN_sim_MC[i][4]->Fill(multOnePlus_sim[i], PTOneMinus_sim[i]);

             // if value of the variable is equal to -1 in REC that means the event was lost by rec event cuts -> we should put it in "miss" to normalize histogram
            
            if((multOne_rec[i]<0)&&(PTOne_rec[i]<0))
            {
                responsePTN[i][0]->Miss(multOne_sim[i], PTOne_sim[i]);
            }
            else
            {
                responsePTN[i][0]->Fill(multOne_rec[i], PTOne_rec[i], multOne_sim[i], PTOne_sim[i]);
            }

            if((multOnePlus_rec[i]<0) && (PTOnePlus_rec[i]<0))
            {
              responsePTN[i][1]->Miss(multOnePlus_sim[i],  PTOnePlus_sim[i]);
            }
            else
            {
              responsePTN[i][1]->Fill(multOnePlus_rec[i], PTOnePlus_rec[i], multOnePlus_sim[i],  PTOnePlus_sim[i]);
            }
            
            if((multOneMinus_rec[i]<0)&&(PTOneMinus_rec[i]<0))
            {
                responsePTN[i][2]->Miss(multOneMinus_sim[i],  PTOneMinus_sim[i]);
            }
            else
            {
              responsePTN[i][2]->Fill(multOneMinus_rec[i], PTOneMinus_rec[i], multOneMinus_sim[i],  PTOneMinus_sim[i]);
            }
            
            if((multOneMinus_rec[i]<0)&&(PTOnePlus_rec[i]<0))
            {
             responsePTN[i][3]->Miss(multOneMinus_sim[i],  PTOnePlus_sim[i]);
            }
            else
            {
                responsePTN[i][3]->Fill(multOneMinus_rec[i], PTOnePlus_rec[i], multOneMinus_sim[i],  PTOnePlus_sim[i]);
            }
            
            if((multOnePlus_rec[i]<0)&&(PTOneMinus_rec[i]<0))
            {
              responsePTN[i][4]->Miss(multOnePlus_sim[i],  PTOneMinus_sim[i]);
            }
            else
            {
              responsePTN[i][4]->Fill(multOnePlus_rec[i], PTOneMinus_rec[i], multOnePlus_sim[i],  PTOneMinus_sim[i]);
            }

        } // end of window loop
    } // end of MC events
    
   // create unfold function (here you can change RooUnfold method: bin-by-bin, RooUnfoldBayes, RooUnfoldSvd) using RM, biased rec data to unfold and the number of iterations in the case of IBU
    
    cout<<"unfold objects created"<<endl;
    // create histogram for unfolded data, perform unfolding and rename unfolded histo
    TH2D* hRecoPTN[NumbOne_winds][charge_combo][n_subsamples];
    TString hist_PTN_name_unf = "unfolded_PTN";
    TString PTN_name_unf[NumbOne_winds][charge_combo][n_subsamples];
    
     cout<<"PTN unfolding and file writing"<<endl;
    
    // create output file and fill it with pure, biased and unfolded data
    TString end = ".root";
    TString outpath = "Unf_PTN_"; // put here the path where to save output root file
    TString Sboot;
    TString Sunf_2d;
    TString Sbins;
    Sboot = to_string(n_subsamples);
    Sunf_2d = to_string(n_unf_iterations_2d);
    Sbins = to_string(n_bins);
    
    TString filenamePTN;
    filenamePTN = outpath+Sboot+"_"+Sunf_2d+"_"+Sbins+end;
    
    TFile *fFilePTN = new TFile(filenamePTN, "RECREATE");
    
    RooUnfoldBayes*   unfoldPTN[NumbOne_winds][charge_combo][n_subsamples];
    for(int m = 0; m<NumbOne_winds;m++)
    {
        sprintf(wind_numb, "_%i", m);
        for(int charge = 0; charge<charge_combo;charge++)
        {
            for (int sub_sample = 0;sub_sample < n_subsamples; sub_sample++)
            {
                unfoldPTN[m][charge][sub_sample] = new RooUnfoldBayes(responsePTN[m][charge], fHistPTN_rec_data[m][charge][sub_sample], n_unf_iterations_2d);
     
                hRecoPTN[m][charge][sub_sample]=(TH2D*) unfoldPTN[m][charge][sub_sample]->Hreco();
                sprintf(sub_numb, "_%i", sub_sample);
                PTN_name_unf[m][charge][sub_sample] = hist_PTN_name_unf+wind_numb+name_charge[charge]+sub_numb;
                hRecoPTN[m][charge][sub_sample]-> SetName(PTN_name_unf[m][charge][sub_sample]);
                cout<<"PTN_name_unf["<<m<<"]["<<charge<<"]["<<sub_sample<<"] = "<<PTN_name_unf[m][charge][sub_sample]<<endl;
                
                delete responsePTN[m][charge];
                delete unfoldPTN[m][charge][sub_sample];
                
                fHistPTN_sim_data[m][charge][sub_sample]->Write();
                fHistPTN_rec_data[m][charge][sub_sample]->Write();
                hRecoPTN[m][charge][sub_sample]->Write();
                
                delete fHistPTN_sim_data[m][charge][sub_sample];
                delete fHistPTN_rec_data[m][charge][sub_sample];
                delete hRecoPTN[m][charge][sub_sample];
            }
        }
    }
    
    fFilePTN->Close();
    cout<<"job is done"<<endl;
    
    return 0;
}
