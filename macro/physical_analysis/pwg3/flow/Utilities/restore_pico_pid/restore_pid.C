#include "../utility.h"
#include "PicoDstPid.cxx"

void restore_pid(TString inFileName, TString outFileName, Double_t energy=11.)
{
  
  Double_t sigM = 4.0, sigE = 4.0, koef = 1.; // n-sigma bands for PID selection
  TString Generator = "URQMD", Tracking = "CF";
  
  TFile * inFile = new TFile(inFileName.Data(),"read");
  TTree * inTree = (TTree*) inFile->Get("mpdsim_reduced");
  
  Long64_t nentries = inTree->GetEntries();

  Long64_t  n_tracks_mpd;
  Float_t   pid_tpc_prob_pion_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
  Float_t   pid_tpc_prob_kaon_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
  Float_t   pid_tpc_prob_proton_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
  Float_t   p_mpd[_MAX_TRACKS];
  Float_t   dEdx_tpc_mpd[_MAX_TRACKS];
  Float_t   tof_mass2_mpd[_MAX_TRACKS];
  Float_t   signed_pt_mpd[_MAX_TRACKS];
  
  int pidFlag = 2, charge = 0;
  
  inTree -> SetBranchAddress("n_tracks_mpd",            &n_tracks_mpd);
  inTree -> SetBranchAddress("pid_tpc_prob_pion_mpd",   pid_tpc_prob_pion_mpd);
  inTree -> SetBranchAddress("pid_tpc_prob_kaon_mpd",   pid_tpc_prob_kaon_mpd);
  inTree -> SetBranchAddress("pid_tpc_prob_proton_mpd", pid_tpc_prob_proton_mpd);
  inTree -> SetBranchAddress("p_mpd",                   p_mpd);
  inTree -> SetBranchAddress("dEdx_tpc_mpd",            dEdx_tpc_mpd);
  inTree -> SetBranchAddress("tof_mass2_mpd",           tof_mass2_mpd);
  inTree -> SetBranchAddress("signed_pt_mpd",           signed_pt_mpd);
  
  TFile * outFile = new TFile(outFileName.Data(),"recreate");
  TTree * outTree = inTree->CloneTree(0);
  
  MpdPid *PID = new MpdPid(sigM, sigE, energy, koef, Generator, Tracking);
  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++){
    if (iEntry % 1000 == 0) std::cout << "Event: " << iEntry << std::endl;
    inTree->GetEntry(iEntry);
    for (Long64_t iTrack=0;iTrack<n_tracks_mpd;iTrack++){
      
      if (signed_pt_mpd[iTrack] > 0) charge = 1;
      else                           charge = -1; 
      
      if ( PID -> FillProbs(p_mpd[iTrack], dEdx_tpc_mpd[iTrack], tof_mass2_mpd[iTrack], charge) ){

        pid_tpc_prob_pion_mpd[iTrack]   = PID-> GetProbPi();
        pid_tpc_prob_kaon_mpd[iTrack]   = PID-> GetProbKa();
        pid_tpc_prob_proton_mpd[iTrack] = PID-> GetProbPr();
        
      }else{
        if ( PID -> FillProbs(p_mpd[iTrack], dEdx_tpc_mpd[iTrack], charge) ){
          pid_tpc_prob_pion_mpd[iTrack]   = PID-> GetProbPi();
          pid_tpc_prob_kaon_mpd[iTrack]   = PID-> GetProbKa();
          pid_tpc_prob_proton_mpd[iTrack] = PID-> GetProbPr();
        }else{
          pid_tpc_prob_pion_mpd[iTrack]   = -1.;
          pid_tpc_prob_kaon_mpd[iTrack]   = -1.;
          pid_tpc_prob_proton_mpd[iTrack] = -1.;
        }
      }
    }
    outTree->Fill();
  }
  
  outTree->AutoSave();
  delete inFile;
  delete outFile;
  
}
