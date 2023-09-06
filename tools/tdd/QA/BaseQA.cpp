#include "BaseQA.h"
#include "FairLogger.h"
ClassImp(BaseQA);

void BaseQA::WriteBaseQA(TString suffix)
{
   TString outputFilename("BaseQA.root");
   if (suffix.Length() > 0) outputFilename = TString("BaseQA_") + suffix + TString(".root");
   LOG(info) << " Writing QA file: " << outputFilename;

   TFile outputFile(outputFilename, "recreate");
   TTree tree("QA", "QA data tree");

   TClonesArray *mcTracks = new TClonesArray("MpdMCTrack");

   TClonesArray &mcTracks_ref = *mcTracks;
   tree.Branch("MC_Tracks", &mcTracks, 256000, 1);

   TClonesArray *tpcTracks     = new TClonesArray("MpdTpcKalmanTrack");
   TClonesArray &tpcTracks_ref = *tpcTracks;
   tree.Branch("TPC_Tracks", &tpcTracks, 256000, 1);

   int eventCount = eventNumber.size();

   for (int event = 0; event < eventCount; ++event) {
      mcTracks_ref.Clear();
      tpcTracks_ref.Clear();
      mcTracks_ref  = *(eventMCTracksArray.at(event));
      tpcTracks_ref = *(eventTpcTracksArray.at(event));
      tree.Fill();
   }

   tree.Write();
   outputFile.Close();
}

void BaseQA::ReadBaseQA(TString suffix, TString dir)
{
   TString inputFilename("");
   if (dir.Length() > 0) inputFilename = dir + TString("/");
   if (suffix.Length() > 0)
      inputFilename += TString("BaseQA_") + suffix + TString(".root");
   else
      inputFilename += TString("BaseQA.root");
   LOG(info) << " Reading QA file: " << inputFilename;

   TFile *f    = new TFile(inputFilename);
   TTree *tree = (TTree *)f->Get("QA");

   TClonesArray *mcTracks = new TClonesArray("MpdMCTrack");
   tree->GetBranch("MC_Tracks")->SetAutoDelete(kFALSE);
   tree->SetBranchAddress("MC_Tracks", &mcTracks);

   TClonesArray *tpcTracks = new TClonesArray("MpdTpcKalmanTrack");
   tree->GetBranch("TPC_Tracks")->SetAutoDelete(kFALSE);
   tree->SetBranchAddress("TPC_Tracks", &tpcTracks);

   int nEvents = tree->GetEntries();
   for (int i = 0; i < nEvents; ++i) {
      mcTracks->Clear();
      tpcTracks->Clear();

      tree->GetEntry(i);
      eventMCTracksArray.push_back(new TClonesArray());
      eventMCTracksArray.back() = (TClonesArray *)mcTracks->Clone();
      eventTpcTracksArray.push_back(new TClonesArray());
      eventTpcTracksArray.back() = (TClonesArray *)tpcTracks->Clone();
   }
}