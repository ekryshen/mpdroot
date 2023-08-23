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