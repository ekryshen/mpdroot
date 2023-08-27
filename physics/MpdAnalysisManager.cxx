#include <iostream>
#include <fstream>

#include "MpdAnalysisManager.h"
#include "TFile.h"
#include <FairRunAna.h>
ClassImp(MpdAnalysisManager);

MpdAnalysisManager::MpdAnalysisManager(const char *name, int nEvents) : fManagerName(name), fEvents(nEvents) {}
void MpdAnalysisManager::AddTask(MpdAnalysisTask *task)
{
   // Add task to the list of tasks to be executed

   if (!task) {
      std::cout << "MpdAnalysisManager: can not add non-existing taks. Please create it first" << std::endl;
   }

   std::cout << "MpdAnalysisManager: adding taks " << task->GetName() << std::endl;

   fTasks.emplace_back(task);
}
void MpdAnalysisManager::Process()
{
   // Do the job:
   //  Initialise taks UserInit
   //  scan data
   //  call Finish method for each task
   //  write outputs
   int nTasks = 0;

   // Open output files (for trees if needed)
   nTasks = fTasks.size();
   std::cout << "Counted " << nTasks << " tasks to process" << std::endl;
   TFile *output[nTasks];

   nTasks = 0;
   for (MpdAnalysisTask *task : fTasks) {

      fOutFileRoot = task->GetOutputName();
      fOutFileRoot = fOutFileRoot + ".root";
      std::cout << "Opening output file for writing " << fOutFileRoot << std::endl;

      output[nTasks] = new TFile(fOutFileRoot, "recreate");
      nTasks++;
   }

   nTasks = 0;
   for (MpdAnalysisTask *task : fTasks) {
      task->UserInit();
      nTasks++;
   }

   if (!CreateChain()) {
      std::cout << "Failed to create input TChain, do nothing" << std::endl;
      return;
   }

   // create geometry and tracking
   // Read geometry field etc
   fChain->GetFile()->Get("FairGeoParSet");
   FairRunAna *ana = new FairRunAna();
   // To propagate tracks
   fKF = MpdKalmanFilter::Instance();
   fKF->Init();

   // Scan data
   int n = fChain->GetEntries();
   if (fEvents != -1 && fEvents < n) n = fEvents;
   std::cout << "MpdAnalysisManager: number of events to process " << n << std::endl;

   for (int iEvent = 0; iEvent < n; iEvent++) {
      if (iEvent % 100 == 0) {
         std::cout << "MpdAnalysisManager: processing event " << iEvent << " of " << n << std::endl;
      }
      fEvent.Clear();
      if (fChain->GetEvent(iEvent)) {
         for (MpdAnalysisTask *task : fTasks) {
            task->ProcessEvent(fEvent);
         }
      } // not broken
   }

   for (MpdAnalysisTask *task : fTasks) {
      task->Finish();
   }

   nTasks = 0;
   for (MpdAnalysisTask *task : fTasks) {

      fOutFileRoot = task->GetOutputName();
      fOutFileRoot = fOutFileRoot + ".root";
      std::cout << "Writing output to file " << fOutFileRoot << std::endl;

      output[nTasks]->cd();
      task->GetOutput()->Write();
      output[nTasks]->Close();
      nTasks++;
   }
}
//=======================================
bool MpdAnalysisManager::CreateChain()
{
   // Create chain and attach branches

   // Create TChain
   fChain = new TChain("mpdsim");
   std::ifstream ifs(fInputFiles);
   std::string   fname;
   if (ifs.is_open()) {
      while (ifs >> fname) {
         std::cout << "Adding to chain: " << fname << std::endl;
         fChain->AddFile(fname.data());
      }
   } else {
      std::cout << "Error opening file " << fInputFiles << std::endl;
      return false;
   }

   ifs.close();

   // Setup branches
   fChain->SetBranchStatus("*", 0); // disable all branches
   if (fBranchList == "*" || fBranchList.Contains("EventHeader")) {
      fChain->SetBranchStatus("EventHeader*", 1);
      fChain->SetBranchAddress("EventHeader.", &(fEvent.fEventHeader));
   }
   if (fBranchList == "*" || fBranchList.Contains("TpcKalmanTrack")) {
      fChain->SetBranchStatus("TpcKalmanTrack*", 1);
      fChain->SetBranchAddress("TpcKalmanTrack", &(fEvent.fTPCKalmanTrack));
   }
   if (fBranchList == "*" || fBranchList.Contains("Vertex")) {
      fChain->SetBranchStatus("Vertex*", 1);
      fChain->SetBranchAddress("Vertex", &(fEvent.fVertex));
   }
   if (fBranchList == "*" || fBranchList.Contains("FfdHit")) {
      fChain->SetBranchStatus("FfdHit*", 1);
      fChain->SetBranchAddress("FfdHit", &(fEvent.fFfdHit));
   }
   if (fBranchList == "*" || fBranchList.Contains("TOFHit")) {
      fChain->SetBranchStatus("TOFHit*", 1);
      fChain->SetBranchAddress("TOFHit", &(fEvent.fTOFHit));
   }
   if (fBranchList == "*" || fBranchList.Contains("TOFMatching")) {
      fChain->SetBranchStatus("TOFMatching*", 1);
      fChain->SetBranchAddress("TOFMatching", &(fEvent.fTOFMatching));
   }
   if (fBranchList == "*" || fBranchList.Contains("EmcDigit")) {
      fChain->SetBranchStatus("EmcDigit*", 1);
      fChain->SetBranchAddress("EmcDigit", &(fEvent.fEmcDigit));
   }
   if (fBranchList == "*" || fBranchList.Contains("EmcCluster")) {
      fChain->SetBranchStatus("EmcCluster", 1);
      fChain->SetBranchAddress("EmcCluster", &(fEvent.fEMCCluster));
   }
   if (fBranchList == "*" || fBranchList.Contains("ZdcDigi")) {
      fChain->SetBranchStatus("ZdcDigi*", 1);
      fChain->SetBranchAddress("ZdcDigi", &(fEvent.fZDCDigit));
   }
   if (fBranchList == "*" || fBranchList.Contains("MCEventHeader")) {
      fChain->SetBranchStatus("MCEventHeader*", 1);
      fChain->SetBranchAddress("MCEventHeader.", &(fEvent.fMCEventHeader));
   }
   if (fBranchList == "*" || fBranchList.Contains("MCTrack")) {
      fChain->SetBranchStatus("MCTrack*", 1);
      fChain->SetBranchAddress("MCTrack", &(fEvent.fMCTrack));
   }
   if (fBranchList == "*" || fBranchList.Contains("MPDEvent")) {
      fChain->SetBranchStatus("MPDEvent*", 1);
      fChain->SetBranchAddress("MPDEvent.", &(fEvent.fMPDEvent));
   }

   std::cout << "Chain created" << std::endl;

   return true;
}
