#include <iostream>
#include <fstream>

#include "MpdAnalysisManager.h"
#include "TFile.h"
#include <FairRunAna.h>
ClassImp(MpdAnalysisManager);

MpdAnalysisManager::MpdAnalysisManager(const char *name) {}
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

   for (MpdAnalysisTask *task : fTasks) {
      task->UserInit();
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
   for (int iEvent = 0; iEvent < n; iEvent++) {
      if (iEvent % (n / 10) == 0) {
         std::cout << "MpdAnalysisManager: processing event " << iEvent << " of " << n << std::endl;
      }
      fEvent.Clear();
      fChain->GetEvent(iEvent);
      for (MpdAnalysisTask *task : fTasks) {
         task->ProcessEvent(fEvent);
      }
   }

   for (MpdAnalysisTask *task : fTasks) {
      task->Finish();
   }

   TFile output(fOutFile, "recreate");
   for (MpdAnalysisTask *task : fTasks) {
      task->GetOutput()->Write(task->GetOutputName());
   }
   output.Close();
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
   if (fBranchList == "*" || fBranchList.Contains("EventHeader")) {
      fChain->SetBranchAddress("EventHeader.", &(fEvent.fEventHeader));
   }
   if (fBranchList == "*" || fBranchList.Contains("Vertex")) {
      fChain->SetBranchAddress("Vertex", &(fEvent.fVertex));
   }
   if (fBranchList == "*" || fBranchList.Contains("MPDEvent")) {
      fChain->SetBranchAddress("MPDEvent.", &(fEvent.fMPDEvent));
      std::cout << "Reading MPDEvent=" << fEvent.fMPDEvent << std::endl;
   }
   if (fBranchList == "*" || fBranchList.Contains("MCEventHeader")) {
      fChain->SetBranchAddress("MCEventHeader.", &(fEvent.fMCEventHeader));
   }
   if (fBranchList == "*" || fBranchList.Contains("MCTrack")) {
      fChain->SetBranchAddress("MCTrack", &(fEvent.fMCTrack));
   }
   if (fBranchList == "*" || fBranchList.Contains("EmcCluster")) {
      fChain->SetBranchAddress("EmcCluster", &(fEvent.fEMCCluster));
   }
   if (fBranchList == "*" || fBranchList.Contains("TpcKalmanTrack")) {
      fChain->SetBranchAddress("TpcKalmanTrack", &(fEvent.fTPCKalmanTrack));
   }
   if (fBranchList == "*" || fBranchList.Contains("TOFHit")) {
      fChain->SetBranchAddress("TOFHit", &(fEvent.fTOFHit));
   }
   if (fBranchList == "*" || fBranchList.Contains("TOFMatching")) {
      fChain->SetBranchAddress("TOFMatching", &(fEvent.fTOFMatching));
   }
   if (fBranchList == "*" || fBranchList.Contains("ZdcDigi")) {
      fChain->SetBranchAddress("ZdcDigi", &(fEvent.fZDCDigit));
   }
   if (fBranchList == "*" || fBranchList.Contains("ELossZdc1Value")) {
      fChain->SetBranchAddress("ELossZdc1Value", &(fEvent.fZDCEloss1Value));
   }
   if (fBranchList == "*" || fBranchList.Contains("ELossZdc2Value")) {
      fChain->SetBranchAddress("ELossZdc2Value", &(fEvent.fZDCEloss1Value));
   }
   if (fBranchList == "*" || fBranchList.Contains("ELossZdc1Histo")) {
      fChain->SetBranchAddress("ELossZdc1Histo", &(fEvent.fZDCEloss1Value));
   }
   if (fBranchList == "*" || fBranchList.Contains("ELossZdc2Histo")) {
      fChain->SetBranchAddress("ELossZdc2Histo", &(fEvent.fZDCEloss1Value));
   }
   std::cout << "Chain created" << std::endl;

   return true;
}