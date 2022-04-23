#ifndef MPDANALYSISMANAGER_H
#define MPDANALYSISMANAGER_H

#include <vector>
#include "MpdEvent.h"
#include "MpdAnalysisTask.h"
#include "TChain.h"
#include "MpdKalmanFilter.h"

class MpdAnalysisManager {

public:
   MpdAnalysisManager();
   MpdAnalysisManager(const char *name);
   ~MpdAnalysisManager() {} // Destructor

   // Add task to perform analyses
   // First added tasks will be performed firs (re-calibration etc. before real analysis)
   void AddTask(MpdAnalysisTask *task);

   // Process inputFileList with tasks added with AddTask method
   void Process();

   // txt file with list of input files
   void InputFileList(const char *filename) { fInputFiles = filename; }

   // Provide coma separated list of branches to be read for particular analyses
   //  * will read all branches (slower)
   void ReadBranches(const char *branchlist = "*") { fBranchList = branchlist; }

   // Where to write output file
   void SetOutput(const char *outputFileName = "histos.root") { fOutFile = outputFileName; }

protected:
   bool CreateChain();

private:
   TString                        fInputFiles;
   TChain                        *fChain = nullptr;
   std::vector<MpdAnalysisTask *> fTasks;
   TString                        fOutFile    = "";
   TString                        fBranchList = "";
   MpdAnalysisEvent               fEvent;
   MpdKalmanFilter               *fKF = nullptr;

   ClassDef(MpdAnalysisManager, 0);
};

#endif
