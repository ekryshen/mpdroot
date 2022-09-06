#include <iostream>

#include "FairTask.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"

#include "TpcAlignmentTask.h"

void run_task(FairTask *task);

void runTpcAlignment()
{
   // ----    Debug option   -------------------------------------------------
   // gDebug = 0;
   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();
   bool vDebugMode{false};
   // Prepare
   TpcAlignmentTask *task = new TpcAlignmentTask(vDebugMode);
   task->SetNumberOfCalibration(1);
   // Run
   run_task(task);
   // End process
   timer.Stop();
   Double_t rtime = timer.RealTime(), ctime = timer.CpuTime();
   cout << "RealTime=" << rtime << " seconds, CpuTime=" << ctime << " seconds\n";
   cout << "Macro finished successfully." << '\n';
}

void run_task(FairTask *task)
{
   // FairRunSim run;
   task->Exec("");
   // run.Init();
   // run.Run();
}
