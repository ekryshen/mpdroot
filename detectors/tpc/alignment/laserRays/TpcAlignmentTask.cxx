#include <iostream>

#include <Rtypes.h>

#include "Alignment.h"
#include "Debug.h"
#include "Run.h"
#include "TpcAlignmentTask.h"

using namespace TpcAlignment;

ClassImp(TpcAlignmentTask);

TpcAlignmentTask::TpcAlignmentTask() : FairTask("TpcAlignmentTask", 1) {}

TpcAlignmentTask::TpcAlignmentTask(bool debugMode) : FairTask("TpcAlignmentTask", 1), vDebugMode(debugMode) {}

TpcAlignmentTask::~TpcAlignmentTask() {}

InitStatus TpcAlignmentTask::Init()
{
   return kSUCCESS;
}

void TpcAlignmentTask::SetNumberOfCalibration(int val)
{
   vNumberOfCalibration = val;
}

void TpcAlignmentTask::Exec(Option_t *opt)
{
   Alignment vAlignment(Debug(false));
   Run       R(vAlignment, Debug(false));

   std::cout << "Load A\n";
   vAlignment.LoadA();

   for (int i = 0; i < vNumberOfCalibration; ++i) {
      std::cout << "Calibration\n";
      R.Calibration("testA.inp");
   }
   std::cout << " Measurement\n";
   R.Measurement("testA.inp", "testA.out");
}

void TpcAlignmentTask::Finish() {}
