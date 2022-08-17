#include "TpcAlignmentTask.h"
#include "Run.h"
#include "Debug.h"
#include "Alignment.h"
#include <iostream>
#include <Rtypes.h>

using namespace std;
using namespace TpcAlignment;

ClassImp(TpcAlignmentTask);

TpcAlignmentTask::TpcAlignmentTask() : FairTask("TpcAlignmentTask", 1) {}

TpcAlignmentTask::TpcAlignmentTask(bool debugMode)
    : FairTask("TpcAlignmentTask", 1), vDebugMode(debugMode) {}

TpcAlignmentTask::~TpcAlignmentTask() {}

InitStatus TpcAlignmentTask::Init() { return kSUCCESS; }

void TpcAlignmentTask::SetNumberOfCalibration(int val) {
  vNumberOfCalibration = val;
}

void TpcAlignmentTask::Exec(Option_t *opt) {
  Alignment vAlignment(Debug(false));
	Run R(vAlignment, Debug(false));

	cout << "Load A\n";
	vAlignment.LoadA();

	for (int i = 0; i < vNumberOfCalibration; ++i)
	{
		cout << "Calibration\n";
		R.Calibration("testA.inp");
	}
	cout << "Measurement\n";
	R.Measurement("testA.inp", "testA.out");
}

void TpcAlignmentTask::Finish() {}