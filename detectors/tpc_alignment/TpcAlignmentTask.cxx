#include "TpcAlignmentTask.h"
#include "Runner.h"
#include "FileHelper.h"
#include "Enums.h"
#include <Rtypes.h>
#include <string>

using namespace std;
using namespace TpcAlignment;

ClassImp(TpcAlignmentTask);

TpcAlignmentTask::TpcAlignmentTask() : 
FairTask("TpcAlignmentTask", 1) 
{}

TpcAlignmentTask::~TpcAlignmentTask() 
{}

InitStatus TpcAlignmentTask::Init()
{
   return kSUCCESS;
}

void TpcAlignmentTask::Exec(Option_t *opt)
{
	int precision{ 100 };
	const int vNumberOfCalibrationIterations{ 6 };
	Runner R;

	string vRaysFile = FileHelper::BuildFilePath(Solution::release, Direction::input, "LaserRays.txt");
	string matrixA = FileHelper::BuildFilePath(Solution::release, Direction::output, "A.out");
	string coeffMR = FileHelper::BuildFilePath(Solution::release, Direction::output, "MR.out");
	R.LoadModelData(vRaysFile);
	R.LoadCorrectionMatrix({}, false);
	R.Calibrate(vNumberOfCalibrationIterations);
	R.SaveAMR2Files(matrixA, coeffMR, precision);
}

void TpcAlignmentTask::Finish() {}