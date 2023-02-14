#include "TpcAlignmentTaskLaserRays.h"

#include "Runner.h"
#include "FileHelper.h"
#include "Enums.h"
#include <Rtypes.h>
#include <string>

using namespace TpcAlignmentLaserRays;

ClassImp(TpcAlignmentTaskLaserRays);

TpcAlignmentTaskLaserRays::TpcAlignmentTaskLaserRays() : FairTask("TpcAlignmentTaskLaserRays", 1) {}

TpcAlignmentTaskLaserRays::~TpcAlignmentTaskLaserRays() {}

InitStatus TpcAlignmentTaskLaserRays::Init()
{
   return kSUCCESS;
}

void TpcAlignmentTaskLaserRays::Exec(Option_t *opt)
{
   int       precision{100};
   const int vNumberOfCalibrationIterations{6};
   Runner    R;

   std::string vRaysFile = FileHelper::BuildFilePath(Solution::release, Direction::input, "LaserRays.txt");
   std::string matrixA   = FileHelper::BuildFilePath(Solution::release, Direction::output, "A.out");
   std::string coeffMR   = FileHelper::BuildFilePath(Solution::release, Direction::output, "MR.out");
   R.LoadModelData(vRaysFile);
   R.LoadCorrectionMatrix({}, false);
   R.Calibrate(vNumberOfCalibrationIterations);
   R.SaveAMR2Files(matrixA, coeffMR, precision);
}

void TpcAlignmentTaskLaserRays::Finish() {}
