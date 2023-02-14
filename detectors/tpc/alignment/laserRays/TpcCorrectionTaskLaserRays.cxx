#include "TpcCorrectionTaskLaserRays.h"
#include "Runner.h"
#include "FileHelper.h"
#include "Enums.h"
#include <Rtypes.h>
#include <string>

using namespace TpcAlignmentLaserRays;

ClassImp(TpcCorrectionTaskLaserRays);

TpcCorrectionTaskLaserRays::TpcCorrectionTaskLaserRays() : FairTask("TpcCorrectionTaskLaserRays", 1) {}

TpcCorrectionTaskLaserRays::~TpcCorrectionTaskLaserRays() {}

InitStatus TpcCorrectionTaskLaserRays::Init()
{
   return kSUCCESS;
}

void TpcCorrectionTaskLaserRays::Exec(Option_t *opt)
{
   Runner R;

   std::string filepathA = FileHelper::BuildFilePath(Solution::release, Direction::input, "CorrectionMatrix.dat");
   R.LoadCorrectionMatrix(filepathA);

   std::string inputTracks = FileHelper::BuildFilePath(Solution::release, Direction::input, "tracks.inp");
   std::string outputTracs = FileHelper::BuildFilePath(Solution::release, Direction::output, "tracks.out");
   R.CorrectTracks(inputTracks, outputTracs);
}

void TpcCorrectionTaskLaserRays::Finish() {}
