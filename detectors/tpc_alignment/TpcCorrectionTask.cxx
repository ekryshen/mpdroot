#include "TpcCorrectionTask.h"
#include "Runner.h"
#include "FileHelper.h"
#include "Enums.h"
#include <Rtypes.h>
#include <string>

using namespace std;
using namespace TpcAlignment;

ClassImp(TpcCorrectionTask);

TpcCorrectionTask::TpcCorrectionTask() : 
FairTask("TpcCorrectionTask", 1) 
{}

TpcCorrectionTask::~TpcCorrectionTask() 
{}

InitStatus TpcCorrectionTask::Init()
{
   return kSUCCESS;
}

void TpcCorrectionTask::Exec(Option_t *opt)
{
	Runner R;

	string filepathA = FileHelper::BuildFilePath(Solution::release, Direction::input, "CorrectionMatrix.dat");
	R.LoadCorrectionMatrix(filepathA);
	
	string inputTracks = FileHelper::BuildFilePath(Solution::release, Direction::input, "tracks.inp");
	string outputTracs = FileHelper::BuildFilePath(Solution::release, Direction::output, "tracks.out");
	R.CorrectTracks(inputTracks, outputTracs);
}

void TpcCorrectionTask::Finish() {}