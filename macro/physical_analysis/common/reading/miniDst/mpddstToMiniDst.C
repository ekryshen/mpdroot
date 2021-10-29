#include <Rtypes.h>
#include <TString.h>

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"
        
/*
 A macro to be used when producing a MiniDst from the standard mpddst
 * arg1 - standard mpddst
 * Output name contains arg1_basename.MiniDst.root
*/

void mpddstToMiniDst(TString arg1 = "") {    
    FairRunAna* fRun = new FairRunAna();

    FairSource* fFileSource = new FairFileSource(arg1);
    fRun->SetSource(fFileSource);
    fRun->SetGenerateRunInfo(kFALSE);
    fRun->SetUseFairLinks(kTRUE);
    
    MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
    fRun->AddTask(kalman);
    
    MpdMiniDstFillTask* miniDst = new MpdMiniDstFillTask(arg1);
    fRun->AddTask(miniDst);

    // Intialise ...
    fRun->Init();

    // Run ...
    fRun->Run(0, 0); // (0, 0) means processing of the whole data set given
}