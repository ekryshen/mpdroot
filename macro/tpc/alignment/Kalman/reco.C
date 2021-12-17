#include <Rtypes.h>
#define EVENTS_TO_PROCESS 50000
//#define MC_OR_DATA  // emulation or data/mc

#if !defined(__CINT__) && !defined(__CLING__)
// ROOT includes
#include "TChain.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"

// Fair includes
#include "FairField.h"
#include "FairFileSource.h"
#include "FairParRootFileIo.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairTask.h"

// MPD includes
#include "MpdGetNumEvents.h"
#include "MpdTpcAlignmentKalman.h"
#ifdef MC_OR_DATA
#include "MpdTpcClusterFinderMlem.h"
//#include "MpdTpcLaserGridPreprocess.h"
#include "MpdTpcDigitizerAZlt.h"
#else
#include "MpdTpcTrackEmulator.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

void reco(Int_t nEvents, Int_t nStartEvent, TString inFile, TString outFile, TString run_type);

int main()
{
#ifdef MC_OR_DATA
    reco(EVENTS_TO_PROCESS, 0,
         "./data/mc/laser.root",
         "./data/dst_mc_laser_misalign.root", "local");
#else
    reco(EVENTS_TO_PROCESS, 0,
         "./data/mc/emulation.root",
         "./data/dst_mc_laser_misalign.root", "local");
#endif

    return 0;
}

#else // __CINT__

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

#endif

void reco(Int_t nEvents = EVENTS_TO_PROCESS, Int_t nStartEvent = 0, TString inFile = "./data/mc/emulation.root",
          TString outFile = "./data/dst_mc_laser_misalign.root", TString run_type = "local")
{
    // ========================================================================
    // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
    //Int_t iVerbose = 0;

    // -----   Timer   --------------------------------------------------------
    TStopwatch timer;
    timer.Start();

    // -----   Digitization run   -------------------------------------------
    FairRunAna *fRun = new FairRunAna();

    FairSource *fFileSource = new FairFileSource(inFile);
    fRun->SetSource(fFileSource);
    fRun->SetOutputFile(outFile);
    fRun->SetGenerateRunInfo(false);
    fRun->SetUseFairLinks(true);
    // ------------------------------------------------------------------------

    // Parameter file
    TString parFile = inFile;

    // -----  Parameter database   --------------------------------------------
    FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
    FairParRootFileIo *parIo1 = new FairParRootFileIo();
    parIo1->open(parFile.Data());
    rtdb->setFirstInput(parIo1);
    rtdb->setOutput(parIo1);
    rtdb->saveOutput();

    // ------------------------------------------------------------------------

#ifdef MC_OR_DATA
    //auto *laserPrep = new MpdTpcLaserGridPreprocess(); // filter muon track ends and (optional) secondaries
    //laserPrep->SetMakeQA(kFALSE);
    //fRun->AddTask(laserPrep);

    auto *tpcDigitizer = new MpdTpcDigitizerAZlt();
    tpcDigitizer->SetPersistence(kFALSE);
    tpcDigitizer->SetOnlyPrimary(kTRUE);//SetDriftVel(0.0055);
    fRun->AddTask(tpcDigitizer);

    auto *tpcClustering = new MpdTpcClusterFinderMlem();
    tpcClustering->SetPersistence(kTRUE);
    fRun->AddTask(tpcClustering);

    auto *tpcAlign = new MpdTpcAlignmentKalman("laser_mc");
#else
    auto *trackEmu = new MpdTpcTrackEmulator("model");
    trackEmu->stddev_set(.1).mean_free_path_set(1.7).hits_resolution_set(1.0);
    fRun->AddTask(trackEmu);

    auto *tpcAlign = new MpdTpcAlignmentKalman("emulation");
#endif

    //tpcAlign->SetAlignParamsFile("./alignment/MPDTpcAlignment.xml");
    tpcAlign->SetPersistence(kTRUE);
    fRun->AddTask(tpcAlign);

    // -----   Intialise   ----------------------------------------------------
    fRun->Init();

    //if (run_type != "proof")
    //    cout << "Field: " << fRun->GetField()->GetBz(0., 0., 0.) << endl;

    if (!nEvents)
        nEvents = MpdGetNumEvents::GetNumROOTEvents((const char *)inFile.Data()) - nStartEvent;

    // -----   Run   ______________--------------------------------------------
    fRun->Run(nStartEvent, nStartEvent + nEvents);

    // -----   Finish   -------------------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    cout << endl << endl;
    cout << "Macro finished successfully." << endl; // marker of successful execution for CDASH
    cout << "Output file is " << outFile << endl;
    cout << "Parameter file is " << parFile << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
}
