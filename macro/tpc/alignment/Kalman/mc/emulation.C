#if !defined(__CINT__) && !defined(__CLING__)
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"

#include "FairCave.h"
#include "FairMagnet.h"
#include "FairParRootFileIo.h"
#include "FairPipe.h"
#include "FairPrimaryGenerator.h"
#include "FairRunSim.h"
#include "FairRuntimeDb.h"
#include "FairTrajFilter.h"
#include "FairBoxGenerator.h"

#include "MpdConstField.h"
#include "MpdFieldMap.h"
//#include "MpdGetNumEvents.h"
#include "MpdMultiField.h"
#include "MpdMultiFieldPar.h"
#include "TpcDetector.h"
#include "MpdTpcLaserGridPreprocess.h"

#include <iostream>
using namespace std;
#endif

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/tpc/tpc_only.C"
#include "macro/mpd/mpdloadlibs.C"

#define GEANT3 // Choose: GEANT3 GEANT4

void emulation(TString outFile = "emulation.root", Int_t nStartEvent = 0, Int_t nEvents = 9999990)
{
    TStopwatch timer;
    timer.Start();
    gDebug = 0;

    FairRunSim *fRun = new FairRunSim();
    // Choose the Geant Navigation System
#ifdef GEANT3
    fRun->SetName("TGeant3");
#else
    fRun->SetName("TGeant4");
#endif

    tpc_only(fRun);

    fRun->SetOutputFile(outFile.Data());

    // Magnetic Field Map - for proper use in the analysis MultiField is necessary here
    MpdMultiField *fField = new MpdMultiField();
    MpdConstField* fMagField = new MpdConstField();

    fMagField->SetField(0., 0., 0.);
    fMagField->SetFieldRegion(-230, 230, -230, 230, -375, 375);
    fField->AddField(fMagField);
    fRun->SetField(fField);
    //cout << "FIELD at (0., 0., 0.) = (" << fMagField->GetBx(0., 0., 0.) << "; " << fMagField->GetBy(0., 0.,  0.) << "; " << fMagField->GetBz(0., 0., 0.) << ")" << endl;

    fRun->SetStoreTraj(kTRUE);


    fRun->Init();

    // Set cuts for storing the trajectories
    FairTrajFilter *trajFilter = FairTrajFilter::Instance();
    trajFilter->SetStepSizeCut(0.01); // 1 cm
    //  trajFilter->SetVertexCut(-2000., -2000., 4., 2000., 2000., 100.);
    trajFilter->SetMomentumCutP(.50); // p_lab > 500 MeV
    //  trajFilter->SetEnergyCut(.2, 3.02); // 0 < Etot < 1.04 GeV

    trajFilter->SetStorePrimaries(kTRUE);
    trajFilter->SetStoreSecondaries(kFALSE);

    FairParRootFileIo *output = new FairParRootFileIo(kTRUE);
    output->open(gFile);

    // Fill the Parameter containers for this run
    FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
    rtdb->setOutput(output);

    MpdMultiFieldPar *Par = (MpdMultiFieldPar*) rtdb->getContainer("MpdMultiFieldPar");
    if (fField)
        Par->SetParameters(fField);
    Par->setInputVersion(fRun->GetRunId(), 1);
    Par->setChanged();
    // Par->printParams();

    rtdb->saveOutput();
    rtdb->print();

    fRun->Run(nEvents);

    timer.Stop();
    Double_t rtime = timer.RealTime(), ctime = timer.CpuTime();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);
    //cout << "Macro finished successfully." << endl; // marker of successful execution for CDASH
}
