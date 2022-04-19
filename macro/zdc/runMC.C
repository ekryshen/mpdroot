// Macro for running Fair  with Geant3  or Geant4 (M. Al-Turany , D. Bertini)
// Modified 22/06/2005 D.Bertini

// Last modified 18/07/2013 P.Batyuk (MPD)

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"

#include "FairRunSim.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTrajFilter.h"
#include "FairUrqmdGenerator.h"
#include "FairPrimaryGenerator.h"
#include "FairCave.h"
#include "FairPipe.h"
#include "FairMagnet.h"

#include "MpdSts.h"
#include "TpcDetector.h"
#include "MpdEtof.h"
#include "MpdFsa.h"
#include "MpdBbc.h"
#include "MpdCpc.h"
#include "MpdTof.h"
#include "MpdStrawendcap.h"
#include "MpdZdc.h"
#include "MpdFfd.h"
#include "MpdGetNumEvents.h"

#include <iostream>
using namespace std;
#endif

// inFile - input file with generator data, default: auau.09gev.mbias.98k.ftn14
// nStartEvent - for compatibility, any number
// nEvents - number of events to transport, default: 1
// outFile - output file with MC data, default: evetest.root
// flag_store_FairRadLenPoint
// FieldSwitcher: 0 - corresponds to the ConstantField (0, 0, 5) kG (It is used by default); 1 - corresponds to the FieldMap ($VMCWORKDIR/input/B-field_v2.dat)

void runMC(TString inFile = "auau.04gev.0_3fm.10k.f14", TString outFile = "$VMCWORKDIR/macro/mpd/test/evetest.root", Int_t nStartEvent = 0, Int_t nEvents = 1, Int_t nSkip = 0,Bool_t flag_store_FairRadLenPoint = kFALSE, Int_t FieldSwitcher = 0) {

#define LAQGSM
//#define URQMD
    TStopwatch timer;
    timer.Start();
    gDebug = 0;

    gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
    mpdloadlibs(1, 1); // load main libraries

    //gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/geometry_v2.C");
    //gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/geometry_stage1.C");
    gROOT->LoadMacro("/lustre/nyx/cbm/users/marina/laqgsm/AuAuss11mb/geometry_stage1.C");

    FairRunSim *fRun = new FairRunSim();

    // Choose the Geant Navigation System
    //fRun->SetName("TGeant3");
     fRun->SetName("TGeant4");
    // fRun->SetGeoModel("G3Native");

    geometry_stage1(fRun, kTRUE); // load mpd geometry

    // Use the experiment specific MC Event header instead of the default one                                                     
    // This one stores additional information about the reaction plane                                                           
    //MpdMCEventHeader* mcHeader = new MpdMCEventHeader();
    FairMCEventHeader* mcHeader = new FairMCEventHeader();
                                                                         
    fRun->SetMCEventHeader(mcHeader);

    // Create and Set Event Generator
    //-------------------------------

    FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
    fRun->SetGenerator(primGen);

    // smearing of beam interaction point
    //primGen->SetBeam(0.0,0.0,0.1,0.1);

    primGen->SetTarget(0.0,24.0);
    primGen->SmearGausVertexZ(kTRUE);

    //primGen->SmearVertexXY(kTRUE);

#ifdef URQMD
    // ------- Urqmd  Generator
    TString hostname = gSystem->HostName(), dataFile;

    if (inFile.Contains("/"))
        dataFile = inFile;
    else {
        dataFile = find_path_to_URQMD_files();
        if ((hostname == "lxmpd-ui.jinr.ru") || (hostname == "lxmpd-ui"))
            dataFile += "auau.09gev.mbias.10k.f14";
        else
            dataFile += inFile;
    }

    if (!CheckFileExist(dataFile)) return;

    MpdUrqmdGenerator* urqmdGen = new MpdUrqmdGenerator(dataFile);
    primGen->AddGenerator(urqmdGen);
    if (nStartEvent > 0) urqmdGen->SkipEvents(nStartEvent);

    // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
    if (nEvents == 0)
        nEvents = MpdGetNumEvents::GetNumURQMDEvents(dataFile.Data()) - nStartEvent;

#else
#ifdef FLUID

    Mpd3fdGenerator* fluidGen = new Mpd3fdGenerator(inFile);
    fluidGen->SkipEvents(0);
    primGen->AddGenerator(fluidGen);
    //if (nStartEvent > 0) urqmdGen->SkipEvents(nStartEvent);

    // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
    //    if (nEvents == 0)
    //        nEvents = MpdGetNumEvents::GetNumURQMDEvents(dataFile.Data()) - nStartEvent;

#else

#ifdef PART
    // ------- Particle Generator
    FairParticleGenerator* partGen =
            new FairParticleGenerator(211, 10, 1, 0, 3, 1, 0, 0);
    primGen->AddGenerator(partGen);

#else
#ifdef ION
    // ------- Ion Generator
    FairIonGenerator *fIongen =
            new FairIonGenerator(79, 197, 79, 1, 0., 0., 25, 0., 0., -1.);
    primGen->AddGenerator(fIongen);

#else
#ifdef BOX
    gRandom->SetSeed(0);
    // ------- Box Generator
    FairBoxGenerator* boxGen = new
            FairBoxGenerator(13, 1); // 13 = muon; 1 = multipl.
    boxGen->SetPRange(0.25, 2.5); // GeV/c //setPRange vs setPtRange
    boxGen->SetPhiRange(0, 360); // Azimuth angle range [degree]
    boxGen->SetThetaRange(0, 180); // Polar angle in lab system range [degree]
    boxGen->SetXYZ(0., 0., 0.); // mm o cm ??
    primGen->AddGenerator(boxGen);

#else
#ifdef HSD
    // ------- HSD/PHSD Generator
    TString dataFile;
~    if (inFile.Contains("/"))
        dataFile = inFile;
    else {
        dataFile = find_path_to_URQMD_files();
        dataFile += "/../../HSD/"; //  nc-farm
        dataFile += inFile;
    }

    if (!CheckFileExist(dataFile)) return;

    MpdPHSDGenerator *hsdGen = new MpdPHSDGenerator(dataFile.Data());
    //hsdGen->SetPsiRP(0.); // set fixed Reaction Plane angle instead of random
    primGen->AddGenerator(hsdGen);
    if (nStartEvent > 0) hsdGen->SkipEvents(nStartEvent);

    // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
    if (nEvents == 0)
        nEvents = MpdGetNumEvents::GetNumPHSDEvents(dataFile.Data()) - nStartEvent;

#else
#ifdef LAQGSM
    // ------- LAQGSM Generator
/*
    TString dataFile;
    if (inFile.Contains("/"))
        dataFile = inFile;
    else {
        dataFile = find_path_to_URQMD_files();
        dataFile += "/../../QGSM/"; //  nc-farm
        dataFile += inFile;
    }
*/

    TString dataFile;
    dataFile = inFile;
    if (!CheckFileExist(dataFile)) return;

    //nEvents = 2;
    //Int_t nSkip=0;
    Bool_t use_collider_system = kTRUE;
    Int_t QGSM_format_ID;

    //inFile = "/hera/cbm/users/marina/laqgsm/AuAuss11mb/AuAuss11mb_1.r12";
    //inFile = "AuAuss11mb_1.r12";
    dataFile = inFile;

    //MpdLAQGSMGenerator* guGen = new MpdLAQGSMGenerator(inFile);
    MpdLAQGSMGenerator* guGen = new MpdLAQGSMGenerator(dataFile.Data(),use_collider_system,QGSM_format_ID);
    primGen->AddGenerator(guGen);
    //if (nStartEvent > 0) guGen->SkipEvents(nStartEvent);
	cout <<"runMC " <<nSkip <<endl;
    guGen->SkipEvents(nSkip);
    //primGen->SetTarget(50.0,24.0);

    // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
    //if (nEvents == 0)
    //nEvents = MpdGetNumEvents::GetNumQGSMEvents(dataFile.Data()) - nStartEvent;

/*
    TString dataFile;
    dataFile = inFile;
    if (!CheckFileExist(dataFile)) return;

    MpdLAQGSMGenerator* guGen = new MpdLAQGSMGenerator(dataFile.Data());
    primGen->AddGenerator(guGen);
    if (nStartEvent > 0) guGen->SkipEvents(nStartEvent);

    // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
    if (nEvents == 0)
        nEvents = MpdGetNumEvents::GetNumQGSMEvents(dataFile.Data()) - nStartEvent;
*/
#endif
#endif
#endif
#endif
#endif
#endif
#endif

    fRun->SetOutputFile(outFile.Data());

    // Magnetic Field Map - for proper use in the analysis MultiField is necessary here
    // --------------------
    MpdMultiField *fField = new MpdMultiField();

    if (FieldSwitcher == 0) {
        MpdConstField *fMagField = new MpdConstField();
        fMagField->SetField(0., 0., 5.); // values are in kG:  1T = 10kG
        fMagField->SetFieldRegion(-230, 230, -230, 230, -375, 375);
        fField->AddField(fMagField);
        fRun->SetField(fField);
        cout << "FIELD at (0., 0., 0.) = (" <<
                fMagField->GetBx(0., 0., 0.) << "; " << fMagField->GetBy(0., 0., 0.) << "; " << fMagField->GetBz(0., 0., 0.) << ")" << endl;
    } else if (FieldSwitcher == 1) {
        MpdFieldMap* fMagField = new MpdFieldMap("B-field_v2", "A");
        fMagField->Init();
        fField->AddField(fMagField);
        fRun->SetField(fField);
        cout << "FIELD at (0., 0., 0.) = (" <<
                fMagField->GetBx(0., 0., 0.) << "; " << fMagField->GetBy(0., 0., 0.) << "; " << fMagField->GetBz(0., 0., 0.) << ")" << endl;
    }

    fRun->SetStoreTraj(kTRUE);
    fRun->SetRadLenRegister(flag_store_FairRadLenPoint); // radiation length manager

    //  MpdTpcDigitizerTask* tpcDigitizer = new MpdTpcDigitizerTask();
    //  tpcDigitizer->SetOnlyPrimary(kTRUE); /// Digitize only primary track
    //  tpcDigitizer->SetMakeQA(kTRUE);  /// SetMakeQA(kTRUE) prepares Quality Assurance Histograms  
    //  tpcDigitizer->SetDiffuse(kFALSE);
    //  tpcDigitizer->SetDebug(kFALSE);
    //  tpcDigitizer->SetDistort(kFALSE);
    //  tpcDigitizer->SetResponse(kFALSE);
    //  tpcDigitizer->SetDistribute(kFALSE);
    //fRun->AddTask(tpcDigitizer);

    fRun->Init();

    // -Trajectories Visualization (TGeoManager Only )
    // -----------------------------------------------

    // Set cuts for storing the trajectories
    FairTrajFilter* trajFilter = FairTrajFilter::Instance();
    trajFilter->SetStepSizeCut(0.01); // 1 cm
    //  trajFilter->SetVertexCut(-2000., -2000., 4., 2000., 2000., 100.);
    trajFilter->SetMomentumCutP(.50); // p_lab > 500 MeV
    //  trajFilter->SetEnergyCut(.2, 3.02); // 0 < Etot < 1.04 GeV

    //trajFilter->SetStorePrimaries(kTRUE);
    trajFilter->SetStorePrimaries(kFALSE);
    trajFilter->SetStoreSecondaries(kFALSE);

    // Fill the Parameter containers for this run
    //-------------------------------------------

    FairRuntimeDb *rtdb = fRun->GetRuntimeDb();

    Bool_t kParameterMerged = kTRUE;
    FairParRootFileIo* output = new FairParRootFileIo(kParameterMerged);
    //AZ output->open(parFile.Data());
    output->open(gFile);
    rtdb->setOutput(output);

    MpdMultiFieldPar* Par = (MpdMultiFieldPar*) rtdb->getContainer("MpdMultiFieldPar");
    if (fField)
        Par->SetParameters(fField);
    Par->setInputVersion(fRun->GetRunId(), 1);
    Par->setChanged();
    // Par->printParams();

    rtdb->saveOutput();
    rtdb->print();

    // Transport nEvents
    // -----------------
    fRun->Run(nEvents);

#ifdef LAQGSM
    TString Pdg_table_name = TString::Format("%s%s%c%s", gSystem->BaseName(dataFile.Data()), ".g", (fRun->GetName())[6], ".pdg_table.dat");
    (TDatabasePDG::Instance())->WritePDGTable(Pdg_table_name.Data());
#endif

    Bool_t file = fRun->GetWriteRunInfoFile();
    timer.Stop();
    Double_t rtime = timer.RealTime(), ctime = timer.CpuTime();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);

    cout << "Macro finished succesfully." << endl;

    gApplication->Terminate();
}
