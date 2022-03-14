#if !defined(__CINT__) && !defined(__CLING__)
#include "TString.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"

#include "FairRun.h"
#include "FairRunAna.h"
#include "FairRunSim.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTrajFilter.h"
#include "FairUrqmdGenerator.h"
#include "FairPrimaryGenerator.h"
#include "MpdLAQGSMGenerator.h"
#include "FairCave.h"
#include "FairPipe.h"
#include "FairMagnet.h"

#include "TpcDetector.h"
#include "MpdEmc.h"
#include "MpdTof.h"
#include "MpdZdc.h"
#include "MpdFfd.h"
#include "MpdMultiField.h"
#include "MpdMultiFieldPar.h"
#include "MpdConstField.h"
#include "MpdFieldMap.h"
#include "MpdGetNumEvents.h"

#include <iostream>
#include <vector>
#endif

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"
#include "macro/mpd/geometry_stage1.C"
//#include "macro/mpd/geometry_v2.C"
R__LOAD_LIBRARY(libMpdGenFactory);
#include "simulation/generators/genFactory/MpdGeneratorsFactory.h"

// inFile - input file with generator data, default: auau.09gev.mbias.98k.ftn14
// nStartEvent - for compatibility, any number
// nEvents - number of events to transport, default: 1
// outFile - output file with MC data, default: evetest.root
// flag_store_FairRadLenPoint
// FieldSwitcher: 0 - corresponds to the ConstantField (0, 0, 5) kG (It is used by default); 1 - corresponds to the FieldMap ($VMCWORKDIR/input/B-field_v2.dat)
void runMC_v2(TString inFile = "auau.04gev.0_3fm.10k.f14.gz", TString outFile = "evetest.root", Int_t nStartEvent = 0, Int_t nEvents = 10,
        Bool_t flag_store_FairRadLenPoint = kFALSE, Int_t FieldSwitcher = 0,
        MpdGeneratorType GenType = MpdGeneratorType::HADGEN /*Choose generator: URQMD VHLLE FLUID PART ION BOX HSD LAQGSM HADGEN CUSTOM*/, Bool_t useGeant4 = false)
{
    TStopwatch timer;
    timer.Start();
    gDebug = 0;

    FairRunSim* fRun = new FairRunSim();
    // Choose the Geant Navigation System
    fRun->SetName(useGeant4 ? "TGeant4" : "TGeant3");

    geometry_stage1(fRun); // load mpd geometry
    //geometry_v2(fRun); // load mpd geometry

    // Use extended MC Event header instead of the default one.
    //MpdMCEventHeader* mcHeader = new MpdMCEventHeader();
    //fRun->SetMCEventHeader(mcHeader);

    // Create and Set Event Generator
    FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
    fRun->SetGenerator(primGen);

    // smearing of beam interaction point
    primGen->SetBeam(0.0,0.0,0.1,0.1);
    primGen->SetTarget(0.0,24.0);
    primGen->SmearGausVertexZ(kTRUE);
    primGen->SmearVertexXY(kTRUE);

    // Use user defined decays https://fairroot.gsi.de/?q=node/57
    fRun->SetUserDecay(kTRUE);

    vector<MpdGenerator> Generators;
    //struct MpdGenerator -> {GeneratorType, FilePath, StartEvent, nEvents}
    if (GenType == MpdGeneratorType::CUSTOM)
    {
        Generators.push_back({MpdGeneratorType::PART, "", 0, 0});
        Generators.push_back({MpdGeneratorType::ION, "", 0, 0});
        //and so on... fill vector with all needed generators
    }
    else
        Generators.push_back({GenType, inFile, nStartEvent, nEvents});

    unique_ptr<MpdGeneratorsFactory> genFactory(new MpdGeneratorsFactory());
    vector<shared_ptr<MpdFactoryMadeGenerator>> fmGens;
    
    for (auto Gen : Generators)
    {
        
        switch (Gen.genType)
        {
            case MpdGeneratorType::URQMD:
            case MpdGeneratorType::VHLLE:
            case MpdGeneratorType::FLUID:
            case MpdGeneratorType::HSD:
            case MpdGeneratorType::LAQGSM:
            {
                shared_ptr<MpdFactoryMadeGenerator> fmgen = genFactory->create(Gen);
                if (fmgen.get() == NULL)
                {
                    cout << "unable to create generator" << endl;
                    exit(1);
                }
                fmGens.push_back(fmgen->getptr());
                shared_ptr<FairGenerator> g = fmgen->GetGeneratorPtr();
                primGen->AddGenerator(g.get());
            }
            break;
            case MpdGeneratorType::PART:
            {
                // ------- Particle Generator
                FairParticleGenerator* partGen = new FairParticleGenerator(211, 10, 1, 0, 3, 1, 0, 0);
                primGen->AddGenerator(partGen);
            }
            break;
            case MpdGeneratorType::ION:
            {
                // ------- Ion Generator
                FairIonGenerator *fIongen = new FairIonGenerator(79, 197, 79, 1, 0., 0., 25, 0., 0., -1.);
                primGen->AddGenerator(fIongen);
            }
            break;
            case MpdGeneratorType::BOX:
            {
                gRandom->SetSeed(0);
                // ------- Box Generator
                FairBoxGenerator* boxGen = new FairBoxGenerator(13, 100); // 13 = muon; 1 = multipl.
                boxGen->SetPRange(0.25, 2.5); // GeV/c //setPRange vs setPtRange
                boxGen->SetPhiRange(0, 360); // Azimuth angle range [degree]
                boxGen->SetThetaRange(0, 180); // Polar angle in lab system range [degree]
                boxGen->SetXYZ(0., 0., 0.); // mm o cm ??
                primGen->AddGenerator(boxGen);
            }
            break;
            case MpdGeneratorType::HADGEN:
            {
                THadgen* hadGen = new THadgen();
                hadGen->SetRandomSeed(clock() + time(0));
                hadGen->SetParticleFromPdgCode(0, 196.9665, 79);
                hadGen->SetEnergy(6.5E3);
                MpdGeneralGenerator* generalHad = new MpdGeneralGenerator(hadGen);
                primGen->AddGenerator(generalHad);
            }
            break;
            case MpdGeneratorType::CUSTOM:
            default:
                cout << "unable to find generator type" << endl;
                exit(1);
            break;
        }
    }

    fRun->SetOutputFile(outFile.Data());

    // Magnetic Field Map - for proper use in the analysis MultiField is necessary here
    MpdMultiField* fField = new MpdMultiField();

    if (FieldSwitcher == 0) {
        MpdConstField* fMagField = new MpdConstField();
        fMagField->SetField(0., 0., 5.); // values are in kG:  1T = 10kG
        fMagField->SetFieldRegion(-230, 230, -230, 230, -375, 375);
        fField->AddField(fMagField);
        fRun->SetField(fField);
        cout << "FIELD at (0., 0., 0.) = (" <<
                fMagField->GetBx(0., 0., 0.) << "; " << fMagField->GetBy(0., 0., 0.) << "; " << fMagField->GetBz(0., 0., 0.) << ")" << endl;
    }
    else if (FieldSwitcher == 1) {
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
    //  fRun->AddTask(tpcDigitizer);

    fRun->Init();

    // -Trajectories Visualization (TGeoManager Only)
    // Set cuts for storing the trajectories
    FairTrajFilter* trajFilter = FairTrajFilter::Instance();
    trajFilter->SetStepSizeCut(0.01); // 1 cm
    //  trajFilter->SetVertexCut(-2000., -2000., 4., 2000., 2000., 100.);
    trajFilter->SetMomentumCutP(.50); // p_lab > 500 MeV
    //  trajFilter->SetEnergyCut(.2, 3.02); // 0 < Etot < 1.04 GeV

    trajFilter->SetStorePrimaries(kTRUE);
    trajFilter->SetStoreSecondaries(kFALSE);

    // Fill the Parameter containers for this run
    FairRuntimeDb* rtdb = fRun->GetRuntimeDb();

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
    fRun->Run(nEvents);

    for (auto fmgen : fmGens)
    {
        fmgen->PostActions(fRun);
    }

    timer.Stop();
    Double_t rtime = timer.RealTime(), ctime = timer.CpuTime();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);
    cout << "Macro finished successfully." << endl;     // marker of successful execution for CDASH
}
