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

void laser(TString outFile = "laser.root", Int_t nStartEvent = 0, Int_t nEvents = 5000)
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
    //  geometry_stage1(fRun); // load mpd geometry
    // geometry_v2(fRun); // load mpd geometry

    //fRun->SetUserConfig("g4Config.C");

    // Use extended MC Event header instead of the default one.
    // MpdMCEventHeader* mcHeader = new MpdMCEventHeader();
    // fRun->SetMCEventHeader(mcHeader);

    // Use user defined decays https://fairroot.gsi.de/?q=node/57
    // fRun->SetUserDecay(kTRUE);
    // fRun->SetUserDecay(const TString& Config);

    // Create and Set Event Generator
    FairPrimaryGenerator *primGen = new FairPrimaryGenerator();
    fRun->SetGenerator(primGen);

    // smearing of beam interaction point
    // primGen->SetBeam(0.0, 0.0, 0.0, 0.0);
    // primGen->SetTarget(0.0, 0.0);
    // primGen->SmearGausVertexZ(kFALSE);
    // primGen->SmearVertexZ(kFALSE);
    // primGen->SmearGausVertexXY(kFALSE);
    // primGen->SmearVertexXY(kFALSE);

    /* gRandom->SetSeed(0);   */

    const Double_t Fieldcage_shift_deg = 0.;//8.;
    const Double_t R = 123 - 3. / TMath::Cos(15.*TMath::DegToRad()); // Membrane_outer_holder_R_edge -
    const Double_t spread_ang = 11.;
    const Double_t pos_offset = 30.;

    auto boxGen = [](const Double_t &ph, const Double_t &th, const Double_t &x, const Double_t &y, const Double_t &z) -> FairBoxGenerator*
    {
        FairBoxGenerator *gen = new FairBoxGenerator(13, 1);
        gen->SetPRange(5.5, 5.5);   // GeV/c, setPRange vs setPtRange
        gen->SetPhiRange(ph, ph);   // Azimuth angle range [degree]
        gen->SetThetaRange(th, th); // Polar angle in lab system range [degree]
        gen->SetXYZ(x, y, z);

        return gen;
    };

    /*
    for (Int_t i = 0; i < 4; ++i)
    {
        Double_t angle = i * 90.;
        Double_t x = R * TMath::Sin(angle * TMath::DegToRad());
        Double_t y = R * TMath::Cos(angle * TMath::DegToRad());

        Double_t phi = 360. - 90. * (i + 1);
        Double_t theta = 90.;

        for (Int_t a = 0; a < 2; ++a)
        {
            Double_t z = a ? -4 * pos_offset : 4 * pos_offset;

            primGen->AddGenerator(boxGen(phi, theta, x, y, z));
        }
    } */

    for (Int_t j = -4; j < 4; ++j)
    {
        for (Int_t i = 0; i < 4; ++i)
        {
            for (Int_t a = -2; a < 2; ++a)
            {
                Double_t angle = i * 90. + Fieldcage_shift_deg;
                Double_t x = R * TMath::Sin(angle * TMath::DegToRad());
                Double_t y = R * TMath::Cos(angle * TMath::DegToRad());
                Double_t z = j < 0 ? j * pos_offset : (j + 1) * pos_offset;

                Double_t phi = 360. - 90. * (i + 1) - Fieldcage_shift_deg + (a < 0 ? (a - 1) * spread_ang : (a + 2) * spread_ang);
                Double_t theta = 90.;

                primGen->AddGenerator(boxGen(phi, theta, x, y, z));

                if (a > -2)
                {
                    phi = 360. - 90. * (i + 1) - Fieldcage_shift_deg + spread_ang* a;
                    primGen->AddGenerator(boxGen(phi, theta, x, y,z));
                }
            }
        }
    }

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

    //  MpdTpcDigitizerTask* tpcDigitizer = new MpdTpcDigitizerTask();
    //  tpcDigitizer->SetOnlyPrimary(kTRUE); /// Digitize only primary track
    //  tpcDigitizer->SetMakeQA(kTRUE);  /// SetMakeQA(kTRUE) prepares Quality
    //  Assurance Histograms tpcDigitizer->SetDiffuse(kFALSE);
    //  tpcDigitizer->SetDebug(kFALSE);
    //  tpcDigitizer->SetDistort(kFALSE);
    //  tpcDigitizer->SetResponse(kFALSE);
    //  tpcDigitizer->SetDistribute(kFALSE);
    //  fRun->AddTask(tpcDigitizer);

    //MpdTpcLaserGridPreprocess *laserPrep = new MpdTpcLaserGridPreprocess(); // filter muon track ends and (optional) secondaries
    //laserPrep->SetMakeQA(kTRUE);
    //fRun->AddTask(laserPrep);

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
