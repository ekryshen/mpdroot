#include "TString.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"

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
#include "MpdGetNumEvents.h" // /23142342234

#include <iostream>

#include "commonFunctions.C"
#include "geometry_stage1.C"

// Available generators
enum EGenerators
{
   BOX = 1,
   FLUID,
   HSD,
   ION,
   LAQGSM,
   MCDST,
   PART,
   SMASH,
   UNIGEN,
   URQMD,
   VHLLE
};

// Available VMC types
enum EVMCType
{
   GEANT3 = 1,
   GEANT4
};

// generator - selected generatory type
// vmc - selected virtual Monte Carlo generator
// inFile - input file with generator data, default: auau.09gev.mbias.98k.ftn14
// nStartEvent - for compatibility, any number
// nEvents - number of events to transport, default: 1
// outFile - output file with MC data, default: evetest.root
// flag_store_FairRadLenPoint
// FieldSwitcher: 0 - corresponds to the ConstantField (0, 0, 5) kG (It is used by default); 1 - corresponds to the
// FieldMap ($VMCWORKDIR/input/B-field_v2.dat)
void runMC(EGenerators generator = EGenerators::BOX, EVMCType vmc = EVMCType::GEANT3, Int_t nStartSeed = 0,
           TString inFile = "auau.04gev.0_3fm.10k.f14.gz", TString outFile = "evetest.root", Int_t nStartEvent = 0,
           Int_t nEvents = 2, Bool_t flag_store_FairRadLenPoint = kFALSE, Int_t FieldSwitcher = 0)
{
   TStopwatch timer;
   timer.Start();
   gDebug = 0;

   FairRunSim *fRun = new FairRunSim();
   // Choose the Geant Navigation System
   switch (vmc)
   {
   case EVMCType::GEANT3:
      fRun->SetName("TGeant3");
      break;
   case EVMCType::GEANT4:
      fRun->SetName("TGeant4");
      break;
   default:
      std::cerr << "Error. Set incorrect VMC type.\n";
      return;
   }

   geometry_stage1(fRun); // load mpd geometry
   // geometry_v2(fRun); // load mpd geometry

   // Use extended MC Event header instead of the default one.
   // MpdMCEventHeader* mcHeader = new MpdMCEventHeader();
   // fRun->SetMCEventHeader(mcHeader);

   // Create and Set Event Generator
   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();
   fRun->SetGenerator(primGen);

   // smearing of beam interaction point
   primGen->SetBeam(0.0, 0.0, 0.1, 0.1);
   primGen->SetTarget(0.0, 24.0);
   primGen->SmearGausVertexZ(kTRUE);
   primGen->SmearVertexXY(kTRUE);

   // Use user defined decays https://fairroot.gsi.de/?q=node/57
   fRun->SetUserDecay(kTRUE);
   // Use external decayer
   // fRun->SetPythiaDecayer(TString("$VMCWORKDIR/gconfig/LambdaDecayConfig.C"));

   switch (generator)
   {
   case EGenerators::BOX: // Box generator
   {
      gRandom->SetSeed(0);
      FairBoxGenerator *boxGen = new FairBoxGenerator(13, 100); // 13 = muon; 1 = multipl.
      boxGen->SetPRange(0.25, 2.5);                             // GeV/c, setPRange vs setPtRange
      boxGen->SetPhiRange(0, 360);                              // Azimuth angle range [degree]
      boxGen->SetThetaRange(0, 180);                            // Polar angle in lab system range [degree]
      boxGen->SetXYZ(0., 0., 0.);                               // mm o cm ??
      primGen->AddGenerator(boxGen);
      break;
   }
   case EGenerators::FLUID:
   {
      if (!CheckFileExist(inFile))
         return;
      Mpd3fdGenerator *fluidGen = new Mpd3fdGenerator(inFile);
      if (nStartEvent > 0)
         fluidGen->SkipEvents(nStartEvent);
      // fluidGen->SetPsiRP(0.); // set fixed Reaction Plane angle [rad] instead of random
      // fluidGen->SetProtonNumberCorrection(79./197.); // Z/A Au for Theseus 2018-03-17-bc2a06d
      primGen->AddGenerator(fluidGen);
      break;
   }
   case EGenerators::HSD: // HSD/PHSD generator
   {
      if (!CheckFileExist(inFile))
         return;
      MpdPHSDGenerator *hsdGen = new MpdPHSDGenerator(inFile.Data());
      // hsdGen->SetPsiRP(0.); // set fixed Reaction Plane angle [rad] instead of random
      // hsdGen->WithHyperonPolarization(); // read extended output with hyperon polarizations
      primGen->AddGenerator(hsdGen);
      if (nStartEvent > 0)
         hsdGen->SkipEvents(nStartEvent);
      // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
      if (nEvents == 0)
         nEvents = MpdGetNumEvents::GetNumPHSDEvents(inFile.Data()) - nStartEvent;
      break;
   }
   case EGenerators::ION: // Ion Generator
   {
      FairIonGenerator *fIongen = new FairIonGenerator(79, 197, 79, 1, 0., 0., 25, 0., 0., -1.);
      primGen->AddGenerator(fIongen);
      break;
   }
   case EGenerators::LAQGSM:
   {
      if (!CheckFileExist(inFile))
         return;
      MpdLAQGSMGenerator *guGen = new MpdLAQGSMGenerator(inFile.Data(), kTRUE, 4, 1 + nStartEvent + nEvents);
      // kTRUE - for NICA/MPD, 1+nStartEvent+nEvents - search ions in selected part of file.
      primGen->AddGenerator(guGen);
      if (nStartEvent > 0)
         guGen->SkipEvents(nStartEvent);
      // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
      if (nEvents == 0)
         nEvents = MpdGetNumEvents::GetNumQGSMEvents(inFile.Data()) - nStartEvent;
      break;
   }
   case EGenerators::MCDST:
   {
      if (!CheckFileExist(inFile))
         return;
      MpdMcDstGenerator *mcDstGen = new MpdMcDstGenerator(inFile);
      primGen->AddGenerator(mcDstGen);
      break;
   }
   case EGenerators::PART: // Particle generator
   {
      FairParticleGenerator *partGen = new FairParticleGenerator(211, 10, 1, 0, 3, 1, 0, 0);
      primGen->AddGenerator(partGen);
      break;
   }
   case EGenerators::SMASH:
   {
      MpdSmashGenerator *smashGen = new MpdSmashGenerator(inFile);
      smashGen->SetRandomRP(kTRUE); // default is also kTRUE
      if (nStartEvent > 0)
         smashGen->SkipEvents(nStartEvent);
      primGen->AddGenerator(smashGen);
      // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
      if (nEvents == 0)
         nEvents = smashGen->GetNeventsInTree() - nStartEvent;
      cout << "newRunMC: MpdSmashGenerator: nEvents = " << nEvents << endl;
      break;
   }
   case EGenerators::UNIGEN:
   {
      if (!CheckFileExist(inFile))
         return;
      Bool_t isSpectatorON = kTRUE; // Does unigen tree have fragments (depends on model)?
      MpdUnigenGenerator *uniGen = new MpdUnigenGenerator(inFile, isSpectatorON);
      // uniGen->SetEventPlane(0., 2.*TMath::Pi());
      primGen->AddGenerator(uniGen);
      break;
   }
   case EGenerators::URQMD:
   {
      if (!CheckFileExist(inFile))
         return;
      MpdUrqmdGenerator *urqmdGen = new MpdUrqmdGenerator(inFile);
      // Event plane angle (in degrees) will be generated by uniform distribution from min to max
      Float_t min = 0.0, max = 30.0;
      urqmdGen->SetEventPlane(min * TMath::DegToRad(), max * TMath::DegToRad());
      primGen->AddGenerator(urqmdGen);
      if (nStartEvent > 0)
         urqmdGen->SkipEvents(nStartEvent);
      // if nEvents is equal 0 then all events (start with nStartEvent) of the given file should be processed
      if (nEvents == 0)
         nEvents = MpdGetNumEvents::GetNumURQMDEvents(inFile.Data()) - nStartEvent;
      break;
   }
   case EGenerators::VHLLE:
   {
      if (!CheckFileExist(inFile))
         return;
      MpdVHLLEGenerator *vhlleGen =
          new MpdVHLLEGenerator(inFile, kTRUE); // kTRUE corresponds to hydro + cascade, kFALSE -- hydro only
      vhlleGen->SkipEvents(0);
      vhlleGen->SetEventPlane(0., 2.*TMath::Pi());
      primGen->AddGenerator(vhlleGen);
      break;
   }
   default:
   {
      std::cerr << "Error. Set incorrect generator type.\n";
      return;
   }
   }
   fRun->SetOutputFile(outFile.Data());

   // Magnetic Field Map - for proper use in the analysis MultiField is necessary here
   MpdMultiField *fField = new MpdMultiField();

   if (FieldSwitcher == 0)
   {
      MpdConstField *fMagField = new MpdConstField();
      fMagField->SetField(0., 0., 5.); // values are in kG:  1T = 10kG
      fMagField->SetFieldRegion(-230, 230, -230, 230, -375, 375);
      fField->AddField(fMagField);
      fRun->SetField(fField);
      cout << "FIELD at (0., 0., 0.) = (" << fMagField->GetBx(0., 0., 0.) << "; " << fMagField->GetBy(0., 0., 0.)
           << "; " << fMagField->GetBz(0., 0., 0.) << ")" << endl;
   }
   else if (FieldSwitcher == 1)
   {
      MpdFieldMap *fMagField = new MpdFieldMap("B-field_v2", "A");
      fMagField->Init();
      fField->AddField(fMagField);
      fRun->SetField(fField);
      cout << "FIELD at (0., 0., 0.) = (" << fMagField->GetBx(0., 0., 0.) << "; " << fMagField->GetBy(0., 0., 0.)
           << "; " << fMagField->GetBz(0., 0., 0.) << ")" << endl;
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

   // AZ - Enable decays of unstable particles
   /*
   MpdStack *stack = NULL;
   if (TString(fRun->GetName()).Contains("TGeant3"))
     stack = (MpdStack*) ((TGeant3*)gMC)->GetStack();
   else stack = (MpdStack*) ((TGeant4*)gMC)->GetStack();
   stack->SetDecayUnstable();
   */

   // -Trajectories Visualization (TGeoManager Only)
   // Set cuts for storing the trajectories
   FairTrajFilter *trajFilter = FairTrajFilter::Instance();
   trajFilter->SetStepSizeCut(0.01); // 1 cm
   //  trajFilter->SetVertexCut(-2000., -2000., 4., 2000., 2000., 100.);
   trajFilter->SetMomentumCutP(.50); // p_lab > 500 MeV
   //  trajFilter->SetEnergyCut(.2, 3.02); // 0 < Etot < 1.04 GeV

   trajFilter->SetStorePrimaries(kTRUE);
   trajFilter->SetStoreSecondaries(kFALSE);

   // Fill the Parameter containers for this run
   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *output = new FairParRootFileIo(kParameterMerged);
   // AZ output->open(parFile.Data());
   output->open(gFile);
   rtdb->setOutput(output);

   MpdMultiFieldPar *Par = (MpdMultiFieldPar *)rtdb->getContainer("MpdMultiFieldPar");
   if (fField)
      Par->SetParameters(fField);
   Par->setInputVersion(fRun->GetRunId(), 1);
   Par->setChanged();
   // Par->printParams();

   rtdb->saveOutput();
   rtdb->print();

   // Transport nEvents
   fRun->Run(nEvents);

   if (generator == EGenerators::LAQGSM)
   {
    TString Pdg_table_name = TString::Format("%s%s%c%s", gSystem->BaseName(inFile.Data()), ".g", (fRun->GetName())[6], ".pdg_table.dat");
    (TDatabasePDG::Instance())->WritePDGTable(Pdg_table_name.Data());
   }
   timer.Stop();
   Double_t rtime = timer.RealTime(), ctime = timer.CpuTime();
   printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);
   cout << "Macro finished successfully." << endl; // marker of successful execution for CDASH
}
