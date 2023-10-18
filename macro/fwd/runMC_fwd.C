#include "FairRunSim.h"
#include "FairModule.h"
#include "FairCave.h"
#include "FairModule.h"
#include "FairPrimaryGenerator.h"
#include "FairBoxGenerator.h"
#include "MpdMultiField.h"
#include "MpdFieldMap.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "MpdMultiFieldPar.h"
#include "FairRootManager.h"
void runMC_fwd(TString outFileName = "mc.root", Int_t nEvents = 100){
  FairRunSim *fRun = new FairRunSim();
  fRun->SetName("TGeant4");
  fRun->SetMaterials("media.geo");
  
  // geometry
  FairModule *Cave = new FairCave("CAVE");
  Cave->SetGeometryFileName("cave.geo");
  fRun->AddModule(Cave);
  
  FairModule* Fwd = new MpdFwd("FWD");
  fRun->AddModule(Fwd);
  
  // generator
  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
  fRun->SetGenerator(primGen);
  FairBoxGenerator *boxGen = new FairBoxGenerator(13, 10);
  boxGen->SetPtRange(0.1,1.0001);
  boxGen->SetEtaRange(2.0, 2.0001);
  primGen->AddGenerator(boxGen);
  
  // output file
  fRun->SetSink(new FairRootFileSink(outFileName));
  // field
  MpdMultiField *fField = new MpdMultiField();
  
  //  MpdFieldMap *fMagField = new MpdFieldMap("B-field_v2", "A");
  //  fMagField->Init();
  
  MpdConstField *fMagField = new MpdConstField();
  fMagField->SetField(0., 0., 5.); // values are in kG:  1T = 10kG
  fMagField->SetFieldRegion(-230, 230, -230, 230, -375, 375);
  
  fField->AddField(fMagField);
  fRun->SetField(fField);
  
  // store trajectories
  fRun->SetStoreTraj(kTRUE);
   
  fRun->Init();

  // database  
  FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo *output = new FairParRootFileIo(kTRUE);
  output->open(gFile);
  rtdb->setOutput(output);
  MpdMultiFieldPar *Par = (MpdMultiFieldPar *)rtdb->getContainer("MpdMultiFieldPar");
  Par->SetParameters(fField);
  Par->setInputVersion(fRun->GetRunId(), 1);
  Par->setChanged();
  rtdb->saveOutput();
  rtdb->print();
  
  FairRootManager::Instance()->CreateGeometryFile("geom.root");
  
  fRun->Run(nEvents);
}
