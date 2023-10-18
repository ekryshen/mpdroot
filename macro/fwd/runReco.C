#if !defined(__CINT__) && !defined(__CLING__)
#include "FairRunAna.h"
#include "FairFileSource.h"
#include "FairRootFileSink.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "MpdFwdHitProducer.h"
#include "MpdFwdTrackProducer.h"
#endif

void runReco(TString inFileName = "mc.root", TString outFileName = "reco.root"){
  FairRunAna* fRun = new FairRunAna();
  fRun->SetSource(new FairFileSource(inFileName));
  fRun->SetSink(new FairRootFileSink(outFileName));
  
  FairRuntimeDb     *rtdb   = fRun->GetRuntimeDb();
  FairParRootFileIo *parIo1 = new FairParRootFileIo();
  parIo1->open(inFileName);
  rtdb->setFirstInput(parIo1);
  rtdb->setOutput(parIo1);
  rtdb->saveOutput();

  FairTask* fwdHitProducer = new MpdFwdHitProducer();
  fRun->AddTask(fwdHitProducer);

  FairTask* fwdTrackProducer = new MpdFwdTrackProducer();
  fRun->AddTask(fwdTrackProducer);

  fRun->Init();
  fRun->Run();
}
