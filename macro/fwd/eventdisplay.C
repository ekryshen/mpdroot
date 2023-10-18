#include "TString.h"
#include "MpdEventManager.h"
//#include "FairRunAna.h"

void eventdisplay(TString strSimFile = "mc.root"){
  FairRunAna* fRun = new FairRunAna();
  FairSource* fFileSource = new FairFileSource(strSimFile);

  // set parameter file with MC data and detector geometry
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parIo1 = new FairParRootFileIo();
  parIo1->open(strSimFile); 
  rtdb->setFirstInput(parIo1);

  fRun->SetSource(fFileSource);
  fRun->SetOutputFile("ed.root");

  MpdEventManager* fMan = new MpdEventManager();
  fMan->isOnline = 0;
  
  MpdMCTracks* GeoTrack = new MpdMCTracks("GeoTracks");
  fMan->AddTask(GeoTrack);

  MpdMCPointDraw* FwdPoint = new MpdMCPointDraw("FwdPoint", kOrange, kFullCircle);
  fMan->AddTask(FwdPoint);
  fMan->Init(1,3); // layers
  // fMan->Init(1,4); // single straws
}

