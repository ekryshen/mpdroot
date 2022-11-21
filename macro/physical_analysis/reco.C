// Macro for running reconstruction

#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TString.h"
#include "TStopwatch.h"
#include "TSystem.h"

// Fair includes
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTask.h"
#include "FairField.h"
#include "FairTrackParP.h"

// MPD includes
#include "TpcClusterizerTask.h"
#include "TpcDriftTask.h"
#include "TpcHitFinderTask.h"
#include "MpdKalmanFilter.h"
#include "TpcLheHitsMaker.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdTofHitProducer.h"
#include "MpdTofMatching.h"
#include "TpcSectorGeoAZ.h"

#include <iostream>
using namespace std;
#endif

void reco(TString inFile = "evetest.root")
{

  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 1;

  // Input file (MC events)
  //TString inFile = "mc.root";

  // Parameter file
  //TString parFile = "testparams.root";
  TString parFile = inFile;

  // Output file
  //TString outFile = "test.raw.1251-1500.root";
  TString outFile = "mpddst.root";

//   // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  mpdloadlibs(kTRUE);      // load full set of main libraries

  //gSystem->Load("libMpdData");
  gSystem->Load("libXMLIO");

  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/geometry_v2.C");
  geometry_v2(0x0, kFALSE);
  // ------------------------------------------------------------------------

  // ---  Now choose concrete engines for the different tasks   -------------
  // ------------------------------------------------------------------------

  // In general, the following parts need not be touched
  // ========================================================================

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  // ------------------------------------------------------------------------

  // -----   Digitization run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(inFile);
  //fRun->AddFriend(inFile);
  fRun->SetOutputFile(outFile);
  // ------------------------------------------------------------------------

  // -----  Parameter database   --------------------------------------------
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(parFile.Data());
  //FairParAsciiFileIo* parInput2 = new FairParAsciiFileIo();
  //TString stsDigiFile = gSystem->Getenv("VMCWORKDIR");
  //stsDigiFile += "/parameters/sts/sts_digi_new_standard.par";
  //parInput2->open(stsDigiFile.Data(),"in");
  rtdb->setFirstInput(parInput1);
  //rtdb->setSecondInput(parInput2);

  // fRun->LoadGeometry();   // EL
  // ------------------------------------------------------------------------
  
  // -----  Initialize geometry   --------------------------------------------
  BaseTpcSectorGeo *secGeo = new TpcSectorGeoAZ(); 

//   TpcClusterizerTask* tpcClusterizer = new TpcClusterizerTask();
//   fRun->AddTask(tpcClusterizer);
    
//   TpcDriftTask* tpcDrifter = new TpcDriftTask();
//   fRun->AddTask(tpcDrifter);
  
//   TpcHitFinderTask* tphHitFinderTask = new TpcHitFinderTask();
//   fRun->AddTask(tphHitFinderTask);

  MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
  fRun->AddTask(kalman);

  //FairTask* trackMS = new TpcLheHitsMaker("Hit producer");
  //fRun->AddTask(trackMS);

  MpdTpcHitProducer* hitPr = new MpdTpcHitProducer(*secGeo);
  hitPr->SetModular(0);
  fRun->AddTask(hitPr);

  FairTask* vertZ = new MpdVertexZfinder(*secGeo);
  fRun->AddTask(vertZ);

  FairTask* recoKF = new MpdTpcKalmanFilter(*secGeo, "Kalman filter");
  fRun->AddTask(recoKF);

  FairTask* findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
  fRun->AddTask(findVtx);

  MpdTofHitProducer* tofHit = new MpdTofHitProducer("Hit producer");
  fRun->AddTask(tofHit);

  MpdTofMatching* tofMatch = new MpdTofMatching("TOF matching");
  fRun->AddTask(tofMatch);

  FairTask *tdigi= new MpdZdcDigiProducer("MpdZdcDigiProducer");
  fRun->AddTask(tdigi);

  FairTask* fillDST = new  MpdFillDstTask("MpdDst task");
  fRun->AddTask(fillDST);
  
  // Number of events to process
  Int_t nEvents = 1; // 100; //50; //250; //90;
  
  // -----   Intialise and run   --------------------------------------------
  fRun->Init();
  cout << "Field: " << fRun->GetField()->GetBz(0.,0.,0.) << endl;
  fRun->Run(0, nEvents);
  // ------------------------------------------------------------------------

  // -----   Finish   -------------------------------------------------------

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Output file is "    << outFile << endl;
  cout << "Parameter file is " << parFile << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------

  cout << " Test passed" << endl;
  cout << " All ok " << endl;
  exit(0);
}
