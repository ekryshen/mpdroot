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
#include "TpcLheKalmanFilter.h"
#include "MpdTofHitProducer.h"
#include "MpdTofMatching.h"

#include <iostream>
using namespace std;
#endif

void reco(TString inFile="auau_09gev_0_3fm_0.root",TString  outFile="dst_auau_09gev_0_3fm_0.root",Int_t nEvents = 200)
{

  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 1;

  // Number of events to process
  //  Int_t nEvents = 200; // 100; //50; //250; //90;
  
  // Input file (MC events)
  // Parameter file
  TString parFile = inFile;

  // Output file
  //TString outFile = "test.raw.1251-1500.root";
  //  TString outDir = "./dst";
  //  TString outFile = outDir+"/"+inFile;
  //  TString parDir = "./par";
  //  TString  fPDFtestFlnm = parDir+"/"+inFile;

//   // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  mpdloadlibs(kTRUE);      // load full set of main libraries

  //gSystem->Load("libMpdData");
  gSystem->Load("libXMLIO");

  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/geometry_v1.C");
  geometry_v1(0x0, kFALSE);
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

  fRun->LoadGeometry();
  // ------------------------------------------------------------------------
  
//   TpcClusterizerTask* tpcClusterizer = new TpcClusterizerTask();
//   fRun->AddTask(tpcClusterizer);
    
//   TpcDriftTask* tpcDrifter = new TpcDriftTask();
//   fRun->AddTask(tpcDrifter);
  
//   TpcHitFinderTask* tphHitFinderTask = new TpcHitFinderTask();
//   fRun->AddTask(tphHitFinderTask);

  MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
  fRun->AddTask(kalman);

  FairTask* trackMS = new TpcLheHitsMaker("Hit producer");
  fRun->AddTask(trackMS);

  FairTask* recoKF = new MpdTpcKalmanFilter("Kalman filter");
  fRun->AddTask(recoKF);

  MpdTofHitProducer* tofHit = new MpdTofHitProducer("Hit producer");
  fRun->AddTask(tofHit);

  MpdTofMatching* tofMatch = new MpdTofMatching("TOF matching"); // ,1, "kTRUE");
  fRun->AddTask(tofMatch);

  FairTask* fillDST = new  MpdFillDstTask("MpdDst task");
  fRun->AddTask(fillDST);
  
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
