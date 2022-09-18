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

//void reco_with_zdc (TString inFile = "zdctest.root", Int_t number_of_events=0)
//void reco_with_zdc (TString inFile = "/hera/cbm/users/marina/laqgsm/AuAuss11mb/output_1/job_0000/evetest.root", Int_t number_of_events=0)
//void reco_with_zdc (TString inFile = "evetest_mup_158GeV_mod13_1000ev_oldgeom84mods_oldProcHits.root", Int_t number_of_events=0)
//void reco_with_zdc (TString inFile = "evetest_mup_158GeV_mod24_1000ev_newgeom_oldProcHits.root", Int_t number_of_events=0)
//void reco_with_zdc (TString inFile = "evetest_p_40GeV_mod25_1000ev_newProcHits_FscScint.root", Int_t number_of_events=0)
//void reco_with_zdc (TString inFile = "evetest_2ev_onlyPSD.root", Int_t number_of_events=0)
void reco_with_zdc (TString inFile = "evetest_10ev_new.root", Int_t number_of_events=0)

//void reco_with_zdc (TString inFile = "evetest_40ev_node207.root", Int_t number_of_events=0)
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
  /*
  if (number_of_events)
    TString outFile = Form("mpddst_%devents_%s",number_of_events,inFile.Data());
  else
    TString outFile = Form("mpddst_%s",inFile.Data());
  */
    TString outFile = Form("mpddst.root");

//   // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  //  mpdloadlibs(kTRUE);      // load full set of main libraries
  mpdloadlibs(1,1);      // load full set of libraries

  //gSystem->Load("libMpdData");
  gSystem->Load("libXMLIO");

  // gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/geometry_v1.C");
  // geometry_v1(0x0, kFALSE);
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

   // -----  Initialize geometry   --------------------------------------------
   BaseTpcGeo *secGeo = new TpcSectorGeoAZ(); 

  // fRun->LoadGeometry();   // EL
  // ------------------------------------------------------------------------
  
//   TpcClusterizerTask* tpcClusterizer = new TpcClusterizerTask();
//   fRun->AddTask(tpcClusterizer);
    
//   TpcDriftTask* tpcDrifter = new TpcDriftTask();
//   fRun->AddTask(tpcDrifter);
  
//   TpcHitFinderTask* tphHitFinderTask = new TpcHitFinderTask();
//   fRun->AddTask(tphHitFinderTask);

  MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
  //fRun->AddTask(kalman);

  MpdTpcHitProducer* hitPr = new MpdTpcHitProducer(*secGeo);
  //hitPr->SetModular(0);
  //fRun->AddTask(hitPr);

  // FairTask* trackMS = new TpcLheHitsMaker("Hit producer");
  // fRun->AddTask(trackMS);

  FairTask* vertZ = new MpdVertexZfinder(*secGeo);
  //fRun->AddTask(vertZ);

  FairTask* recoKF = new MpdTpcKalmanFilter(*secGeo,"Kalman filter");
  //fRun->AddTask(recoKF);

  FairTask* findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
  //fRun->AddTask(findVtx);

  MpdTofHitProducer* tofHit = new MpdTofHitProducer("Hit producer");
  //fRun->AddTask(tofHit);

  MpdTofMatching* tofMatch = new MpdTofMatching("TOF matching");
  //fRun->AddTask(tofMatch);

  FairTask *tdigi= new MpdZdcDigiProducer("MpdZdcDigiProducer");
  fRun->AddTask(tdigi);

  FairTask* fillDST = new  MpdFillDstTask("MpdDst task");
  fRun->AddTask(fillDST);
  
  // Number of events to process
  Int_t nEvents = 1; // 100; //50; //250; //90;
  if (number_of_events)
    nEvents = number_of_events;
  
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
