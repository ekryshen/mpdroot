// Macro for running reconstruction

#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TString.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TChain.h"

// Fair includes
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTask.h"
#include "FairField.h"
#include "FairTrackParP.h"

// MPD includes
#include "MpdTpcHitProducer.h"
#include "TpcClearerTask.h"
#include "TpcClusterizerTask.h"
//#include "TpcClusterizerTaskQA.h"
#include "TpcDriftTask.h"
//#include "TpcDriftTaskQA.h"
#include "TpcMWPCTask.h"
//#include "TpcMWPCTaskQA.h"
#include "TpcPadResponseTask.h"
//#include "TpcPadResponseTaskQA.h"
#include "TpcADCTask.h"
//#include "MpdTpcClusterFinderTaskQA.h"
#include "TpcDistributor.h"
#include "TpcHitFinderTask.h"
#include "MpdKalmanFilter.h"
#include "TpcLheHitsMaker.h"
#include "MpdVertexZfinder.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdKfPrimaryVertexFinder.h"
#include "MpdTofHitProducer.h"
#include "MpdEtofHitProducer.h"
// #include "MpdEctTrackFinderTpc.h"
// undefined #include "MpdEctTrackFinderTof.h"
#include "MpdTofMatching.h"
#include "MpdEtofMatching.h"
#include "MpdFillDstTask.h"
#include "MpdGetNumEvents.h"
#include "TpcSectorGeoAZ.h"

#include <iostream>
using namespace std;
#endif

// inFile - input file with MC data, default: evetest.root
// nStartEvent - number (start with zero) of first event to process, default: 0
// nEvents - number of events to process, 0 - all events of given file will be proccessed, default: 1
// outFile - output file with reconstructed data, default: mpddst.root
void recoTPC(TString inFile = "$VMCWORKDIR/macro/mpd/evetest.root", TString outFile = "mpddst.root", Int_t nStartEvent = 0, Int_t nEvents = 1)
{
  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 1;

  // Parameter file
  TString parFile = inFile;

  // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  mpdloadlibs(kTRUE);       // load full set of main libraries

  gSystem->Load("libXMLIO");

  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/geometry_v2.C");
  geometry_v2(0x0, kFALSE);
  // ------------------------------------------------------------------------

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
  rtdb->setFirstInput(parInput1);

  // fRun->LoadGeometry();  // EL

  // ------------------------------------------------------------------------

  // -----  Initialize geometry   --------------------------------------------
  BaseTpcSectorGeo *secGeo = new TpcSectorGeoAZ(); 

  MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
  fRun->AddTask(kalman);

  //FairTask* trackMS = new TpcLheHitsMaker("Hit producer");
  //fRun->AddTask(trackMS);

//  TpcDistributor* tpcDistributor = new TpcDistributor(10000, kTRUE, kTRUE);
//  fRun->AddTask(tpcDistributor);

  MpdTpcHitProducer* hitPr = new MpdTpcHitProducer(*secGeo);
  hitPr->SetModular(0);
  fRun->AddTask(hitPr);

  FairTask* vertZ = new MpdVertexZfinder(*secGeo);
  fRun->AddTask(vertZ);

  FairTask* recoKF = new MpdTpcKalmanFilter(*secGeo, "Kalman filter");
  fRun->AddTask(recoKF);

  FairTask* findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
  fRun->AddTask(findVtx);

  // TOF hit producers
  // MpdTofHitProducer* tofHit = new MpdTofHitProducer("Hit producer");
  // fRun->AddTask(tofHit);

  // TOF matching
  // MpdTofMatching* tofMatch = new MpdTofMatching("TOF matching");
  // fRun->AddTask(tofMatch);

  FairTask *tdigi= new MpdZdcDigiProducer("MpdZdcDigiProducer");
  fRun->AddTask(tdigi);

  MpdFillDstTask* fillDST = new MpdFillDstTask("MpdDst task");
  fRun->AddTask(fillDST);

  // -----   Intialise   ----------------------------------------------------
  fRun->Init();
  cout << "Field: " << fRun->GetField()->GetBz(0.,0.,0.) << endl;

  // if nEvents is equal 0 then all events of the given file starting with "nStartEvent" should be processed
  if (nEvents == 0)
      nEvents = MpdGetNumEvents::GetNumROOTEvents(inFile.Data()) - nStartEvent;

  // -----   Run   ______________--------------------------------------------
  fRun->Run(nStartEvent, nStartEvent+nEvents);
  // ------------------------------------------------------------------------

  // -----   Finish   -------------------------------------------------------

  delete fRun;

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Output file is "    << outFile << endl;
  cout << "Parameter file is " << parFile << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << "Macro finished succesfully." << endl;
  cout << endl;
  // ------------------------------------------------------------------------
}
