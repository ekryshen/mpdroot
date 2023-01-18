#include <Rtypes.h>
// ROOT includes
#include "TString.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TChain.h"

// Fair includes
#include "FairFileSource.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTask.h"
#include "FairField.h"
#include "FairRunAna.h"

// MPD includes
#include "MpdTpcHitProducer.h"
#include "MpdTpcDigitizerAZlt.h"
#include "MpdTpcClusterFinderAZ.h"
#include "MpdTpcClusterFinderMlem.h"
#include "MpdKalmanFilter.h"
#include "MpdVertexZfinder.h"
#include "MpdTpcTracker.h"
#include "MpdKfPrimaryVertexFinder.h"
#include "MpdFfdHitProducer.h"
#include "MpdTofHitProducer.h"
#include "MpdEtofHitProducer.h"
#include "MpdTofMatching.h"
#include "MpdZdcDigiProducer.h"
#include "MpdEtofMatching.h"
#include "MpdFillDstTask.h"
#include "MpdGetNumEvents.h"
#include "MpdEmcHitCreation.h"
#include "MpdPid.h"
#include "TpcSectorGeoAZ.h"
#include "TpcClusterHitFinderFast.h"

#include <iostream>

// Choose: UseMlem (MLEM clusterhitfinder)
//         UseHitProducer (simple hit producer without digitizer)
//         comment out the line for Fast clusterhitfinder
//#define UseMlem 

#include "commonFunctions.C"

// Macro for running reconstruction:
// inFile - input file with MC data, default: evetest.root
// nStartEvent - number (start with zero) of first event to process, default: 0
// nEvents - number of events to process, 0 - all events of given file will be proccessed, default: 1
// outFile - output file with reconstructed data, default: mpddst.root
void trackingActs(TString inFile = "evetest.root", TString outFile = "mpddst.root", Int_t nStartEvent = 0,
                  Int_t nEvents = 10)
{

   if (!CheckFileExist(inFile)) return;
   // ========================================================================
   // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
   Int_t iVerbose = 0;

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();

   // -----   Digitization run   -------------------------------------------
   FairRunAna *fRun;
   fRun = new FairRunAna();

   FairSource *fFileSource = new FairFileSource(inFile);
   fRun->SetSource(fFileSource);
   fRun->SetOutputFile(outFile);
   fRun->SetGenerateRunInfo(false);
   fRun->SetUseFairLinks(true);
   // ------------------------------------------------------------------------

   // Parameter file
   TString parFile = inFile;

   // -----  Parameter database   --------------------------------------------
   FairRuntimeDb     *rtdb   = fRun->GetRuntimeDb();
   FairParRootFileIo *parIo1 = new FairParRootFileIo();
   parIo1->open(parFile.Data());
   rtdb->setFirstInput(parIo1);
   rtdb->setOutput(parIo1);
   rtdb->saveOutput();
   // ------------------------------------------------------------------------

   // -----  Initialize geometry   --------------------------------------------
   BaseTpcSectorGeo *secGeo = new TpcSectorGeoAZ();

   // ------------------------------------------------------------------------

   MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
   fRun->AddTask(kalman);

#ifdef UseHitProducer
   MpdTpcHitProducer *hitPr = new MpdTpcHitProducer(*secGeo);
   hitPr->SetModular(0);
   fRun->AddTask(hitPr);
#else
   // MpdTpcDigitizerAZ* tpcDigitizer = new MpdTpcDigitizerAZ(*secGeo);
   MpdTpcDigitizerAZlt *tpcDigitizer = new MpdTpcDigitizerAZlt(*secGeo);
   tpcDigitizer->SetPersistence(kTRUE);
   fRun->AddTask(tpcDigitizer);
 #ifdef UseMlem
   MpdTpcClusterFinderMlem *tpcClusAZ = new MpdTpcClusterFinderMlem(*secGeo);
 #else
   TpcClusterHitFinderFast *tpcClusAZ = new TpcClusterHitFinderFast(*secGeo);
 #endif 
   fRun->AddTask(tpcClusAZ);
#endif

   FairTask *vertZ = new MpdVertexZfinder(*secGeo);
   fRun->AddTask(vertZ);

   MpdTpcTracker *recoKF = new MpdTpcTracker(*secGeo);
#ifndef UseHitProducer
//   recoKF->UseTpcHit(kFALSE); // do not use hits from the hit producer
#endif
   fRun->AddTask(recoKF);

//FIXME:   FairTask *findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
//FIXME:   fRun->AddTask(findVtx);

//FIXME:   MpdFfdHitProducer *ffdHit = new MpdFfdHitProducer("FFDHitProducer");
//FIXME:   fRun->AddTask(ffdHit);

   // TOF hit producers
//FIXME:   MpdTofHitProducer *tofHit = new MpdTofHitProducer("Hit producer");
//FIXME:   tofHit->SetTimeResolution(0.080);
//FIXME:   fRun->AddTask(tofHit);

   /*
       MpdEtofHitProducer* etofHitProd = new MpdEtofHitProducer("ETOF HitProducer");
       fRun->AddTask(etofHitProd);

       // Endcap tracking
       FairTask* tpcECT = new MpdEctTrackFinderTpc(); // undefined class
       tpcECT->SetVerbose(iVerbose);
       fRun->AddTask(tpcECT);

       MpdEctTrackFinderCpc* tofECT = new MpdEctTrackFinderCpc(); // undefined class
       tofECT->SetVerbose(iVerbose);
       tofECT->SetTpc(kTRUE);
       fRun->AddTask(tofECT);
   */

   // TOF matching
//FIXME:   MpdTofMatching *tofMatch = new MpdTofMatching("TOF matching");
//FIXME:   fRun->AddTask(tofMatch);

   // ETOF matching
   // MpdEtofMatching* etofMatch = new MpdEtofMatching("ETOF matching");
   // fRun->AddTask(etofMatch);

   // FairTask *emcHP = new MpdEmcHitCreation();
   // fRun->AddTask(emcHP);

//FIXME:   FairTask *tdigi = new MpdZdcDigiProducer("MpdZdcDigiProducer");
//FIXME:   fRun->AddTask(tdigi);

   // MpdPidRefitTrackTask* trRefit = new MpdPidRefitTrackTask("Track PID and Refit");
   // fRun->AddTask(trRefit);

   MpdFillDstTask *fillDST = new MpdFillDstTask("MpdDst task");
   //fRun->AddTask(fillDST);

   MpdMiniDstFillTask *miniDst = new MpdMiniDstFillTask(outFile);
   //fRun->AddTask(miniDst);

   // -----   Intialise   ----------------------------------------------------
   fRun->Init();
   cout << "Field: " << fRun->GetField()->GetBz(0., 0., 0.) << endl;

   // if nEvents is equal 0 then all events of the given file starting with "nStartEvent" should be processed
   if (nEvents == 0) nEvents = MpdGetNumEvents::GetNumROOTEvents((char *)inFile.Data()) - nStartEvent;

   // -----   Run   ______________--------------------------------------------
   fRun->Run(nStartEvent, nStartEvent + nEvents);

   // -----   Finish   -------------------------------------------------------
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl << endl;
   cout << "Macro finished successfully." << endl; // marker of successful execution for CDASH
   cout << "Output file is " << outFile << endl;
   cout << "Parameter file is " << parFile << endl;
   cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
   cout << endl;
}
