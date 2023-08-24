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
#include "MpdKalmanFilter.h"
#include "MpdVertexZfinder.h"
#include "MpdTpcKalmanFilter.h"
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
#include "TpcClusterHitFinderMlem.h"
#include "QA_TpcClusterHitFinder.h"

#include <iostream>

#include "commonFunctions.C"

// available clustering / hit finder modules for TPC
enum ETpcClustering {
   HITPRODUCER, // simple hit producer without digitizer
   MLEM,        // default Mlem module (AZ)
   FAST
};

// Macro for running reconstruction:
// inFile - input file with MC data, default: evetest.root
// outFile - output file with reconstructed data, default: mpddst.root
// nStartEvent - number (start with zero) of first event to process, default: 0
// nEvents - number of events to process, 0 - all events of given file will be proccessed, default: 1
// qaSetting - enum variable defining which QA data are collected
void runReco(TString inFile = "evetest.root", TString outFile = "mpddst.root", Int_t nStartEvent = 0,
             Int_t nEvents = 10, ETpcClustering tpcClusteringModule = ETpcClustering::MLEM,
             EQAMode qaSetting = EQAMode::OFF)
{

   if (!CheckFileExist(inFile)) return;
   // ========================================================================
   // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
   Int_t iVerbose = 0;

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();

   // -----   set QA Engine Mode   -------------------------------------------
   BaseQA *qaObject;
   switch (qaSetting) {
   case EQAMode::OFF: {
      qaObject = nullptr;
      break;
   }
   case EQAMode::BASIC: {
      qaObject = new BaseQA();
      break;
   }
   case EQAMode::TPCCLUSTERHITFINDER: {
      qaObject = new QA_TpcClusterHitFinder();
      break;
   }
   default: {
      std::cerr << "Error. You've set the non-existing QA Engine mode.\n";
      return;
   }
   }

   // -----   Digitization run   ---------------------------------------------
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

   if (tpcClusteringModule == ETpcClustering::HITPRODUCER) {
      MpdTpcHitProducer *hitPr = new MpdTpcHitProducer(*secGeo);
      hitPr->SetModular(0);
      fRun->AddTask(hitPr);
   } else {
      MpdTpcDigitizerAZlt *tpcDigitizer = new MpdTpcDigitizerAZlt(*secGeo);
      tpcDigitizer->SetPersistence(kTRUE);
      fRun->AddTask(tpcDigitizer);

      AbstractTpcClusterHitFinder *tpcClus;
      switch (tpcClusteringModule) {
      case ETpcClustering::MLEM: {
         tpcClus = new TpcClusterHitFinderMlem(*secGeo, qaObject);
         break;
      }
      case ETpcClustering::FAST: {
         tpcClus = new TpcClusterHitFinderFast(*secGeo, qaObject);
         break;
      }
      default: {
         std::cerr << "Error. You've set the non-existing Tpc Clustering Module.\n";
         return;
      }
      }
      fRun->AddTask(tpcClus);
   }

   FairTask *vertZ = new MpdVertexZfinder(*secGeo);
   fRun->AddTask(vertZ);

   MpdTpcKalmanFilter *recoKF = new MpdTpcKalmanFilter(*secGeo, "Kalman filter");
   if (tpcClusteringModule != ETpcClustering::HITPRODUCER)
      recoKF->UseTpcHit(kFALSE); // do not use hits from the Tpc hit producer

   fRun->AddTask(recoKF);

   FairTask *findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
   fRun->AddTask(findVtx);

   MpdFfdHitProducer *ffdHit = new MpdFfdHitProducer("FFDHitProducer");
   fRun->AddTask(ffdHit);

   // TOF hit producers
   MpdTofHitProducer *tofHit = new MpdTofHitProducer("Hit producer");
   tofHit->SetTimeResolution(0.080);
   fRun->AddTask(tofHit);

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
   MpdTofMatching *tofMatch = new MpdTofMatching("TOF matching");
   fRun->AddTask(tofMatch);

   // ETOF matching
   // MpdEtofMatching* etofMatch = new MpdEtofMatching("ETOF matching");
   // fRun->AddTask(etofMatch);

   // FairTask *emcHP = new MpdEmcHitCreation();
   // fRun->AddTask(emcHP);

   FairTask *tdigi = new MpdZdcDigiProducer("MpdZdcDigiProducer");
   fRun->AddTask(tdigi);

   // MpdPidRefitTrackTask* trRefit = new MpdPidRefitTrackTask("Track PID and Refit");
   // fRun->AddTask(trRefit);

   MpdFillDstTask *fillDST = new MpdFillDstTask(qaObject, "MpdDst task");
   fRun->AddTask(fillDST);

   MpdMiniDstFillTask *miniDst = new MpdMiniDstFillTask(outFile);
   fRun->AddTask(miniDst);

   // -----   Initialize   ----------------------------------------------------
   fRun->Init();
   cout << "Field: " << fRun->GetField()->GetBz(0., 0., 0.) << endl;

   // if nEvents is equal 0 then all events of the given file starting with "nStartEvent" should be processed
   if (nEvents == 0) nEvents = MpdGetNumEvents::GetNumROOTEvents((char *)inFile.Data()) - nStartEvent;

   // -----   Run   ______________--------------------------------------------
   fRun->Run(nStartEvent, nStartEvent + nEvents);

   // -----   QA Engine Output  -------------------------------------------
   if (qaObject) qaObject->WriteToFile();

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
