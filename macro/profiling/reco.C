// ROOT includes
#include "TString.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TChain.h"

// Fair includes
#include "FairRunAna.h"
#include "FairFileSource.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTask.h"
#include "FairField.h"

// MPD includes
#include "MpdTpcHitProducer.h"
#include "MpdTpcClusterFinderTask.h"
#include "MpdTpcDigitizerAZ.h"
#include "MpdTpcClusterFinderAZ.h"
#include "MpdTpcClusterFinderMlem.h"
#include "MpdKalmanFilter.h"
#include "MpdVertexZfinder.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdKfPrimaryVertexFinder.h"
#include "MpdTofHitProducer.h"
#include "MpdEtofHitProducer.h"
// #include "MpdEctTrackFinderTpc.h"
// non-existing #include "MpdEctTrackFinderTof.h"
// non-existing #include "MpdEctTrackFinderCpc.h"
#include "MpdTofMatching.h"
#include "MpdZdcDigiProducer.h"
#include "MpdEtofMatching.h"
#include "MpdFillDstTask.h"
#include "MpdGetNumEvents.h"
#include "MpdEmcHitCreation.h"
#include "TpcSectorGeoAZ.h"

#include <iostream>
using namespace std;

bool CheckFileExist(TString fileName)
{
   gSystem->ExpandPathName(fileName);
   if (gSystem->AccessPathName(fileName.Data()) == true) {
      cout << endl << "no specified file: " << fileName << endl;
      return false;
   }

   return true;
}

#define Mlem // Choose: Mlem HitProducer

// Macro for running reconstruction:
// inFile - input file with MC data, default: evetest.root
// nStartEvent - number (start with zero) of first event to process, default: 0
// nEvents - number of events to process, 0 - all events of given file will be proccessed, default: 1
// outFile - output file with reconstructed data, default: mpddst.root
void reco(TString inFile = "$VMCWORKDIR/macro/mpd/evetest.root", TString outFile = "mpddst.root", Int_t nStartEvent = 0,
          Int_t nEvents = 1)
{
   if (!CheckFileExist(inFile)) return;

   // ========================================================================
   // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
   Int_t iVerbose = 0;

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();
   // ------------------------------------------------------------------------

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

   MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
   fRun->AddTask(kalman);

#ifdef Mlem
   MpdTpcDigitizerAZ *tpcDigitizer = new MpdTpcDigitizerAZ(*secGeo);
   tpcDigitizer->SetPersistence(kFALSE);
   fRun->AddTask(tpcDigitizer);
#endif

   //  MpdTpcClusterFinderTask *tpcClusterFinder = new MpdTpcClusterFinderTask();
   //  tpcClusterFinder->SetDebug(kFALSE);
   //  tpcClusterFinder->SetMakeQA(kTRUE);
   //  tpcClusterFinder->SetCalcResiduals(kFALSE);
   //  fRun->AddTask(tpcClusterFinder);

#ifdef Mlem
   MpdTpcClusterFinderMlem *tpcClusAZ = new MpdTpcClusterFinderMlem(*secGeo);
   fRun->AddTask(tpcClusAZ);
#else
   MpdTpcHitProducer *hitPr = new MpdTpcHitProducer(*secGeo);
   hitPr->SetModular(0);
   fRun->AddTask(hitPr);
#endif

   FairTask *vertZ = new MpdVertexZfinder(*secGeo);
   fRun->AddTask(vertZ);

   MpdTpcKalmanFilter *recoKF = new MpdTpcKalmanFilter(*secGeo,"Kalman filter");
#ifdef Mlem
   recoKF->UseTpcHit(kFALSE); // do not use hits from the hit producer
#endif
   fRun->AddTask(recoKF);

   FairTask *findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
   fRun->AddTask(findVtx);

   // TOF hit producers
   MpdTofHitProducer *tofHit = new MpdTofHitProducer("Hit producer");
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

   FairTask *emcHP = new MpdEmcHitCreation();
   fRun->AddTask(emcHP);

   FairTask *tdigi = new MpdZdcDigiProducer("MpdZdcDigiProducer");
   fRun->AddTask(tdigi);

   // MpdPidRefitTrackTask* trRefit = new MpdPidRefitTrackTask("Track PID and Refit");
   // fRun->AddTask(trRefit);

   MpdFillDstTask *fillDST = new MpdFillDstTask("MpdDst task");
   fRun->AddTask(fillDST);

   // -----   Intialise   ----------------------------------------------------
   fRun->Init();
   cout << "Field: " << fRun->GetField()->GetBz(0., 0., 0.) << endl;

   // if nEvents is equal 0 then all events of the given file starting with "nStartEvent" should be processed
   if (nEvents == 0) nEvents = MpdGetNumEvents::GetNumROOTEvents((char *)inFile.Data()) - nStartEvent;

   // -----   Run   ______________--------------------------------------------
   fRun->Run(nStartEvent, nStartEvent + nEvents);
   // ------------------------------------------------------------------------

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
   // ------------------------------------------------------------------------
}

int main(int argc, char **arg)
{
   reco();
}
