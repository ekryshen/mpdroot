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
#include "AbstractQA.h"
//#include "TpcLaserGridPreprocess.h"
//#include "TpcLaserGridVelocityMakeMap.h"  //include causes rebuilding dictionary

#include <iostream>

// Choose: UseMlem (MLEM clusterhitfinder)
//         UseHitProducer (simple hit producer without digitizer)
//         comment out the line for Fast clusterhitfinder
#define UseMlem 

#include "commonFunctions.C"

// Macro for running reconstruction:
// inFile - input file with MC data, default: evetest.root
// outFile - output file with reconstructed data, default: mpddst.root
// nStartEvent - number (start with zero) of first event to process, default: 0
// nEvents - number of events to process, 0 - all events of given file will be proccessed, default: 1
// qaSetting - stored in static qaEngineMode enum variable to generate desired QA plots
void runReco_laserMakeDriftVelMap(TString inFile = "evetest_laser_tpc.root", TString outFile = "mpddst_laser_makeMap.root", Int_t nStartEvent = 0,
             Int_t nEvents = 20, EQAMode qaSetting = EQAMode::OFF)
{

   if (!CheckFileExist(inFile)) return;
   // ========================================================================
   // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
   Int_t iVerbose = 0;

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();

   // -----   set QA Engine Mode   -------------------------------------------
   AbstractQA::qaEngineMode = qaSetting; 

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

   TpcLaserGridPreprocess *laserPrep = new TpcLaserGridPreprocess(); // filter muon track ends and (optional) secondaries
   laserPrep->SetMakeQA(kTRUE);
   fRun->AddTask(laserPrep);

   MpdTpcDigitizerAZlt *tpcDigitizer = new MpdTpcDigitizerAZlt(*secGeo);
   tpcDigitizer->SetPersistence(kFALSE);
   fRun->AddTask(tpcDigitizer);
 #ifdef UseMlem
   TpcClusterHitFinderMlem *tpcClus = new TpcClusterHitFinderMlem(*secGeo);
   tpcClus->SetPersistence(kTRUE);
 #else
   TpcClusterHitFinderFast *tpcClus = new TpcClusterHitFinderFast(*secGeo);
 #endif 
   fRun->AddTask(tpcClus);
   
   TpcLaserGridVelocityMakeMap * tpcDriftVelMap = new TpcLaserGridVelocityMakeMap(*secGeo);
   //tpcDriftVelMap->WriteParamsToDB(0, 0);
   tpcDriftVelMap->SetParamFileName("MPDTpcVelocityMap.xml");
   tpcDriftVelMap->SetFilterBackground(kTRUE);
   tpcDriftVelMap->SetMakeQA(kTRUE);
   fRun->AddTask(tpcDriftVelMap);

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
