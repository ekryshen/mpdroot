// Macro for running physic analysis

#if !defined(__CINT__) && !defined(__CLING__)
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
#include "MpdRoInvMassTask.h"

#include <iostream>
using namespace std;
#endif

void phys_RoInvMass(TString inFile = "mpddst.root")
{

  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 3;

  // Input file (MC events)
  //TString inFile = "mc.root";

  // Parameter file
  TString parFile = "$VMCWORKDIR/macro/mpd/evetest.root";
  //TString parFile = inFile;

  // Output file
  //TString outFile = "test.raw.1251-1500.root";
  TString outFile = "mpdphys.root";

  //gSystem->Load("libXMLIO");

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
  
  FairTask* physics = new  MpdRoInvMassTask("MpdPhysics task");
  fRun->AddTask(physics);
  
  // Number of events to process
  Int_t nEvents = 1; // 100; //50; //250; //90;
  
  // -----   Intialise and run   --------------------------------------------
  fRun->Init();
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
