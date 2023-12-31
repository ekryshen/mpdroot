void runKF(TString inFile = "../mpd/evetest.root")
{

  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 1;


  // Parameter file
  //TString parFile = "testparams.root";
  TString parFile = inFile;

  // Output file
  //TString outFile = "test.raw.1251-1500.root";
  TString outFile = "tpc.reco.root";

//   // ----  Load libraries   -------------------------------------------------
//   gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
//   basiclibs();
//   gSystem->Load("libGeoBase");
//   gSystem->Load("libParBase");
//   gSystem->Load("libBase");
//   gSystem->Load("libMCStack");
//   gSystem->Load("libField");
//   gSystem->Load("libGen");
//   gSystem->Load("libPassive");

//   gSystem->Load("libKalman");
//   gSystem->Load("libLHETrack");

  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  mpdloadlibs(kTRUE);                 // load full set of main libraries

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
  // fRun->LoadGeometry();  // EL
  // ------------------------------------------------------------------------
  
  MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
  fRun->AddTask(kalman); 

  //FairTask* trackMS = new TpcLheHitsMaker("Hit producer");
  //fRun->AddTask(trackMS);

  MpdTpcHitProducer* hitPr = new MpdTpcHitProducer();
  hitPr->SetModular(0);
  fRun->AddTask(hitPr);

  FairTask* vertZ = new MpdVertexZfinder();
  fRun->AddTask(vertZ);

  FairTask* recoKF = new MpdTpcKalmanFilter("Kalman filter");
  fRun->AddTask(recoKF);
  
  FairTask* findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
  fRun->AddTask(findVtx);
  
  
       // -----   TOF  producers   ---------------------------------------------
  MpdTofHitProducer* tofHitProd = new MpdTofHitProducer("TOF HitProducer", iVerbose, true);
  tofHitProd->SetParamFlnm("tof.geo.par.xml");
  fRun->AddTask(tofHitProd);
  
  MpdTofMatching* tofMatching = new MpdTofMatching("TOF MatchingProducer", iVerbose, true);
  tofMatching->SetParamFlnm("tof.geo.par.xml");
  fRun->AddTask(tofMatching);
   // -----   TOF  producers   ---------------------------------------------  
  
  
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
