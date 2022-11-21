void runRecoCpc(TString inFile = "mc.1k.root")
{

  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 1;

  // Input file (MC events)
  //TString inFile = "mc.root";
  //TString inFile = "mc_10.root";
  TString mcFile = TString(gSystem->Getenv("MCFILE"));
  if (mcFile != "") inFile = mcFile;

  // Parameter file
  //TString parFile = "testparams.root";
  TString parFile = inFile;

  // Output file
  //TString outFile = "test.raw.1251-1500.root";
  TString outFile = "tpc.reco.1k.root";
  TString recoFile = TString(gSystem->Getenv("RECOFILE"));
  if (recoFile != "") outFile = recoFile;

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
 // gROOT->ProcessLine(".x $HOME/mpd/loadlibs.C");
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
  // ------------------------------------------------------------------------
  

  // -----  Initialize geometry   --------------------------------------------
  BaseTpcSectorGeo *secGeo = new TpcSectorGeoAZ(); 

  MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
  fRun->AddTask(kalman);
  //kalman->SetNumer(0);

  MpdTpcDigitizerAZ* tpcDigitizer = new MpdTpcDigitizerAZ(*secGeo);
  tpcDigitizer->SetPersistence(kFALSE);
  //fRun->AddTask(tpcDigitizer);

  /*
  MpdTpcClusterFinderTask *tpcClusterFinder = new MpdTpcClusterFinderTask();
  //  tpcClusterFinder->SetDebug(kFALSE);
  //  tpcClusterFinder->SetMakeQA(kTRUE);
  //  tpcClusterFinder->SetCalcResiduals(kFALSE);
  fRun->AddTask(tpcClusterFinder);
  */

  MpdTpcClusterFinderAZ *tpcClusAZ = new MpdTpcClusterFinderAZ(*secGeo);
  //MpdTpcClusterFinderMlem *tpcClusAZ = new MpdTpcClusterFinderMlem(*secGeo);
  tpcClusAZ->SetPersistence();
  //fRun->AddTask(tpcClusAZ);

  MpdTpcHitProducer* hitPr = new MpdTpcHitProducer(*secGeo);
  //hitPr->SetModular(1);
  //hitPr->SetModular(0);
  fRun->AddTask(hitPr);

  FairTask* vertZ = new MpdVertexZfinder(*secGeo);
  fRun->AddTask(vertZ);

  MpdTpcKalmanFilter* recoKF = new MpdTpcKalmanFilter(*secGeo, "Kalman filter");
  //recoKF->UseTpcHit(kFALSE); // do not use hits from the hit producer
  fRun->AddTask(recoKF);
  
  MpdKfPrimaryVertexFinder* findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
  //findVtx->SetConstrFlag();
  fRun->AddTask(findVtx);

  // TOF hit producers
  MpdTofHitProducer* tofHit = new MpdTofHitProducer("TOF Hit producer");
  fRun->AddTask(tofHit);

  MpdEtofHitProducer* etofHitProd = new MpdEtofHitProducer("ETOF HitProducer");
  fRun->AddTask(etofHitProd);

  // Endcap tracking
  /* undefined class
  FairTask* tpcECT = new MpdEctTrackFinderTpc();
  fRun->AddTask(tpcECT);
  */

  //MpdEctTrackFinderTof* tofECT = new MpdEctTrackFinderTof();
  /* undefined class
  MpdEctTrackFinderCpc* tofECT = new MpdEctTrackFinderCpc();
  tofECT->SetTpc(kTRUE);
  fRun->AddTask(tofECT);
  */

  // TOF matching
  MpdTofMatching* tofMatch = new MpdTofMatching("TOF matching");
  fRun->AddTask(tofMatch);

  // ETOF matching
  MpdEtofMatching* etofMatch = new MpdEtofMatching("ETOF matching");
  fRun->AddTask(etofMatch);
  

  FairTask* fillDST = new  MpdFillDstTask("MpdDst task");
  fRun->AddTask(fillDST);

  // Number of events to process
  Int_t nEvents = 1000; //55; //50; //250; //90;
  
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

}
