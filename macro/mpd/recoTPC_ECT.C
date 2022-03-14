void recoTPC_ECT()
{

  // ========================================================================
  // Verbosity level (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 1;

  // Input file (MC events)
  TString inFile = "mc.root";
  //TString inFile = "mc_10_etof.root";
  //TString inFile = "mc_etof.root";
  //TString inFile = "mu2.root";
  //TString inFile = "urqmd_100_etof.root";
  //TString inFile = "mu2_5.root";

  // Parameter file
  //TString parFile = "testparams.root";
  TString parFile = inFile;

  // ECT hits file
  TString hitFile = "";
  //TString hitFile = "ect.hits_urqmd.root";
  //TString hitFile = "ect_10.hits.root";

  // TPC reco file
  TString tpcFile = "tpc.reco.root";
  //TString tpcFile = "tpc_10.reco.root";
  //TString tpcFile = "urqmd_tpc.reco.root";

  // Output file
  TString outFile = "ect.reco.root";

  // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
  basiclibs();
  gSystem->Load("libGeoBase");
  gSystem->Load("libParBase");
  gSystem->Load("libBase");
  gSystem->Load("libMCStack");
  gSystem->Load("libField");
  gSystem->Load("libGen");
  gSystem->Load("libPassive");
  gSystem->Load("libtpc");
  gSystem->Load("libStrawendcap");
  gSystem->Load("libKalman");
  gSystem->Load("libLHETrack");
  gSystem->Load("libTof");
  gSystem->Load("libEtof");
//gSystem->Load("libZdc");
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
  if (hitFile != "") fRun->AddFriend(hitFile);
  fRun->AddFriend(tpcFile);
  //fRun->AddFriend(inFile);
  fRun->SetOutputFile(outFile);
  // ------------------------------------------------------------------------

  // -----  Parameter database   --------------------------------------------
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(parFile.Data());
  //CbmParAsciiFileIo* parInput2 = new CbmParAsciiFileIo();
  //TString stsDigiFile = gSystem->Getenv("VMCWORKDIR");
  //stsDigiFile += "/parameters/sts/sts_digi_new_standard.par";
  //parInput2->open(stsDigiFile.Data(),"in");
  rtdb->setFirstInput(parInput1);
  //rtdb->setSecondInput(parInput2);
  // fRun->LoadGeometry();  // EL
  // ------------------------------------------------------------------------
  
  MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
  fRun->AddTask(kalman);

  //CbmTask* trackMS = new TpcLheHitsMaker("Tracking routine");
  //fRun->AddTask(trackMS);

  /* unknown class FairTask* reco = new MpdEctTrackFinderTpc();
  fRun->AddTask(reco); */
  
  // Number of events to process
  Int_t nEvents = 100; //1; //50; //250; //90;
  
  // -----   Intialise and run   --------------------------------------------
  fRun->Init();
  cout << "Field: " << fRun->GetField()->GetBz(0.,0.,0.) << endl;
  //CbmRootManager::Instance()->GetInChain()->Print();
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
