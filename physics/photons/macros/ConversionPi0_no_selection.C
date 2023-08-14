void ConversionPi0_no_selection(int nEvents = -1)
{

   gSystem->Load("libZdc.so");
   gSystem->Load("libMpdPhysics.so");
   gSystem->Load("libEmc.so");
   gSystem->Load("libMpdPhotons.so");

   MpdAnalysisManager man("ManagerPi0", nEvents);
   man.InputFileList("list.txt");
   man.ReadBranches("*");
   man.SetOutput("histos.root");

   MpdCentralityAll pCentr("pCentr", "pCentr");
   man.AddTask(&pCentr);

   MpdEventPlaneAll pEP("pEP", "pEP");
   man.AddTask(&pEP);

   MpdTrackPidMaker pPID("pPID", "pPID");
   man.AddTask(&pPID);

   MpdV0Maker pV0maker("pV0maker", "pV0maker");
   man.AddTask(&pV0maker);

   MpdConvPi0 pi0_no_selection("pi0_no_selection", "pi0_no_selection");
   man.AddTask(&pi0_no_selection);

   man.Process();
}
