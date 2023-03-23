void RunAnalyses(int nEvents = -1){

  //gROOT->LoadMacro("mpdloadlibs.C");
  //gROOT->ProcessLine("mpdloadlibs()");

   gSystem->Load("libEmc.so") ;
   gSystem->Load("libMpdPhysics.so") ;
   gSystem->Load("libMpdPhotons.so") ;

   MpdAnalysisManager man("ManagerAnal", nEvents) ;
   man.InputFileList("list.txt") ;
   man.ReadBranches("*") ; 
   man.SetOutput("histos.root") ;
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;
   
//   MpdConvPi0 pDef("pi0Def","ConvDef") ; //name, parametes file
//   man.AddTask(&pDef) ;

   MpdPairKK pKK("pKK","pKK") ;
   man.AddTask(&pKK) ;

   man.Process() ;

}
