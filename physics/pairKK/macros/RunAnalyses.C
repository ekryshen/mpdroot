void RunAnalyses(){

  gROOT->LoadMacro("mpdloadlibs.C");
  gROOT->ProcessLine("mpdloadlibs()");

   MpdAnalysisManager man("ManagerAnal") ;
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
