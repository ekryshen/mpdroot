void RunAnalyses(int nEvents = -1){

   gSystem->Load("libEmc.so") ;
   gSystem->Load("libMpdPhysics.so") ;
   gSystem->Load("libMpdPhotons.so") ;

   MpdAnalysisManager man("ManagerAnal", nEvents) ;
   man.InputFileList("list.txt") ;
   man.ReadBranches("*") ; 
   man.SetOutput("histos.root") ;
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;
   
   MpdEventPlaneAll pEP("pEP","pEP") ;
   man.AddTask(&pEP) ;

   man.Process() ;

}
