void ConversionPi0(){

   gSystem->Load("libEmc.so") ;
   gSystem->Load("libMpdPhysics.so") ;
   gSystem->Load("libMpdPhotons.so") ;


   MpdAnalysisManager man("ManagerPi0") ;
   man.InputFileList("list.txt") ;
   man.ReadBranches("*") ; 
   man.SetOutput("histos.root") ;
   
   
   MpdConvPi0 pDef("pi0Def","ConvDef") ; //name, parametes file
   man.AddTask(&pDef) ;

   MpdConvPi0 pMinE("pi0MinE","ConvEmin") ;
   man.AddTask(&pMinE) ;

   man.Process() ;

}
