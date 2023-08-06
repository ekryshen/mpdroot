bool CheckFileExist(TString fileName){
    gSystem->ExpandPathName(fileName);
    if (gSystem->AccessPathName(fileName.Data()) == true)
    {
        cout<<endl<<"no specified file: "<<fileName<<endl;
        return false;
    }                

    return true;
}


void RunAnalyses(int nEvents = -1, TString inFileList = "list.txt"){

//  gROOT->LoadMacro("mpdloadlibs.C");
//  gROOT->ProcessLine("mpdloadlibs()");

   gSystem->Load("libZdc.so") ;
   gSystem->Load("libEmc.so") ;
   gSystem->Load("libMpdPhotons.so") ;
   gSystem->Load("libMpdPhysics.so") ;

   MpdAnalysisManager man("ManagerAnal", nEvents) ;
   if (!CheckFileExist(inFileList)) return;
   man.InputFileList(inFileList) ;
   man.ReadBranches("*") ; 
   man.SetOutput("histos.root") ;
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;

   MpdEventPlaneAll pEP("pEP","pEP") ;
   man.AddTask(&pEP) ;

   MpdTrackPidMaker pPID("pPID","pPID") ;
   man.AddTask(&pPID) ;
   
//   MpdConvPi0 pDef("pi0Def","ConvDef") ; //name, parametes file
//   man.AddTask(&pDef) ;

//   MpdGlobalPolarizationRECO pGlobalPol("pGlobalPolRECO","pGlobalPolRECO","selection","omega2") ;
//   man.AddTask(&pGlobalPol) ;

   MpdPairKK pKK("pKK","pKK") ;
   man.AddTask(&pKK) ;

   MpdPairPK pPK("pPK","pPK") ;
   man.AddTask(&pPK) ;

   MpdPairPiK pPiK("pPiK","pPiK") ;
   man.AddTask(&pPiK) ;

   MpdPairPiPi pPiPi("pPiPi","pPiPi") ;
   man.AddTask(&pPiPi) ;

   MpdPairPiLambda pPiLambda("pPiLambda","pPiLambda") ;
   man.AddTask(&pPiLambda) ;

   MpdPairPiKs pPiKs("pPiKs","pPiKs") ;
   man.AddTask(&pPiKs) ;

   man.Process() ;

}

