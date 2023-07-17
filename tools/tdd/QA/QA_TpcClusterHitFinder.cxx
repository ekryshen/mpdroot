#include "QA_TpcClusterHitFinder.h"
#include "FairLogger.h"

ClassImp(QA_TpcClusterHitFinder);

void QA_TpcClusterHitFinder::WriteToFile()
{
   WriteBaseQA(moduleNameSuffix);
   WriteTpcClusterHitFinderQA();
}

void QA_TpcClusterHitFinder::WriteTpcClusterHitFinderQA()
{
   TString outputFilename = TString("QA_TpcClusterHitFinder_") + moduleNameSuffix + TString(".root");
   LOG(info) << " Writing QA file: " << outputFilename;
}
