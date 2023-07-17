#include "BaseQA.h"
#include "FairLogger.h"
ClassImp(BaseQA);

void BaseQA::WriteBaseQA(TString suffix)
{
   TString outputFile("BaseQA.root");
   if (suffix.Length() > 0) outputFile = TString("BaseQA_") + suffix + TString(".root");
   LOG(info) << " Writing QA file: " << outputFile;
}