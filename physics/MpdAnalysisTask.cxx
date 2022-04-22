#include "MpdAnalysisTask.h"

ClassImp(MpdAnalysisTask);


MpdAnalysisTask::MpdAnalysisTask(const char *name, const char *outputName):
fTaskName(name),fOutputName(outputName)
{
  fOutputList=nullptr;
}
MpdAnalysisTask::~MpdAnalysisTask(){
  if(fOutputList){
  	delete fOutputList;
  	fOutputList=nullptr;
  }
}