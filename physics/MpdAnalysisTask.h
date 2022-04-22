#ifndef MPDANALYSISTASK_H
#define MPDANALYSISTASK_H

#include "MpdEvent.h"
#include "MpdAnalysisEvent.h"
#include "TList.h"

class MpdAnalysisTask {

public:

  MpdAnalysisTask() = default;
  MpdAnalysisTask(const char *name, const char *outputName="taskName");
  virtual ~MpdAnalysisTask() ;  // Destructor

  //Methods to be re-defined in real class
  virtual void  UserInit(){}      //Method called before scan e.g. to create histograms and trees to fill 
  virtual void  ProcessEvent(MpdAnalysisEvent &event){}  //Process one event
  virtual void  Finish() {}       // Method called after scan and before writing output

  //List of output object to be written to output file
  TList * GetOutput(){ return fOutputList;}
  const char * GetName(){return fTaskName.Data();}
  const char * GetOutputName(){return fOutputName.Data();}
protected:
  TString fTaskName ;
  TString fOutputName ;
  TList * fOutputList = nullptr;   ///< List of output objects to be written to output file   

  ClassDef(MpdAnalysisTask,0);
};

#endif
