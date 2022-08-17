#ifndef TPCALIGNMENTTASK_HH
#define TPCALIGNMENTTASK_HH

#include "FairTask.h"

class TpcAlignmentTask : public FairTask {
private:
  bool vDebugMode{false};
  int vNumberOfCalibration{6};

public:
  TpcAlignmentTask();
  TpcAlignmentTask(bool debugMode);

  virtual ~TpcAlignmentTask();

  virtual InitStatus Init();
  virtual void Exec(Option_t *opt = "");
  virtual void Finish();
  void SetNumberOfCalibration(int val=6);

public:
  ClassDef(TpcAlignmentTask, 1);
};

#endif