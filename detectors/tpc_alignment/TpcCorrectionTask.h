#ifndef TPCCORRECTIONTASK_HH
#define TPCCORRECTIONTASK_HH

#include "FairTask.h"

class TpcCorrectionTask : public FairTask {

public:
   TpcCorrectionTask();
   virtual ~TpcCorrectionTask();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt = "");
   virtual void       Finish();

public:
   ClassDef(TpcCorrectionTask, 1);
};

#endif