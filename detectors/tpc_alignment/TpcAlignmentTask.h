#ifndef TPCALIGNMENTTASK_HH
#define TPCALIGNMENTTASK_HH

#include "FairTask.h"

class TpcAlignmentTask : public FairTask {
private:

public:
   TpcAlignmentTask();
   virtual ~TpcAlignmentTask();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt = "");
   virtual void       Finish();

public:
   ClassDef(TpcAlignmentTask, 1);
};

#endif