#ifndef TPC_ALIGNMENT_TASK_LASER_RAYS_HH
#define TPC_ALIGNMENT_TASK_LASER_RAYS_HH

#include "FairTask.h"

class TpcAlignmentTaskLaserRays : public FairTask {
private:
public:
   TpcAlignmentTaskLaserRays();
   virtual ~TpcAlignmentTaskLaserRays();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt = "");
   virtual void       Finish();

public:
   ClassDef(TpcAlignmentTaskLaserRays, 1);
};
#endif
