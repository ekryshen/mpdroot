#ifndef TPC_CORRECTION_TASK_LASER_RAYS_HH
#define TPC_CORRECTION_TASK_LASER_RAYS_HH

#include "FairTask.h"

class TpcCorrectionTaskLaserRays : public FairTask {

public:
   TpcCorrectionTaskLaserRays();
   virtual ~TpcCorrectionTaskLaserRays();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt = "");
   virtual void       Finish();

public:
   ClassDef(TpcCorrectionTaskLaserRays, 1);
};
#endif
