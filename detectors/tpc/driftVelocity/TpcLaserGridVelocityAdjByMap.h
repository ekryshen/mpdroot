#ifndef TPCLASERGRIDVELOCITYADJBYMAP_H
#define TPCLASERGRIDVELOCITYADJBYMAP_H

#include <FairTask.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TString.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "BaseTpcSectorGeo.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <thread>
#include <mutex>

#include "TpcLaserGridConsts.h"

class TpcLaserGridVelocityAdjByMap : public FairTask, TpcLaserGridConsts {
public:
   TpcLaserGridVelocityAdjByMap(BaseTpcSectorGeo &sectorGeo);
   virtual ~TpcLaserGridVelocityAdjByMap();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   virtual void       Finish();

   void ReadParamsFromDB(size_t pID, size_t rID, bool rfDB = kTRUE)
   {
      periodID   = pID;
      runID      = rID;
      readFromDB = rfDB;
   }
   void SetParamFileName(std::string path) { mapFile = path; }
   void UseZY(bool opt = kTRUE);
   // void SetMakeQA(bool opt = kTRUE) { makeQA = opt; }

private:
   TString           inputBranchName;
   TClonesArray     *hits;
   BaseTpcSectorGeo *secGeo;
   std::string       mapFile;

   double driftDistance;
   double driftTheorVelocity;
   int    sectorsNum;

   // std::vector<double> y_layers;//(loc_y_offset.size() + 2);

   int                                           events_num;
   std::vector<std::vector<std::vector<double>>> map_vel_harm;
   std::vector<std::vector<std::vector<double>>> map_corr_coeff;

   std::vector<std::vector<gsl_spline *>> corr_coeff_z;
   std::vector<std::vector<gsl_spline *>> vel_z;

   std::vector<gsl_interp_accel *> accel;

   int                      threadPoolSize = 24;
   std::vector<std::thread> threadPool;
   std::mutex               getTpcHit_lock;
   std::mutex               sectGeo_lock;
   void                     ThreadProcessing(int hitIdxLow, int hitIdxHigh, gsl_interp_accel *accel);

   // bool makeQA = kFALSE;

   bool   readFromDB = kFALSE;
   size_t periodID   = 0;
   size_t runID      = 0;
   int    ReadMap(bool useDB);

public:
   ClassDef(TpcLaserGridVelocityAdjByMap, 0)
};

#endif
