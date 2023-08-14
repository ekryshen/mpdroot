#ifndef TPCLASERGRIDVELOCITYMAKEMAP_H
#define TPCLASERGRIDVELOCITYMAKEMAP_H

#include "FairTask.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

#include "BaseTpcSectorGeo.h"

#include <cmath>
#include <memory>
#include <vector>
#include <thread>
#include <mutex>

#include "TpcLaserGridConsts.h"

static double fit_func_gaussian(double x, double a, double b, double c)
{
   const double z = (x - b) / c;
   return a * std::exp(-0.5 * z * z);
}

class TpcLaserGridVelocityMakeMap : public FairTask, TpcLaserGridConsts {
public:
   struct TpcLaserGridDistribution {
      double              minval;
      double              maxval;
      std::vector<double> distr;

      TpcLaserGridDistribution() : distr(0), minval(0), maxval(0) {}
      TpcLaserGridDistribution(double min, double max, int nbins)
         : // distr(nbins),
           minval(min), maxval(max)
      {
         distr.resize(nbins);
      }
      virtual ~TpcLaserGridDistribution() {}

      virtual int FindBin(double val)
      {
         assert(maxval > minval);
         assert(distr.size() > 0);

         const double dval = DeltaVal();
         int          idx  = (val - minval) / dval;
         if (idx < 0 || idx >= distr.size()) return -1;
         return idx;
      }

      virtual bool Add(double val)
      {
         assert(maxval > minval);
         assert(distr.size() > 0);

         int idx = FindBin(val);
         if (idx < 0) return false;
         distr[idx] += 1.;
         return true;
      }

      virtual double DeltaVal() const { return fabs(maxval - minval) / distr.size(); }

      virtual double GetParam(int ibin) { return minval + (ibin + 0.5) * DeltaVal(); }

      virtual std::pair<std::vector<double>, std::vector<double>> GetData(int lowbin, int highbin)
      {
         assert(maxval > minval);
         assert(highbin > lowbin);
         assert(distr.size() > lowbin);
         assert(distr.size() > highbin);
         assert(distr.size() > 0);

         int size = highbin - lowbin;

         std::vector<double> vals(size), params(size);

         const double dval  = DeltaVal();
         double       param = GetParam(lowbin);
         for (auto &p : params) {
            p = param;
            param += dval;
         }

         std::copy(distr.begin() + lowbin, distr.begin() + highbin, vals.begin());

         return std::make_pair(params, vals);
      }

      std::vector<int> findPeaks(int m = 30, double paramLowThreshold = 45.)
      {
         std::vector<double> shape(distr.size());
         for (int i = 1; i < distr.size() - 1; i++) {
            shape[i] = (distr[i] - distr[i - 1]) * (distr[i + 1] - distr[i]);
         }

         std::vector<int> pks;
         for (int i = 1; i < shape.size() - 1; i++) {
            if (shape[i] < 0) {
               int z        = i - m + 1;
               z            = (z > 0) ? z : 1;
               int w        = i + m + 1;
               w            = (w < distr.size()) ? w : distr.size() - 1;
               bool is_peak = true;
               for (int j = z; j < w; j++) {
                  if (distr[j] > distr[i]) {
                     is_peak = false;
                     break;
                  }
               }
               if (is_peak) {
                  pks.push_back(i + 1);
               }
            }
         }

         std::vector<int> filteredPeaks;
         for (int i = 0; i < pks.size(); ++i) {
            // std::cout << "Peak pos raw " << pks[i] << std::endl;
            if (pks[i] > 0 && pks[i] < distr.size() - 1)
               if (distr[pks[i] - 1] > paramLowThreshold || distr[pks[i]] > paramLowThreshold ||
                   distr[pks[i] + 1] > paramLowThreshold) {
                  filteredPeaks.push_back(pks[i]);
                  // std::cout << "Peak pos " << pks[i] << std::endl;
               }
         }

         return filteredPeaks;
      }

      std::vector<int> findPeaksSimple(double paramLowThreshold = 45.)
      {
         std::vector<int> pks;
         bool             inPeak = false;
         for (int i = 0; i < distr.size(); ++i) {
            if (!inPeak && distr[i] > paramLowThreshold) {
               pks.push_back(i);
               inPeak = true;
               // std::cout << "Peak pos " << i << std::endl;
            }
            if (distr[i] <= paramLowThreshold) inPeak = false;
         }
         return pks;
      }

      // ClassDef(TpcLaserGridDistribution, 0)
   };

public:
   TpcLaserGridVelocityMakeMap(BaseTpcSectorGeo &sectorGeo);
   virtual ~TpcLaserGridVelocityMakeMap();

   virtual InitStatus Init() override;
   virtual void       Exec(Option_t *opt) override;
   virtual void       Finish() override;

   void SetParamFileName(std::string path) { mapFile = path; }
   void WriteParamsToDB(size_t pID, size_t rID, bool wtDB = kTRUE)
   {
      periodID  = pID;
      runID     = rID;
      writeToDB = wtDB;
   }
   void SetFilterBackground(bool opt = kTRUE, double factor = 1.0)
   {
      useFilter     = opt;
      lin_fit_bmult = factor;
   }
   void UseZY(bool opt = kTRUE);
   void SetMakeQA(bool opt = kTRUE) { makeQA = opt; }

private:
   void SectorProcessing(const int secID, std::vector<TpcLaserGridDistribution> *hitZDistrRef,
                         std::vector<TpcLaserGridDistribution> *timeDriftDistrRef,
                         std::vector<std::vector<double>> *zPosVisRef, std::vector<std::vector<int>> *zPosOkRef,
                         std::vector<std::vector<double>> *drTimeVisRef, std::vector<std::vector<int>> *drTimeOkRef);

private:
   TString                                            inputBranchName;
   TClonesArray                                      *hits;
   BaseTpcSectorGeo                                  *secGeo;
   std::vector<std::vector<TpcLaserGridDistribution>> hit_z_distr;
   std::vector<std::vector<TpcLaserGridDistribution>> time_drift_distr;
   std::string                                        mapFile;

   double driftDistance;
   double driftTheorVelocity;
   int    sectorsNum;

   int                                           events_num;
   std::vector<std::vector<std::vector<double>>> map_vel;
   std::vector<std::vector<std::vector<double>>> map_corr_coeff;
   std::vector<std::vector<std::vector<int>>>    map_zPos_ok;
   std::vector<std::vector<std::vector<double>>> map_zPos_vis;
   std::vector<std::vector<std::vector<int>>>    map_drTime_ok;
   std::vector<std::vector<std::vector<double>>> map_drTime_vis;

   TH1D *vel_distr;
   TH1D *z_corr_coef_distr;
   TH1D *timeOffset;
   TH1D *zOffset;

   bool useFilter = kFALSE;

   bool makeQA = kFALSE;

   double lin_fit_bmult = 1.0;

   const int           wl_data_nsize = 64;
   const int           wl_data_ncut  = 15;
   void                wl_clear(std::vector<double> &data);
   std::pair<int, int> wl_filter(TpcLaserGridDistribution &distr, const int peak_bin); // return filtered bin range
   int                 FindPeak(TpcLaserGridDistribution &d, const int low_ibin, const int high_ibin);

   bool   writeToDB = kFALSE;
   size_t periodID  = 0;
   size_t runID     = 0;
   void   WriteMap(bool useDB);

   Int_t                    threadPoolSize = 24;
   std::vector<std::thread> threadPool;
   // std::mutex z_corr_coef_distr_lock;
   // std::mutex vel_distr_lock;

   clock_t clStart  = 0;
   clock_t clFinish = 0;
   clock_t clAll    = 0;

   ClassDef(TpcLaserGridVelocityMakeMap, 0)
};

#endif
