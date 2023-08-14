#include "TpcLaserGridVelocityMakeMap.h"
#include "MpdTpcHit.h"

#include <FairRootManager.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
//#include <TObject.h>
//#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
//#include <TSpline.h>
#include <TVector3.h>
#include <TXMLEngine.h>

#include "TpcPoint.h"
#include "TpcGas.h"
#include "TaskHelpers.h"

#include "TpcLaserGridVelocityDbUtils.h"
#include "curve_fit.h"
//#include "spline.h"

#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include <boost/math/interpolators/barycentric_rational.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>

TpcLaserGridVelocityMakeMap::TpcLaserGridVelocityMakeMap(BaseTpcSectorGeo &sectorGeo)
   : FairTask("TPC laser grid drift velocity map calculartor"), hits(nullptr),
     // z_layers(num_layers),
     events_num(0), vel_distr(nullptr), z_corr_coef_distr(nullptr), secGeo(&sectorGeo)
{
   inputBranchName = "TpcRecPoint";
   mapFile         = "MPDTpcVelocityMap.xml";

   z_layers.resize(num_layers);
   for (int j = 0; j < num_layers; ++j) {
      z_layers[j] = (j + 1) * z_pos_offset;
      // std::cout << "Z[" << j << "] = " << z_layers[j] << std::endl;
   }
}

//---------------------------------------------------------------------------
TpcLaserGridVelocityMakeMap::~TpcLaserGridVelocityMakeMap()
{
   /*for(int s = 0; s < sectorsNum; ++s)
   {
   }*/
}

//---------------------------------------------------------------------------
void TpcLaserGridVelocityMakeMap::UseZY(bool opt)
{
   if (opt) {
      loc_y_offset = {10., 20., 35., 85.};
   } else {
      loc_y_offset = {85.};
   }
}

//---------------------------------------------------------------------------
InitStatus TpcLaserGridVelocityMakeMap::Init()
{
   FairRootManager *ioman = FairRootManager::Instance();

   if (!ioman) {
      std::cout << "\n-E- [TpcLaserGridVelocityMakeMap::Init]: RootManager not instantiated!" << std::endl;
      return kFATAL;
   }
   hits = (TClonesArray *)ioman->GetObject(inputBranchName);

   driftDistance = secGeo->GetDriftLength(); // 163.879 cm
   sectorsNum    = secGeo->GetSectorCount(); // 24

   map_vel.resize(sectorsNum,
                  std::vector<std::vector<double>>(loc_y_offset.size(), std::vector<double>(num_layers + 3, 0.)));
   map_corr_coeff.resize(
      sectorsNum, std::vector<std::vector<double>>(loc_y_offset.size(), std::vector<double>(num_layers + 3, 0.)));

   map_zPos_vis.resize(sectorsNum,
                       std::vector<std::vector<double>>(num_layers, std::vector<double>(loc_y_offset.size(), 0)));
   map_zPos_ok.resize(sectorsNum, std::vector<std::vector<int>>(num_layers, std::vector<int>(loc_y_offset.size(), 0)));

   map_drTime_vis.resize(sectorsNum,
                         std::vector<std::vector<double>>(num_layers, std::vector<double>(loc_y_offset.size(), 0)));
   map_drTime_ok.resize(sectorsNum,
                        std::vector<std::vector<int>>(num_layers, std::vector<int>(loc_y_offset.size(), 0)));

   std::string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   TpcGas *gas        = new TpcGas(tpcGasFile, 130);
   driftTheorVelocity = gas->VDrift() * 1000; // 5.5 cm/mus
   delete gas;

   int    n_bins = ceil(driftDistance * 10);
   double length = n_bins * 0.1;
   hit_z_distr.resize(sectorsNum, std::vector<TpcLaserGridDistribution>(loc_y_offset.size(),
                                                                        TpcLaserGridDistribution(0., length, n_bins)));
   double time = length / driftTheorVelocity;
   time_drift_distr.resize(sectorsNum, std::vector<TpcLaserGridDistribution>(
                                          loc_y_offset.size(), TpcLaserGridDistribution(0., time, n_bins)));

   // vel_distr = new TH1D("tpc_vel_distr", "Distribution of electron drift velocities; Velocity [cm/#mus]; Beams
   // [count]", 600, 5.2, 5.8); z_corr_coef_distr = new TH1D("tpc_z_corr_coef_distr", "Distribution of z correction
   // coefficient; coef; Beams [count]", 4000, 0.9, 1.1);
   vel_distr         = new TH1D("vel_distr", "Velocity; Velocity [cm/#mus]; Beams [count]", 100, 5.25, 5.75);
   z_corr_coef_distr = new TH1D("z_corr_coef", "Z-correction coefficient; coef; Beams [count]", 100, 0.9, 1.1);
   timeOffset        = new TH1D("timeOffset", "Time offset; Time [#mus]; Beams [count]", 110, -2., 2.);
   zOffset           = new TH1D("zOffset", "Z offest; Z[cm]; Beams [count]", 110, -5.5, 5.5);

   threadPool.reserve(threadPoolSize);

   std::cout << "-I- TpcLaserGridVelocityMakeMap: Initialization successful." << std::endl;
   return kSUCCESS;
}

//---------------------------------------------------------------------------

void TpcLaserGridVelocityMakeMap::Exec(Option_t *opt)
{
   clStart = clock();

   events_num++;
   std::cout << "TpcLaserGridVelocityMakeMap::Exec started, event #" << events_num << std::endl;

   int nHits = hits->GetEntriesFast();
   std::cout << "Number of TPC hits is " << nHits << std::endl;

   for (int s = 0; s < sectorsNum; ++s)
      for (int y = 0; y < loc_y_offset.size(); ++y) {
         // hit_z_distr[s][y].distr.clear();
         // time_drift_distr[s][y].distr.clear();
         std::fill(hit_z_distr[s][y].distr.begin(), hit_z_distr[s][y].distr.end(), 0.);
         std::fill(time_drift_distr[s][y].distr.begin(), time_drift_distr[s][y].distr.end(), 0.);
      }

   for (int h = 0; h < nHits; ++h) {
      MpdTpcHit *hit = (MpdTpcHit *)hits->UncheckedAt(h);
      if (hit == NULL) continue;
      TVector3 gpos, lpos;
      hit->Position(gpos);
      hit->LocalPosition(lpos);
      double t = (hit->GetDriftTime()) / 1000.; // ns -> mus
      // double t = (hit->GetDriftTime() + 48.) / 1000.; // ns -> mus
      int sectNo = secGeo->SectorNumberFromGlobal(gpos);
      if (sectNo > 23 || sectNo < 0) std::cout << "WARN Sect No " << sectNo << std::endl;
      for (int y = 0; y < loc_y_offset.size(); ++y) {
         if (lpos.Y() < loc_y_offset[y]) {
            hit_z_distr[sectNo][y].Add(abs(gpos.Z()));
            time_drift_distr[sectNo][y].Add(abs(t));
            break;
         }
      }
   }

   if (makeQA) {
      if (events_num == 1) {
         toDirectory("QA/TPC/Laser");
         TH1D *g = new TH1D("z_data_test", "", hit_z_distr[2][0].distr.size(), hit_z_distr[2][0].minval,
                            hit_z_distr[2][0].maxval);
         for (int i = 0; i < hit_z_distr[2][0].distr.size(); i++) {
            g->SetBinContent(i + 1, hit_z_distr[2][0].distr[i]);
         }
         g->Write();
         g = new TH1D("time_data_test", "", time_drift_distr[2][0].distr.size(), time_drift_distr[2][0].minval,
                      time_drift_distr[2][0].maxval);
         for (int i = 0; i < time_drift_distr[2][0].distr.size(); i++) {
            g->SetBinContent(i + 1, time_drift_distr[2][0].distr[i]);
         }
         g->Write();
         gFile->cd();
      }
   }

   std::cout << "TpcLaserGridVelocityMakeMap::Exec Lasers Hits sorted" << std::endl;

   if (useFilter) {
      std::cout << "TpcLaserGridVelocityMakeMap::Exec Lasers Hits filter in use" << std::endl;
   }

   /*for(int s = 0; s < sectorsNum; ++s)
   {
       threadPool.emplace_back(&TpcLaserGridVelocityMakeMap::SectorProcessing, this, s, &(hit_z_distr[s]),
   &(time_drift_distr[s]), &(map_zPos_vis[s]), &(map_zPos_ok[s]), &(map_drTime_vis[s]), &(map_drTime_ok[s]));

       if (threadPool.size() == threadPoolSize || s == sectorsNum)
       {
           for (int i = 0; i < threadPoolSize; i++)
           {
               threadPool[i].join();
           }
           threadPool.clear();
       }
   }*/

   for (int s = 0; s < sectorsNum; ++s) {
      SectorProcessing(s, &(hit_z_distr[s]), &(time_drift_distr[s]), &(map_zPos_vis[s]), &(map_zPos_ok[s]),
                       &(map_drTime_vis[s]), &(map_drTime_ok[s]));
   }

   clFinish = clock();
   clAll    = clAll + (clFinish - clStart);

   std::cout << "TpcLaserGridVelocityMakeMap::Exec finished" << std::endl;
}

void TpcLaserGridVelocityMakeMap::SectorProcessing(const int secID, std::vector<TpcLaserGridDistribution> *hitZDistrRef,
                                                   std::vector<TpcLaserGridDistribution> *timeDriftDistrRef,
                                                   std::vector<std::vector<double>>      *zPosVisRef,
                                                   std::vector<std::vector<int>>         *zPosOkRef,
                                                   std::vector<std::vector<double>>      *drTimeVisRef,
                                                   std::vector<std::vector<int>>         *drTimeOkRef)
{
   std::vector<TpcLaserGridDistribution> &hitZDistr      = *hitZDistrRef;
   std::vector<TpcLaserGridDistribution> &timeDriftDistr = *timeDriftDistrRef;
   std::vector<std::vector<double>>      &zPosVis        = *zPosVisRef;
   std::vector<std::vector<int>>         &zPosOk         = *zPosOkRef;
   std::vector<std::vector<double>>      &drTimeVis      = *drTimeVisRef;
   std::vector<std::vector<int>>         &drTimeOk       = *drTimeOkRef;

   std::vector<std::vector<int>> timePeaks, zPeaks;
   zPeaks.resize(loc_y_offset.size(), std::vector<int>());
   timePeaks.resize(loc_y_offset.size(), std::vector<int>());
   for (int y = 0; y < loc_y_offset.size(); ++y) {
      // std::cout << "TpcLaserGridVelocityMakeMap::SectorProcessing find peaks " << secID << " y-layer " << y <<
      // std::endl; zPeaks[y] = {300, 600, 900, 1200}; timePeaks[y] = {498, 798, 1098, 1398}; zPeaks[y] =
      // hitZDistr[y].findPeaks(); timePeaks[y] = timeDriftDistr[y].findPeaks();
      zPeaks[y]    = hitZDistr[y].findPeaksSimple();
      timePeaks[y] = timeDriftDistr[y].findPeaksSimple();
      // std::cout << "TpcLaserGridVelocityMakeMap::SectorProcessing z-peaks " << zPeaks[y].size() << " t-peaks " <<
      // timePeaks[y].size() << std::endl;
   }
   for (int y = 0; y < loc_y_offset.size(); ++y) {
      // std::cout << "TpcLaserGridVelocityMakeMap::SectorProcessing process sector " << secID << " y-layer " << y <<
      // std::endl;
      if (zPeaks[y].size() == num_layers) {
         for (int l = 0; l < num_layers; ++l) {
            // int lowBinHit = hitZDistr[y].FindBin(z_layers[l] - eps_z);
            // int highBinHit = hitZDistr[y].FindBin(z_layers[l] + eps_z);

            int zEpsBins   = hitZDistr[y].FindBin(hitZDistr[y].minval + eps_z) - 1;
            int lowBinHit  = zPeaks[y][l] - zEpsBins > 0 ? zPeaks[y][l] - zEpsBins : 0;
            int highBinHit = zPeaks[y][l] + zEpsBins < hitZDistr[y].distr.size() ? zPeaks[y][l] + zEpsBins
                                                                                 : hitZDistr[y].distr.size() - 1;

            // std::cout << "lhBin = " << lowBinHit << ":" << highBinHit << std::endl;

            if (useFilter) {
               std::tie(lowBinHit, highBinHit) = wl_filter(hitZDistr[y], zPeaks[y][l]);
            }

            // std::cout << "lhBin = " << lowBinHit << ":" << highBinHit << std::endl;

            std::vector<double> xs, ys;
            std::vector<double> z_avg;
            int                 z_ok, z_info;

            std::tie(xs, ys) = hitZDistr[y].GetData(lowBinHit, highBinHit);
            // std::cout << "!!! xs = " << xs.size() << " ys = " << ys.size() << "  cbin = " << zPeaks[y][l] <<
            // std::endl;
            std::tie(z_avg, z_ok, z_info) =
               curve_fit(fit_func_gaussian, {150.0, hitZDistr[y].GetParam(zPeaks[y][l]), 0.2}, xs, ys);
            // std::cout << "!!!" << std::endl;

            if (z_info == 1 && z_ok == 0 && z_avg[0] > 40.) {
               zPosVis[l][y] += z_avg[1] * z_avg[1];
               ++zPosOk[l][y];
            }
         }
      } else {
         std::cout << "WARN:TpcLaserGridVelocityMakeMap::SectorProcessing ccf skip sector " << secID << " y-layer " << y
                   << std::endl;
      }

      if (timePeaks[y].size() == num_layers) {
         for (int l = 0; l < num_layers; ++l) {
            // int lowBinTime = timeDriftDistr[y].FindBin((fDriftDistance - z_layers[l] - eps_z) / fDriftTheorVelocity);
            // int highBinTime = timeDriftDistr[y].FindBin((fDriftDistance - z_layers[l] + eps_z) /
            // fDriftTheorVelocity);

            int timeEpsBins = timeDriftDistr[y].FindBin(timeDriftDistr[y].minval + eps_z / driftTheorVelocity) - 1;
            int lowBinTime  = timePeaks[y][l] - timeEpsBins > 0 ? timePeaks[y][l] - timeEpsBins : 0;
            int highBinTime = timePeaks[y][l] + timeEpsBins < timeDriftDistr[y].distr.size()
                                 ? timePeaks[y][l] + timeEpsBins
                                 : timeDriftDistr[y].distr.size() - 1;

            if (useFilter) {
               std::tie(lowBinTime, highBinTime) = wl_filter(timeDriftDistr[y], timePeaks[y][l]);
            }

            std::vector<double> xs, ys;
            std::vector<double> t_avg;
            int                 t_ok, t_info;

            std::tie(xs, ys) = timeDriftDistr[y].GetData(lowBinTime, highBinTime);
            std::tie(t_avg, t_ok, t_info) =
               curve_fit(fit_func_gaussian, {150.0, timeDriftDistr[y].GetParam(timePeaks[y][l]), 0.1}, xs, ys);

            if (t_info == 1 && t_ok == 0 && t_avg[0] > 40.) {
               drTimeVis[l][y] += 1. / t_avg[1];
               ++drTimeOk[l][y];
            }
         }
      } else {
         std::cout << "WARN:TpcLaserGridVelocityMakeMap::SectorProcessing vel skip sector " << secID << " y-layer " << y
                   << std::endl;
      }
   }
}

void TpcLaserGridVelocityMakeMap::Finish()
{
   std::cout << "TpcLaserGridVelocityMakeMap::Finish" << std::endl;
   for (int s = 0; s < sectorsNum; ++s)
      for (int y = 0; y < loc_y_offset.size(); ++y) {
         for (int l = 0; l < num_layers; ++l) {
            map_zPos_vis[s][l][y] = std::sqrt(map_zPos_vis[s][l][y] / map_zPos_ok[s][l][y]);
            // std::cout << "s = " << s << " y = " << y << " l = " << l << " zpos = " << map_zPos_vis[s][l][y] <<
            // std::endl;
            map_drTime_vis[s][l][y] = map_drTime_ok[s][l][y] / map_drTime_vis[s][l][y];
         }
      }

   for (int s = 0; s < sectorsNum; ++s)
      for (int y = 0; y < loc_y_offset.size(); ++y) {
         std::vector<double> x_spl(num_layers - 1), y_spl(num_layers - 1);
         std::vector<double> vel_spl(num_layers - 1);
         for (int l = 0; l < num_layers - 1; ++l) {
            x_spl[l] = abs(z_layers[l + 1] - z_layers[l]) / 2. + abs(z_layers[l]);
            y_spl[l] = abs(z_layers[l + 1] - z_layers[l]) / (map_zPos_vis[s][l + 1][y] - map_zPos_vis[s][l][y]);
            // std::cout << "s = " << s << " y = " << y << " l = " << l << " x_spl = " << x_spl[l] << " y_spl = " <<
            // y_spl[l] << std::endl;
            vel_spl[l] = abs(z_layers[l + 1] - z_layers[l]) /
                         (map_drTime_vis[s][num_layers - 1 - l][y] - map_drTime_vis[s][num_layers - 2 - l][y]);
         }
         /*tk::spline sp_xy(x_spl, y_spl);
         tk::spline sp_xVel(x_spl, vel_spl);

         map_corr_coeff[s][y][0] = sp_xy(0.);
         map_vel[s][y][0] = sp_xVel(0.);
         for (int l = 0; l < num_layers; ++l)
         {
             map_corr_coeff[s][y][l + 1] = sp_xy(z_layers[l]);
             map_vel[s][y][l + 1] = sp_xVel(z_layers[l]);
         }
         map_corr_coeff[s][y][num_layers + 1] = sp_xy(driftDistance);
         map_vel[s][y][num_layers + 1] = sp_xVel(driftDistance);*/

         boost::math::barycentric_rational<double> extr_xy(std::move(x_spl), std::move(y_spl), 1);
         map_corr_coeff[s][y][0] = extr_xy(0.);
         for (int l = 0; l < num_layers; ++l) {
            map_corr_coeff[s][y][l + 1] = extr_xy(z_layers[l]);
         }
         map_corr_coeff[s][y][num_layers + 1] = extr_xy(driftDistance);
         double corrCoeffLastInterval =
            extr_xy((driftDistance - z_layers[num_layers - 1]) / 2. + z_layers[num_layers - 1]);
         x_spl = extr_xy.return_x();
         y_spl = extr_xy.return_y();

         boost::math::barycentric_rational<double> extr_xVel(std::move(x_spl), std::move(vel_spl), 1);
         map_vel[s][y][0] = extr_xVel(0.);
         for (int l = 0; l < num_layers; ++l) {
            map_vel[s][y][l + 1] = extr_xVel(z_layers[l]);
         }
         map_vel[s][y][num_layers + 1] = extr_xVel(driftDistance);
         double velLastInterval = extr_xVel((driftDistance - z_layers[num_layers - 1]) / 2. + z_layers[num_layers - 1]);
         x_spl                  = extr_xVel.return_x();
         vel_spl                = extr_xVel.return_y();

         // trigger delay offset
         // double corrCoeffLastInterval = sp_xy((driftDistance - z_layers[num_layers - 1]) / 2. + z_layers[num_layers -
         // 1]); std::cout << "s = " << s << " y = " << (driftDistance - z_layers[num_layers - 1]) / 2. +
         // z_layers[num_layers - 1] << " corrCoeffLastInterval = " << corr_coeff_spl << std::endl;
         // map_corr_coeff[s][y][num_layers + 2] = z_layers[num_layers - 1] / corrCoeffLastInterval -
         // map_zPos_vis[s][num_layers - 1][y];
         map_corr_coeff[s][y][num_layers + 2] = z_layers[num_layers - 1] - map_zPos_vis[s][num_layers - 1][y];

         // double velLastInterval = sp_xVel((driftDistance - z_layers[num_layers - 1]) / 2. + z_layers[num_layers -
         // 1]);
         map_vel[s][y][num_layers + 2] =
            (driftDistance - z_layers[num_layers - 1]) / velLastInterval - map_drTime_vis[s][0][y];
      }

   WriteMap(writeToDB);

   std::cout << "TpcLaserGridVelocityMakeMap work time = " << ((Float_t)clAll) / CLOCKS_PER_SEC << std::endl;

   if (makeQA) {
      for (int s = 0; s < sectorsNum; ++s)
         for (int y = 0; y < loc_y_offset.size(); ++y) {
            for (int l = 0; l < num_layers + 2; ++l) {
               vel_distr->Fill(map_vel[s][y][l]);
               z_corr_coef_distr->Fill(map_corr_coeff[s][y][l]);
            }

            timeOffset->Fill(map_vel[s][y][num_layers + 2]);
            zOffset->Fill(map_corr_coeff[s][y][num_layers + 2]);
         }

      toDirectory("QA/TPC/Laser");
      z_corr_coef_distr->Write(0, kOverwrite);
      timeOffset->Write(0, kOverwrite);
      zOffset->Write(0, kOverwrite);

      TF1 fit_vel("fit_vel", "gaus", 5.3, 5.7);
      fit_vel.SetParLimits(0, 10., 1000.);
      fit_vel.SetParLimits(1, 5.3, 5.6);
      fit_vel.SetParLimits(2, 0.005, 0.04);
      // fit_vel.SetParameter(2, 0.0025);
      vel_distr->Fit(&fit_vel, "BLMQR");
      // vel_distr->Fit("gaus", "WWQR", "", 5.3, 5.7);
      vel_distr->Write(0, kOverwrite);
   }
   gFile->cd();

   // fHits->Compress();

   std::cout << "TpcLaserGridVelocityMakeMap::Finished" << std::endl;
}

void TpcLaserGridVelocityMakeMap::WriteMap(bool useDB)
{
   if (useDB) {
      TpcLaserGridVelocityDbUtils::setTpcVelocityMap(periodID, runID, map_corr_coeff);
      return;
   }

   if (mapFile.empty()) {
      Error("TpcLaserGridVelocityMakeMap::WriteFile", "file not specified");
      return;
   }

   TXMLEngine      *xml      = new TXMLEngine();
   XMLNodePointer_t mainnode = xml->NewChild(0, 0, "MPDTpcDriftVelocityMap");

   for (int iSect = 0; iSect < sectorsNum; ++iSect) {
      XMLNodePointer_t sector = xml->NewChild(mainnode, 0, "TpcSector");
      xml->NewAttr(sector, 0, "n", (std::to_string(iSect)).c_str());

      for (int ya = 0; ya < loc_y_offset.size(); ++ya) {
         XMLNodePointer_t yarea = xml->NewChild(sector, 0, "yarea");
         xml->NewAttr(yarea, 0, "n", (std::to_string(ya)).c_str());

         for (int zl = 0; zl < num_layers + 3; ++zl) {
            XMLNodePointer_t zlayer = xml->NewChild(yarea, 0, "zlayer");
            xml->NewAttr(zlayer, 0, "n", (std::to_string(zl)).c_str());
            xml->NewChild(zlayer, 0, "vel", (std::to_string(map_vel[iSect][ya][zl])).c_str());
            xml->NewChild(zlayer, 0, "coeff", (std::to_string(map_corr_coeff[iSect][ya][zl])).c_str());
         }
      }
   }

   XMLDocPointer_t xmldoc = xml->NewDoc();
   xml->DocSetRootElement(xmldoc, mainnode);
   xml->SaveDoc(xmldoc, mapFile.c_str());
   xml->FreeDoc(xmldoc);
   delete xml;
}

void TpcLaserGridVelocityMakeMap::wl_clear(std::vector<double> &data)
{
   std::vector<double> abs_coeffs(wl_data_nsize);
   std::vector<size_t> sort_indexes(wl_data_nsize);

   gsl_wavelet           *wl    = gsl_wavelet_alloc(gsl_wavelet_haar_centered, 2);
   gsl_wavelet_workspace *wl_ws = gsl_wavelet_workspace_alloc(wl_data_nsize);

   gsl_wavelet_transform_forward(wl, data.data(), 1, wl_data_nsize, wl_ws);

   for (int i = 0; i < wl_data_nsize; ++i) abs_coeffs[i] = fabs(data[i]);

   gsl_sort_index(sort_indexes.data(), abs_coeffs.data(), 1, wl_data_nsize);

   for (int i = 0; (i + wl_data_ncut) < wl_data_nsize; ++i) data[sort_indexes[i]] = 0;

   gsl_wavelet_transform_inverse(wl, data.data(), 1, wl_data_nsize, wl_ws);

   double ar_avg = 0.;
   int    n      = 0;
   for (int i = 0; i < wl_data_nsize; i++)
      if (fabs(data[i]) < 10.) {
         ar_avg += fabs(data[i]);
         n++;
      }
   ar_avg /= (double)n;
   // std::cout << "avd_val = " << ar_avg << std::endl;
   ar_avg *= lin_fit_bmult;

   for (int i = 0; i < wl_data_nsize; ++i) data[i] = (fabs(data[i]) - ar_avg) > 0. ? fabs(data[i]) - ar_avg : 0.;

   gsl_wavelet_free(wl);
   gsl_wavelet_workspace_free(wl_ws);
}

std::pair<int, int> TpcLaserGridVelocityMakeMap::wl_filter(TpcLaserGridDistribution &distr, const int peak_bin)
{
   std::vector<double> wl_data(wl_data_nsize);

   int low_ibin_insert = peak_bin - (wl_data_nsize / 2);
   if (low_ibin_insert < 0) low_ibin_insert = 0;
   std::copy(distr.distr.begin() + low_ibin_insert, distr.distr.begin() + low_ibin_insert + wl_data_nsize,
             wl_data.begin());

   wl_clear(wl_data);

   std::copy(wl_data.begin(), wl_data.begin() + wl_data_nsize, distr.distr.begin() + low_ibin_insert);

   return std::make_pair(low_ibin_insert, low_ibin_insert + wl_data_nsize);
}

int TpcLaserGridVelocityMakeMap::FindPeak(TpcLaserGridDistribution &d, const int low_ibin, const int high_ibin)
{
   std::vector<double>::iterator it_max_val = std::max_element(d.distr.begin() + low_ibin, d.distr.begin() + high_ibin);
   return std::distance(d.distr.begin(), it_max_val);
}

ClassImp(TpcLaserGridVelocityMakeMap)
