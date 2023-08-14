#include "TpcLaserGridVelocityAdjByMap.h"
#include "MpdTpcHit.h"

#include <FairRootManager.h>
//#include <TObject.h>
//#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TXMLEngine.h>

#include <boost/math/interpolators/barycentric_rational.hpp>

#include "TpcLaserGridVelocityDbUtils.h"

#include "TpcGas.h"
#include "TaskHelpers.h"

#include <cmath>
#include <iostream>
#include <string>

TpcLaserGridVelocityAdjByMap::TpcLaserGridVelocityAdjByMap(BaseTpcSectorGeo &sectorGeo)
   : FairTask("TPC laser grid drift velocity map adjuster"), hits(nullptr), secGeo(&sectorGeo), events_num(0)
{
   inputBranchName = "TpcRecPoint";
   mapFile         = "MPDTpcVelocityMap.xml";
}

//---------------------------------------------------------------------------
TpcLaserGridVelocityAdjByMap::~TpcLaserGridVelocityAdjByMap()
{
   for (int s = 0; s < sectorsNum; ++s) {
      for (auto ccz : corr_coeff_z[s]) gsl_spline_free(ccz);
      for (auto vz : vel_z[s]) gsl_spline_free(vz);

      gsl_interp_accel_free(accel[s]);
   }
}

//---------------------------------------------------------------------------
void TpcLaserGridVelocityAdjByMap::UseZY(bool opt)
{
   if (opt) {
      loc_y_offset = {10., 20., 35., 85.};
   } else {
      loc_y_offset = {85.};
   }
}

//---------------------------------------------------------------------------
InitStatus TpcLaserGridVelocityAdjByMap::Init()
{
   FairRootManager *ioman = FairRootManager::Instance();

   if (!ioman) {
      std::cout << "\n-E- [TpcLaserGridVelocityAdjByMap::Init]: RootManager not instantiated!" << std::endl;
      return kFATAL;
   }
   hits = (TClonesArray *)ioman->GetObject(inputBranchName);

   driftDistance = secGeo->GetDriftLength(); // 163.879 cm
   sectorsNum    = secGeo->GetSectorCount();

   map_vel_harm.resize(sectorsNum,
                       std::vector<std::vector<double>>(loc_y_offset.size(), std::vector<double>(num_layers + 3, 0.)));
   map_corr_coeff.resize(
      sectorsNum, std::vector<std::vector<double>>(loc_y_offset.size(), std::vector<double>(num_layers + 3, 0.)));

   int status = ReadMap(readFromDB);

   std::string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   TpcGas *gas        = new TpcGas(tpcGasFile, 130);
   driftTheorVelocity = gas->VDrift() * 1000; // 5.5 cm/mus
   delete gas;

   z_layers.resize(num_layers + 2);
   z_layers[0]              = 0.;
   z_layers[num_layers + 1] = driftDistance;
   for (int j = 0; j < num_layers; ++j) {
      z_layers[j + 1] = (j + 1) * z_pos_offset;
   }
   // std::cout<<z_layers[0]<<" "<<z_layers[1]<<" "<<z_layers[2]<<" "<<z_layers[3]<<" "<<z_layers[4]<<"
   // "<<z_layers[5]<<std::endl;
   /*y_layers.resize(loc_y_offset.size() + 2);
   y_layers[0] = -0.5;
   y_layers[loc_y_offset.size() + 1] = loc_y_offset[loc_y_offset.size() - 1];
   y_layers[1] = (loc_y_offset[0] - y_layers[0]) / 2.;
   for (int i = 2; i < loc_y_offset.size() + 1; ++i)
   {
       y_layers[i] = loc_y_offset[i - 2] + (loc_y_offset[i - 1] - loc_y_offset[i - 2]) / 2.;
   }*/
   // std::cout<<y_layers.size()<<" "<<y_layers[0]<<" "<<y_layers[1]<<" "<<y_layers[2]<<" "<<y_layers[3]<<"
   // "<<y_layers[4]<<" "<<y_layers[5]<<std::endl;

   corr_coeff_z.resize(sectorsNum, std::vector<gsl_spline *>(loc_y_offset.size(), nullptr));
   vel_z.resize(sectorsNum, std::vector<gsl_spline *>(loc_y_offset.size(), nullptr));
   for (int s = 0; s < sectorsNum; ++s) {
      for (int y = 0; y < loc_y_offset.size(); ++y) {
         corr_coeff_z[s][y] = gsl_spline_alloc(gsl_interp_steffen, z_layers.size());
         vel_z[s][y]        = gsl_spline_alloc(gsl_interp_steffen, z_layers.size());
         gsl_spline_init(corr_coeff_z[s][y], z_layers.data(), map_corr_coeff[s][y].data(), z_layers.size());
         gsl_spline_init(vel_z[s][y], z_layers.data(), map_vel_harm[s][y].data(), z_layers.size());
      }
   }

   accel.resize(threadPoolSize, nullptr);
   for (int th = 0; th < threadPoolSize; ++th) accel[th] = gsl_interp_accel_alloc();

   threadPool.reserve(threadPoolSize);

   std::cout << "-I- TpcLaserGridVelocityAdjByMap: Initialization successful." << std::endl;
   return kSUCCESS;
}

//---------------------------------------------------------------------------

void TpcLaserGridVelocityAdjByMap::Exec(Option_t *opt)
{
   events_num++;
   std::cout << "TpcLaserGridVelocityAdjByMap::Exec started, event #" << events_num << std::endl;

   int nHits = hits->GetEntriesFast();
   std::cout << "Number of TPC hits is " << nHits << std::endl;

   ThreadProcessing(0, nHits, accel[0]);
   /*double hitPortion = nHits / threadPoolSize;
   int low = 0, high = floor(hitPortion);
   for(int th = 0; th < threadPoolSize; ++th)
   {
       high = th == threadPoolSize - 1 ? nHits : high;
       //std::cout << "!!! th " << th << " low=" << low << " high=" << high << std::endl;
       threadPool.emplace_back(&TpcLaserGridVelocityAdjByMap::ThreadProcessing, this, low, high, accel[th]);
       low += hitPortion;
       high += hitPortion;
   }
   for (int th = 0; th < threadPoolSize; ++th)
   {
       threadPool[th].join();
   }
   threadPool.clear();*/

   std::cout << "TpcLaserGridVelocityAdjByMap::Exec finished" << std::endl;
}

void TpcLaserGridVelocityAdjByMap::ThreadProcessing(int hitIdxLow, int hitIdxHigh, gsl_interp_accel *accel)
{
   std::vector<double> y_layers(loc_y_offset.size());
   y_layers[0] = (loc_y_offset[0] - (-.5)) / 2.;
   for (int i = 1; i < loc_y_offset.size(); ++i) {
      y_layers[i] = loc_y_offset[i - 1] + (loc_y_offset[i] - loc_y_offset[i - 1]) / 2.;
   }

   for (int h = hitIdxLow; h < hitIdxHigh; ++h) {
      getTpcHit_lock.lock();
      MpdTpcHit *hit = (MpdTpcHit *)hits->UncheckedAt(h);
      getTpcHit_lock.unlock();
      if (hit == NULL) continue;
      TVector3 gpos, lpos;
      hit->Position(gpos);
      hit->LocalPosition(lpos);
      int sectNo = secGeo->SectorNumberFromGlobal(gpos);

      double z = lpos.Z();
      if (abs(z) < z_layers[0] || abs(z) > z_layers[z_layers.size() - 1]) continue;

      double corr_coeff_val = 1.;
      int    yArea          = 0;
      if (loc_y_offset.size() > 1) {
         double              y = lpos.Y();
         std::vector<double> cc(loc_y_offset.size());
         for (int yi = 0; yi < loc_y_offset.size(); ++yi) {
            cc[yi] = gsl_spline_eval(corr_coeff_z[sectNo][yi], abs(z), accel);
            if (y < loc_y_offset[yi]) yArea++;
         }
         // std::cout<<cc[0]<<" "<<cc[1]<<" "<<cc[2]<<" "<<cc[3]<<" "<<cc[4]<<" "<<cc[5]<<std::endl;

         boost::math::barycentric_rational<double> corr_coeff_y(std::move(y_layers), std::move(cc),
                                                                y_layers.size() > 3 ? 3 : y_layers.size() - 1);
         corr_coeff_val = corr_coeff_y(y);
         y_layers       = corr_coeff_y.return_x();
         cc             = corr_coeff_y.return_y();
      } else {
         corr_coeff_val = gsl_spline_eval(corr_coeff_z[sectNo][yArea], abs(z), accel);
      }

      // double z_corr = driftDistance - (driftDistance - abs(gpos.Z())) * corr_coeff_val;
      double z_corr =
         driftDistance -
         (driftDistance - (abs(gpos.Z()) + map_corr_coeff[sectNo][yArea][num_layers + 2])) * corr_coeff_val;

      hit->SetZ(std::copysign(z_corr, z));
      z      = hit->GetLocalZ();
      z_corr = driftDistance - (driftDistance - abs(z)) * corr_coeff_val;
      hit->SetLocalZ(std::copysign(z_corr, z));
   }
}

void TpcLaserGridVelocityAdjByMap::Finish()
{
   /*if (fMakeQA)
   {
       toDirectory("QA/TPC/Laser");
   }
   gFile->cd();*/

   // fHits->Compress();

   std::cout << "TpcLaserGridVelocityAdjByMap::Finished" << std::endl;
}

int TpcLaserGridVelocityAdjByMap::ReadMap(bool useDB)
{
   if (useDB) {
      std::vector<std::vector<std::vector<double>>> map_buf;
      if (TpcLaserGridVelocityDbUtils::getTpcVelocityMap(periodID, runID, map_buf)) {
         map_corr_coeff = std::move(map_buf);

         for (int iSect = 0; iSect < sectorsNum; ++iSect) {
            std::cout << "sect " << iSect << std::endl;
            for (int zl = 0; zl < num_layers; ++zl) {
               std::cout << "zlay " << zl << std::endl;
               std::cout << "ya";
               for (int ya = 0; ya < loc_y_offset.size(); ++ya) {
                  std::cout << ":" << map_corr_coeff[iSect][ya][zl];
               }
               std::cout << std::endl;
            }
         }

         return 0;
      } else {
         Error("TpcLaserGridVelocityAdjByMap::ReadMap", "DB read error");
         return 6;
      }
   }

   if (mapFile.empty()) {
      Error("TpcLaserGridVelocityAdjByMap::ReadMap", "file not specified");
      return 1;
   }

   TXMLEngine     *xml    = new TXMLEngine;
   XMLDocPointer_t xmldoc = xml->ParseFile(mapFile.c_str());
   if (xmldoc == NULL) {
      Error("TpcLaserGridVelocityAdjByMap::ReadMap", "XML Engine parse failed");
      return 2;
   }

   XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);

   for (int iSect = 0; iSect < sectorsNum; ++iSect) {
      XMLNodePointer_t sector = iSect == 0 ? xml->GetChild(mainnode) : xml->GetNext(sector);
      if (sector == NULL) {
         Error("TpcLaserGridVelocityAdjByMap::ReadMap", "XML Engine sector failed");
         return 3;
      }

      XMLAttrPointer_t ns_ptr = xml->GetFirstAttr(sector);
      double           ns     = std::stod(xml->GetAttrValue(ns_ptr));

      for (int y = 0; y < loc_y_offset.size(); ++y) {
         XMLNodePointer_t yarea = (y == 0) ? xml->GetChild(sector) : xml->GetNext(yarea);
         if (yarea == NULL) {
            Error("TpcLaserGridVelocityAdjByMap::ReadMap", "XML Engine yarea failed");
            return 4;
         }

         XMLAttrPointer_t ya_ptr = xml->GetFirstAttr(yarea);
         double           ya     = std::stod(xml->GetAttrValue(ya_ptr));

         for (int l = 0; l < num_layers + 3; ++l) {
            XMLNodePointer_t zlayer = (l == 0) ? xml->GetChild(yarea) : xml->GetNext(zlayer);
            if (zlayer == NULL) {
               Error("TpcLaserGridVelocityAdjByMap::ReadMap", "XML Engine layer failed");
               return 5;
            }

            XMLAttrPointer_t zl_ptr = xml->GetFirstAttr(zlayer);
            double           zl     = std::stod(xml->GetAttrValue(zl_ptr));

            XMLNodePointer_t zParam = xml->GetChild(zlayer);
            while (zParam != NULL) {
               std::string paramName(xml->GetNodeName(zParam));
               double      paramVal = std::stod(xml->GetNodeContent(zParam));

               if (paramName == "vel")
                  map_vel_harm[ns][ya][zl] = paramVal;
               else if (paramName == "coeff")
                  map_corr_coeff[ns][ya][zl] = paramVal;

               zParam = xml->GetNext(zParam);
            }
         }
      }
   }

   xml->FreeDoc(xmldoc);
   delete xml;

   /*for (int iSect = 0; iSect < sectorsNum; ++iSect)
   {
       std::cout << "sect " << iSect << std::endl;
       for (int zl = 0; zl < num_layers; ++zl)
       {
           std::cout << "zlay " << zl << std::endl;
           std::cout << "ya";
           for (int ya = 0; ya < loc_y_offset.size(); ++ya)
           {
               std::cout <<":"<< map_corr_coeff[iSect][ya][zl];
           }
           std::cout << std::endl;
       }
   }*/

   return 0;
}

ClassImp(TpcLaserGridVelocityAdjByMap)
