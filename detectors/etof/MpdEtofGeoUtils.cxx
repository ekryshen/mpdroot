//------------------------------------------------------------------------------------------------------------------------
#include <assert.h>
#include <iostream>

#include <TGeoManager.h>
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"

#include "FairRootManager.h"
#include "FairLogger.h"

#include "MpdTofPoint.h"
#include "MpdEtofGeoUtils.h"
using namespace std;
MpdEtofGeoUtils *MpdEtofGeoUtils::instance = nullptr;
//------------------------------------------------------------------------------------------------------------------------
MpdEtofGeoUtils::MpdEtofGeoUtils() {}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofGeoUtils::FindNeighborStrips(Double_t thresh, bool forced)
{
   /*	static bool founded = false;
      if( !forced &&  founded) return; // neighbor strips already founded and filled

      size_t NR = 0, NL= 0;
      const LStrip *strip2; double  distance;

      for(auto it1 = mStrips.begin(), itEnd1 = mStrips.end(); it1 != itEnd1 ; it1++) // cycle1 by strips
      {
         LStrip *strip1 = &(it1->second);

         for(auto it2 = mStrips.begin(), itEnd2 = mStrips.end(); it2 != itEnd2 ; it2++) // cycle2 by strips
         {
            strip2 = &(it2->second);

            // CATION: Ckeck  only left and right sides(one row strips NOW)
            distance = strip1->Distance(LStrip::kRight, *strip2);
            if(h1)  h1->Fill(distance);
            if(distance < thresh) // CAUTION: constant depends on the geometry layout(see h1TestDistance histo)
            {
               strip1->neighboring[LStrip::kRight] = strip2->volumeUID;
               if(h2) h2->Fill(strip1->stripID, strip2->stripID);
            }

            distance = strip1->Distance(LStrip::kLeft, *strip2);
            if(h1)  h1->Fill(distance);
            if(distance < thresh) // CAUTION: constant depends on the geometry layout(see h1TestDistance histo)
            {
               strip1->neighboring[LStrip::kLeft] = strip2->volumeUID;
               if(h2) h2->Fill( strip2->stripID, strip1->stripID);
            }

         }// cycle2 by strips

      }// cycle1 by strips

      founded = true;

      cout<<" [MpdEtofGeoUtils::FindNeighborStrips] Neighbor strips: left = "<<NL<<", right = "<<NR<<endl;
   */
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofGeoUtils::ParseTGeoManager(bool forced)
{
   assert(gGeoManager);
   // static const double degree180 = 3.14159265359; // 180 degree to rad

   if (!forced && !mStrips.empty()) return; // already parsed and filled
   /*
      mStrips.clear();

      int nBoxes = 0, nStrips = 0, nGases = 0;
      const char* pathToTOF[] = {"/cave_1/etof1_1", "/cave_1/etof1_2"};

      for(size_t side = 0, sectorID = 1; side <= 1 ; side++,  sectorID++)
      {
         gGeoManager->cd(pathToTOF[side]);

         TGeoNode *boxNode, *stripNode, *gasNode;
         Int_t volumeUID, boxID, stripID, gasID; // module[1,...,30], pad [1,...,131]
         vector<intervalType> vDetectorsR, vDetectorsPhi;

         Double_t *X0Y0Z0 = new Double_t[3]; X0Y0Z0[0] = X0Y0Z0[1] = X0Y0Z0[2] = 0.; // center of sensetive detector
         Double_t  *local = new Double_t[3], master[3],  dX, dY, dZ;

         TObjArray *array = gGeoManager->GetCurrentVolume()->GetNodes();
         TIterator *it1 = array->MakeIterator();

         while( (boxNode = (TGeoNode*) it1->Next()) )		// BOXES
            {
               TString PATH1 = TString(pathToTOF[side]) + "/" + boxNode->GetName(); boxID = boxNode->GetNumber();
   nBoxes++;
            ///cout<<"\n MODULE: "<<boxNode->GetName()<<", copy# "<<boxID<<" path= "<<PATH1.Data();

                  TIterator *it2 = boxNode->GetNodes()->MakeIterator();
                  while( (stripNode = (TGeoNode*) it2->Next()) )		// STRIPS
            {
               TString stripName = stripNode->GetName();
               TString PATH2 = PATH1 + "/" + stripName; nStrips++;
               TString stripid = stripName; stripid.Remove(0, 7); stripid.Remove(stripid.Last('_')); stripID =
   stripid.Atoi();
               ///cout<<"\n\t PAD: "<<stripName.Data()<<", copy# "<<stripID<<" path= "<<PATH2.Data();

               TIterator *it3 = stripNode->GetNodes()->MakeIterator();
               while( (gasNode = (TGeoNode*) it3->Next()) )	// GASES
                  {
                     TString detName = gasNode->GetName();
                     if(!detName.Contains("gas")) continue;

                     TString PATH3 = PATH2  + "/" + detName;  nGases++;
                  TString gasid = stripName; gasid.Remove(0, 7); gasid.Remove(gasid.Last('_')); gasID = gasid.Atoi();
                  ///cout<<"\n\t\t GAS: "<<detName.Data()<<", copy# "<<gasID<<" path= "<<PATH3.Data();

                  gGeoManager->cd(PATH3);

                  TGeoMatrix *matrix = gGeoManager->GetCurrentMatrix();	// calculate global TGeoHMatrix for current node
                  matrix->LocalToMaster(X0Y0Z0, master);			// 0.0.0. --> MRS

                  TGeoBBox *box = (TGeoBBox*) gGeoManager->GetCurrentNode()->GetVolume()->GetShape();
                  dX = box->GetDX(); dY = box->GetDY(); dZ = box->GetDZ();

                  volumeUID = MpdTofPoint::GetSuid(sectorID, boxID, 1, stripID); // FIXME: now ONE detector per
   box(detector=1)
                  ///cout<<"\n uid="<<volumeUID<<", center= "<<master[0]<<", "<<master[1]<<", "<<master[2]<<", dXYZ=
   "<<dX<<", "<<dY<<", "<<dZ;

                  LStrip stripData(volumeUID, sectorID, boxID, 1, stripID);
                  stripData.center.SetXYZ(master[0], master[1], master[2]);

                  // point A
                  local[0] = -dX;	local[1] = -dY; local[2] = +dZ;
                  matrix->LocalToMaster(local, master);
                  stripData.A.SetXYZ(master[0], master[1], master[2]);

                  // point B
                  local[0] = +dX;	local[1] = -dY; local[2] = +dZ;
                  matrix->LocalToMaster(local, master);
                  stripData.B.SetXYZ(master[0], master[1], master[2]);

                  // point C
                  local[0] = +dX;	local[1] = -dY; local[2] = -dZ;
                  matrix->LocalToMaster(local, master);
                  stripData.C.SetXYZ(master[0], master[1], master[2]);

                  // point D
                  local[0] = -dX;	local[1] = -dY; local[2] = -dZ;
                  matrix->LocalToMaster(local, master);
                  stripData.D.SetXYZ(master[0], master[1], master[2]);

                  stripData.InitCenterPerp();

                  if(h2)
                  {
                     h2->Fill(stripData.A.Perp(), stripData.A.Phi());
                     h2->Fill(stripData.B.Perp(), stripData.B.Phi());
                     h2->Fill(stripData.C.Perp(), stripData.C.Phi());
                     h2->Fill(stripData.D.Perp(), stripData.D.Phi());
                  }
   ///stripData.Dump("\n ETOFstrip:");
                  bool IsUniqueUID = mStrips.insert(make_pair(volumeUID, stripData)).second;
   assert(IsUniqueUID);

                  // Add strip location interval to intervalTree
                  LRectangle * pStrip = new LRectangle(stripData);
                  pStrip->InitCenterPerp();

                  double Rmin, Rmax, Phimin, Phimax;
                  pStrip->GetRPhiRanges(Rmin, Rmax, Phimin, Phimax);

                  vDetectorsR.push_back(intervalType(Rmin, Rmax, pStrip));

                  if(Phimax - Phimin > degree180) // segment overlap the boundary [-180.,+180[
                  {
                     vDetectorsPhi.push_back(intervalType(-degree180, Phimin, pStrip));
                     vDetectorsPhi.push_back(intervalType(Phimax,  degree180, pStrip));
                     ///cout<<"\n OVERLAPED   phimin="<<Phimin<<", phimax="<<Phimax;
                  }
                  else 	vDetectorsPhi.push_back(intervalType(Phimin, Phimax, pStrip));

                  } // GASES

            } // STRIPS

            } // MODULES

            mDetectorsR[side] = intervalTreeType(vDetectorsR);
            mDetectorsPhi[side] = intervalTreeType(vDetectorsPhi);
      }

         FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "[MpdEtofGeoUtils::ParseTGeoManager] Regions= %d, modules= %d,
   strips= %d. ", 2, nBoxes, nStrips);
   */
}
//------------------------------------------------------------------------------------------------------------------------
const LStrip *MpdEtofGeoUtils::FindStrip(Int_t UID)
{
   auto cit = mStrips.find(UID);
   assert(cit != mStrips.end());
   return &(cit->second);
}
//----------------------------------------------------------------------------------------------------------------------