////////////////////////////////////////////////////////////////
//                                                            //
//  MpdEmcGeoParams                                           //
//  EMC geometry in MpdEmcHitCreation, v03                    //
//  Author List : Martemianov M., 2019 		              //
//                                                            //
////////////////////////////////////////////////////////////////

#include "MpdEmcGeoParams.h"
#include "FairParamList.h"
#include "TObjArray.h"
#include <iostream>

#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoPgon.h"
#include "TGeoTube.h"

ClassImp(MpdEmcGeoParams)

   MpdEmcGeoParams ::MpdEmcGeoParams(const char *name, const char *title, const char *context)
   : FairParGenericSet(name, title, context)
{

   fGeoSensNodes = new TObjArray();
   fGeoPassNodes = new TObjArray();
}

MpdEmcGeoParams::MpdEmcGeoParams()
   : fSectorNumber(25), fCrateNumber(6), fModuleTypes(8), fTowerNumber(38400 / 2), fTowerNumberXY(2), fTowerNumberZ(64),
     fTowerLength(41.55)
{

   TString pathChamber = "/cave_1/emcChamber_0";
   gGeoManager->cd(pathChamber + "/emcChH_0");

   // Total barell sizes

   length = 2. * ((TGeoTube *)gGeoManager->GetCurrentVolume()->GetShape())->GetDz();
   rMin   = ((TGeoTube *)gGeoManager->GetCurrentVolume()->GetShape())->GetRmin();
   rMax   = ((TGeoTube *)gGeoManager->GetCurrentVolume()->GetShape())->GetRmax();

   TGeoMatrix     *rowMatrix;
   TString         pathTower, pathCh;
   TGeoNode       *towerNode;
   const Double_t *fRotMatrix;
   const Double_t *fTransMatrix;
   Double_t        master[3], phiAngle, rotAngle, fAngleCrate;

   Int_t fTowerInCrate = fTowerNumberXY * fTowerNumberZ;
   Int_t towerNumber = -1, towerXYNum, towerZNum, iTowerElem = 0;
   Int_t iSector, iCrate, iModule = 0, iTowerNum = 0, iTowerNode = 0, iTowerRange = 0;

   pathCh = pathChamber + "/emcChH_0";

   TGeoVolume *emcRow = (TGeoVolume *)gGeoManager->GetVolume("emcCrate");
   fAngleCrate        = ((TGeoPgon *)emcRow->GetShape())->GetDphi();

   for (Int_t iChamber = 0; iChamber < 2; iChamber++) {
      if (iChamber == 1) pathCh = pathChamber + "/emcChH_1";

      for (Int_t iTower = 0; iTower < fTowerNumber; iTower++) {

         iSector = iTower / (fCrateNumber * fTowerInCrate);
         iCrate  = (iTower - iSector * fCrateNumber * fTowerInCrate) / fTowerInCrate;

         if (iTowerNode == fTowerInCrate) iTowerNode = 0;
         if (iTowerRange % 2 == 0) {
            iTowerRange = 0;
            iTowerNum++;
         }
         if (iTowerNum == fTowerNumberZ + 1) iTowerNum = 1;
         iModule = iTowerNode / (2 * fModuleTypes);
         pathTower.Form("%s/emcSector_%d/emcCrate_%d/emcModule%d_0/emc_box%d_%d/", pathCh.Data(), iSector, iCrate,
                        iModule, iTowerNum, iTowerNode);

         towerXYNum = iSector * fCrateNumber * fTowerNumberXY + 2 * iCrate + iTowerNode % 2;
         towerZNum  = fTowerNumberZ + iTowerNum - 1;
         if (iChamber == 1) towerZNum = fTowerNumberZ - iTowerNum;
         towerNumber = towerXYNum * 1000 + towerZNum;
         iTowerRange++;
         iTowerNode++;

         pointID.insert(pair<Int_t, Int_t>(towerNumber, iTowerElem));

         gGeoManager->cd(pathTower);

         rowMatrix = gGeoManager->GetCurrentMatrix();
         towerNode = gGeoManager->GetCurrentNode();
         rowMatrix->LocalToMaster(towerNode->GetMatrix()->GetTranslation(), master);
         fRotMatrix   = rowMatrix->GetRotationMatrix();
         fTransMatrix = rowMatrix->GetTranslation();
         phiAngle     = ATan2(fRotMatrix[3], fRotMatrix[0]) * RadToDeg();
         if (phiAngle < 0) phiAngle += 360;
         phiTower.push_back(phiAngle);
         thetaTower.push_back(ACos(fRotMatrix[8]) * RadToDeg());
         xTower.push_back(fTransMatrix[0]);
         yTower.push_back(fTransMatrix[1]);
         zTower.push_back(fTransMatrix[2]);
         rhoTower.push_back(sqrt(fTransMatrix[0] * fTransMatrix[0] + fTransMatrix[1] * fTransMatrix[1] +
                                 fTransMatrix[2] * fTransMatrix[2]));
         iTowerElem++;
      }
   }

   for (Int_t iCrt = 0; iCrt < fSectorNumber * fCrateNumber; iCrt++) {
      rotAngle = fAngleCrate * (iCrt + 2);
      if (rotAngle > 360) rotAngle = rotAngle - 360.;
      phiRow.push_back(rotAngle);
   }
}

MpdEmcGeoParams::~MpdEmcGeoParams(void) {}

void MpdEmcGeoParams::clear(void)
{
   if (fGeoSensNodes) delete fGeoSensNodes;
   if (fGeoPassNodes) delete fGeoPassNodes;
}

void MpdEmcGeoParams::putParams(FairParamList *l)
{
   if (!l) return;
   l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
   l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}

Bool_t MpdEmcGeoParams::getParams(FairParamList *l)
{
   if (!l) return kFALSE;
   if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   if (!l->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;
   return kTRUE;
}
