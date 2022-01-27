/*************************************************************************************
 *
 *         Class MpdZdcDigiScheme
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  16-May-2008
 *
 ************************************************************************************/

#include "MpdZdcDigiScheme.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include <iostream>

using std::cout;
using std::endl;

MpdZdcDigiScheme *MpdZdcDigiScheme::fInstance    = 0;
Bool_t            MpdZdcDigiScheme::fInitialized = 0;
Int_t             MpdZdcDigiScheme::fRefcount    = 0;
Bool_t            MpdZdcDigiScheme::fYesPsd      = 0;

static Int_t kZDC = 1; //  // hard-coded basic ZDC detector number

// -------------------------------------------------------------------------

MpdZdcDigiScheme::MpdZdcDigiScheme() : Nx(0), Ny(0), Nz(0)
{
   //  Nx=0;  Ny=0; Nz=0;
   Nx_psd         = 0;
   Ny_psd         = 0;
   Nz_psd         = 0;
   fZdcDigiPar    = 0;
   fZdcPsdDigiPar = 0;
}

// -------------------------------------------------------------------------

MpdZdcDigiScheme::~MpdZdcDigiScheme()
{
   fRefcount--;
   if (!fRefcount) {
      delete this;
      fInstance = nullptr;
   }
}
// -------------------------------------------------------------------------

MpdZdcDigiScheme *MpdZdcDigiScheme::Instance()
{
   if (!fInstance) fInstance = new MpdZdcDigiScheme();
   fRefcount++;
   return fInstance;
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::Init(MpdZdcGeoPar *geoPar, MpdZdcDigiPar *digiPar, Bool_t pAddPsd, Int_t pVerbose)
{
   if (!fInitialized) {

      if (!geoPar) {
         cout << "-W- MpdZdcDigiScheme::Init: "
              << "No geometry parameters available!" << endl;
         return kFALSE;
      }
      if (!digiPar) {
         cout << "-W-  MpdZdcDigiScheme::Init: "
              << "No digitization parameters available!" << endl;
         //      return kFALSE;
      }
      fZdcDigiPar = digiPar;

      TObjArray *sensNodes = geoPar->GetGeoSensitiveNodes();
      if (!sensNodes) {
         cout << "-W-  MpdZdcDigiScheme::Init: "
              << "No sensitive nodes available!" << endl;
         return kFALSE;
      }

      fPasNodes = geoPar->GetGeoPassiveNodes();
      if (!fPasNodes) {
         cout << "-W-  MpdZdcDigiScheme::Init: "
              << "No passive nodes available!" << endl;
         return kFALSE;
      }

      fInitialized = kTRUE;

      AddNodes(sensNodes, kZDC, pAddPsd, pVerbose);
      CalcDimensions(kZDC, Nx, Ny, Nz);

      if (pVerbose) cout << endl << "-W-  MpdZdcDigiScheme::Init: finished." << endl;
   }
   return kTRUE;
}

// -------------------------------------------------------------------------

MpdZdcVolInfo_t *MpdZdcDigiScheme::CreateVolInfoElement(FairGeoNode *nod, Int_t pVerbose)
{
   if (!nod) return nullptr;

   static TString  root_name_copy1 = "";
   static Double_t volData[6]      = {0, 0, 0, 0, 0, 0};
   Int_t           i; //,j;

   TString   shape_name = nod->getShapePointer()->GetName();
   TString   root_name  = nod->getRootVolume()->GetName();
   TString   tmp        = "          ";
   TGeoBBox *shape      = 0;

   FairGeoVector pos = nod->getLabTransform()->getTranslation();

   for (i = 0; i < 3; i++) volData[i] = pos.getValues(i); // [cm]

   if (root_name != root_name_copy1) {

      root_name_copy1 = root_name;

      shape = (TGeoBBox *)nod->getRootVolume()->GetShape();

      if (shape_name.Contains("BOX")) {
         volData[3] = shape->GetDX(); // [cm]
         volData[4] = shape->GetDY(); // [cm]
         volData[5] = shape->GetDZ(); // [cm]
      } else {
         if (shape_name.Contains("PGON")) {
            volData[3] = nod->getPoint(2)->getZ() / 10.; // [cm]
            volData[4] = volData[3];
            volData[5] = shape->GetDZ(); // [cm]
         }
      }
   }

   MpdZdcVolInfo_t *volInfo = new MpdZdcVolInfo_t;
   for (i = 0; i < 6; i++) volInfo->push_back(volData[i]);

   if (pVerbose > 1) {
      if (pVerbose > 2) tmp = "  root: " + root_name + ", shape: " + shape_name;
      tmp += "X,Y,Z,   Dx,Dy,Dz: ";
      cout << tmp << volData[0] << "," << volData[1] << "," << volData[2] << ",   " << volData[3] << "," << volData[4]
           << "," << volData[5] << endl;
   }
   return volInfo;
}

// -------------------------------------------------------------------------

MpdZdcVolId_t *MpdZdcDigiScheme::CreateVolElement(FairGeoNode *nod, Int_t nodeNumber, MpdZdcDigiId_t *right,
                                                  Int_t pGlobalDetectorNumber, Int_t pVerbose)
{
   if (!nod) return nullptr;

   FairGeoNode *nod0, *nod1;
   TString      mother_name, tmp;

   if (!fPasNodes) return nullptr;

   nod0 = (FairGeoNode *)fPasNodes->At(0);
   nod1 = (FairGeoNode *)fPasNodes->At(1);

   if ((!nod0) || (!nod1)) return nullptr;

   mother_name = nod->getMother();

   if (mother_name == nod1->GetName()) {

      (*right).push_back(nod0->getCopyNo() + pGlobalDetectorNumber - 1);
      (*right).push_back(nod1->getCopyNo());
      (*right).push_back(nodeNumber);

      MpdZdcVolId_t *left = new MpdZdcVolId_t;
      (*left).push_back(nod0->getCopyNo());
      (*left).push_back(nod1->getCopyNo());
      (*left).push_back(nod->getMCid());
      (*left).push_back(nod->getCopyNo());

      if (pVerbose > 1) {
         tmp = "  mother: " + mother_name + ", me: " + (nod->getName()) + "    ";
         cout << tmp << (*left)[0] << "," << (*left)[1] << "," << (*left)[2] << "," << (*left)[3] << " : "
              << (*right)[0] << "," << (*right)[1] << "," << (*right)[2] << endl;
      }

      return left;

   } else {
      cout << "-E-  MpdZdcDigiScheme::CreateVolInfoElement:  Strange for me node: " << nod->GetName()
           << "  Node number:" << nodeNumber << "  Mother:" << mother_name << endl;
      return nullptr;
   }
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::AddNodes(TObjArray *sensNodes, Int_t pGlobalDetectorNumber, Bool_t pAddPsd, Int_t pVerbose)
{
   Int_t           nNodes = sensNodes->GetEntriesFast();
   FairGeoNode *   nod    = 0;
   Int_t           nodeNumber; //,nodeCopyNo,nodeVolumeId, chanId2=0, chanId1=0;
   MpdZdcVolId_t * left1, *left2;
   MpdZdcDigiId_t *right1, *right2;

   if (pVerbose) {
      cout << "-W-  MpdZdcDigiScheme::AddNodes: started:"
           << "GlobalDetectorNumber:" << pGlobalDetectorNumber << "  nNodes:" << nNodes << endl;
   }

   for (nodeNumber = 0; nodeNumber < nNodes; nodeNumber++) {
      nod = (FairGeoNode *)sensNodes->At(nodeNumber);
      if (nod) {

         right1 = new MpdZdcDigiId_t;
         left1  = CreateVolElement(nod, nodeNumber, right1, pGlobalDetectorNumber, pVerbose);

         fVolToDigiIdMap[*left1] = *right1;

         left2        = new MpdZdcDigiId_t((*left1).begin(), (*left1).end());
         (*left2)[0]  = (*left1)[0] + 1;
         right2       = new MpdZdcDigiId_t(right1->begin(), right1->end());
         (*right2)[0] = (*right1)[0] + 1;

         fVolToDigiIdMap[*left2] = *right2;

         CreateVolCopyElements(left1, right1);

         if (pAddPsd) {
            MpdZdcVolInfo_t *volInfo1 = CreateVolInfoElement(nod, pVerbose);
            if (volInfo1) {
               fDigiToVolInfoMap[*right1] = volInfo1;
               MpdZdcVolInfo_t *volInfo2  = new MpdZdcVolInfo_t(*volInfo1);
               (*volInfo2)[2]             = -(*volInfo2)[2]; // Z
               fDigiToVolInfoMap[*right2] = volInfo2;

               CreateVolInfoCopyElements(right1, volInfo1);
            }
         }
      } else {
         cout << "-W-  MpdZdcDigiScheme::AddNodes: "
              << "Node number " << nodeNumber << " from " << nNodes << " not found!" << endl;
         return kFALSE;
      }
   }
   return kTRUE;
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::CreateVolCopyElements(MpdZdcVolId_t *left, MpdZdcDigiId_t *right)
{
   MpdZdcVolId_t * left1, *left2;
   MpdZdcDigiId_t *right1, *right2;

   if (!fPasNodes) return kFALSE;

   FairGeoNode *nod1 = (FairGeoNode *)fPasNodes->At(fPasNodes->GetEntries() - 2);

   if (!nod1) return kFALSE;

   Int_t moduleID, nModules; // {MotherMotherCopyNo, MotherCopyNo, VolumeId, CopyNo},{DetectorID, ModuleID, ChannelID}

   TString lm_name = nod1->getName();
   TString last_module_number(&(lm_name[lm_name.Last('#') + 1]));
   nModules = last_module_number.Atoi();

   for (moduleID = 1; moduleID < nModules; moduleID++) {

      left1        = new MpdZdcDigiId_t((*left).begin(), (*left).end());
      (*left1)[1]  = (*left)[1] + moduleID;
      right1       = new MpdZdcDigiId_t(right->begin(), right->end());
      (*right1)[1] = (*right)[1] + moduleID;

      left2        = new MpdZdcDigiId_t((*left1).begin(), (*left1).end());
      (*left2)[0]  = (*left1)[0] + 1;
      right2       = new MpdZdcDigiId_t(right1->begin(), right1->end());
      (*right2)[0] = (*right1)[0] + 1;

      fVolToDigiIdMap[*left1] = *right1;
      fVolToDigiIdMap[*left2] = *right2;
   }

   return kTRUE;
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::CreateVolInfoCopyElements(MpdZdcDigiId_t *right, MpdZdcVolInfo_t *volInfo)
{

   if (!fPasNodes) return kFALSE;

   FairGeoNode *nod1 = (FairGeoNode *)fPasNodes->At(fPasNodes->GetEntries() - 2);

   if (!nod1) return kFALSE;

   Int_t moduleID, nModules;

   TString lm_name = nod1->getName();
   TString last_module_number(&(lm_name[lm_name.Last('#') + 1]));
   nModules = last_module_number.Atoi();

   for (moduleID = 1; moduleID < nModules; moduleID++) {

      nod1 = (FairGeoNode *)fPasNodes->At(fPasNodes->GetEntries() - 2 - nModules + moduleID);

      FairGeoVector pos = nod1->getLabTransform()->getTranslation();

      MpdZdcVolInfo_t *volInfo1 = new MpdZdcVolInfo_t((*volInfo).begin(), (*volInfo).end());
      (*volInfo1)[0]            = pos.getValues(0); // X [cm]
      (*volInfo1)[1]            = pos.getValues(1); // Y [cm]

      MpdZdcDigiId_t *right1 = new MpdZdcDigiId_t(right->begin(), right->end());
      (*right1)[1]           = (*right)[1] + moduleID;

      MpdZdcDigiId_t *right2 = new MpdZdcDigiId_t(right1->begin(), right1->end());
      (*right2)[0]           = (*right1)[0] + 1;

      MpdZdcVolInfo_t *volInfo2 = new MpdZdcVolInfo_t(*volInfo1);
      (*volInfo2)[2]            = -(*volInfo2)[2]; // Z [cm]

      fDigiToVolInfoMap[*right1] = volInfo1;
      fDigiToVolInfoMap[*right2] = volInfo2;
   }

   return kTRUE;
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::InitPsd(MpdZdcPsdGeoPar *geoPar, MpdZdcPsdDigiPar *digiPar, Int_t pVerbose)
{
   if (!fYesPsd) {
      if (!fInitialized) {
         cout << "-W- MpdZdcDigiScheme::InitPsd: "
              << "Init basic ZDC first!" << endl;
         return kFALSE;
      }

      fYesPsd = kTRUE;

      if (!geoPar) {
         cout << "-W- MpdZdcDigiScheme::InitPsd: "
              << "No geometry parameters available!" << endl;
         return kFALSE;
      }
      if (!digiPar) {
         cout << "-W-  MpdZdcDigiScheme::InitPsd: "
              << "No digitization parameters available!" << endl;
         return kFALSE;
      }
      fZdcPsdDigiPar = digiPar;

      TObjArray *sensNodes;

      sensNodes = geoPar->GetGeoSensitiveNodes();
      if (!sensNodes) {
         cout << "-W-  MpdZdcDigiScheme::InitPsd: "
              << "No sensitive nodes available!" << endl;
         return kFALSE;
      }

      AddNodes(sensNodes, kZDC + 2, kTRUE, pVerbose); // kZDC+2 - hard-coded global ZDC_PSD detector number
      CalcDimensions(kZDC + 2, Nx_psd, Ny_psd, Nz_psd);
      fYesPsd = kTRUE;
   }
   return kTRUE;
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::GetVolCenterXYZ(MpdZdcDigiId_t *pDigiID, Double_t &x, Double_t &y, Double_t &z)
{
   if (!pDigiID) return kFALSE;
   if (fDigiToVolInfoMap.find(*pDigiID) == fDigiToVolInfoMap.end())
      return kFALSE;
   else {
      MpdZdcVolInfo_t *volInfo = fDigiToVolInfoMap[*pDigiID];
      if (volInfo) {
         x = (*volInfo)[0];
         y = (*volInfo)[1];
         z = (*volInfo)[2];
         return kTRUE;
      } else
         return kFALSE;
   }
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::GetVolDxDyDz(MpdZdcDigiId_t *pDigiID, Double_t &Dx, Double_t &Dy, Double_t &Dz)
{
   if (!pDigiID) return kFALSE;
   if (fDigiToVolInfoMap.find(*pDigiID) == fDigiToVolInfoMap.end())
      return kFALSE;
   else {
      MpdZdcVolInfo_t *volInfo = fDigiToVolInfoMap[*pDigiID];
      if (volInfo) {
         Dx = (*volInfo)[3];
         Dy = (*volInfo)[4];
         Dz = (*volInfo)[5];
         return kTRUE;
      } else
         return kFALSE;
   }
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::IsVolumeExist(MpdZdcVolId_t *pVolId)
{
   if (!pVolId)
      return kFALSE;
   else
      return (!(fVolToDigiIdMap.find(*pVolId) == fVolToDigiIdMap.end()));
}

// -------------------------------------------------------------------------

MpdZdcDigiId_t MpdZdcDigiScheme::GetDigiId(MpdZdcVolId_t *pVolId)
{
   static const MpdZdcDigiId_t not_found(4, -1);

   if (IsVolumeExist(pVolId))
      return fVolToDigiIdMap[*pVolId];
   else
      return not_found;
}

// -------------------------------------------------------------------------

Int_t MpdZdcDigiScheme::GetDetectorID(MpdZdcVolId_t *pVolId)
{
   MpdZdcDigiId_t digiID = GetDigiId(pVolId);
   //  return digiID.first;
   return digiID[0];
}

// -------------------------------------------------------------------------

Int_t MpdZdcDigiScheme::GetChannelID(MpdZdcVolId_t *pVolId)
{
   MpdZdcDigiId_t digiID = GetDigiId(pVolId);
   //  return digiID.second;
   return digiID[2];
}

// -------------------------------------------------------------------------

MpdZdcVolInfo_t *MpdZdcDigiScheme::GetVolInfo(MpdZdcVolId_t *pVolId)
{
   if (IsVolumeExist(pVolId)) {

      MpdZdcDigiId_t pDigiID = GetDigiId(pVolId);

      if (fDigiToVolInfoMap.find(pDigiID) == fDigiToVolInfoMap.end())
         return nullptr;
      else
         return fDigiToVolInfoMap[pDigiID];
   } else
      return nullptr;
}

// -------------------------------------------------------------------------

void MpdZdcDigiScheme::PrintVolume(Int_t volID, Int_t copyNo, Int_t copyNoMother, Int_t copyNoMotherMother)
{
   Int_t         content[] = {copyNoMotherMother, copyNoMother, volID, copyNo};
   MpdZdcVolId_t pVolId(content, content + sizeof(content) / sizeof(Int_t));

   MpdZdcDigiId_t pDigiID = GetDigiId(&pVolId);

   cout << " MpdZdc Volume: " << copyNoMotherMother << "," << copyNoMother << "," << volID << "," << copyNo
        << "   DigiID: " << pDigiID[0] << "," << pDigiID[1] << "," << pDigiID[2];

   MpdZdcVolInfo_t *pVolInfo = GetVolInfo(&pVolId);

   if (pVolInfo)
      cout << "    X,Y,Z [cm]: " << (*pVolInfo)[0] << "," << (*pVolInfo)[1] << "," << (*pVolInfo)[2]
           << "    Dx,Dy,Dz [cm]: " << (*pVolInfo)[3] << "," << (*pVolInfo)[4] << "," << (*pVolInfo)[5];

   cout << endl;
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::CalcDimensions(Int_t pGlobalDetectorNumber, Int_t &nx, Int_t &ny, Int_t &nz)
{

   if (fDigiToVolInfoMap.empty()) return kFALSE;

   Bool_t                                                result = kFALSE;
   std::map<MpdZdcDigiId_t, MpdZdcVolInfo_t *>::iterator it;
   std::map<Double_t, Int_t>                             xmap, ymap, zmap;
   Double_t                                              x, y, z;
   nx = ny = nz = 0;

   for (it = fDigiToVolInfoMap.begin(); it != fDigiToVolInfoMap.end(); ++it) {

      //    cout << ". ";

      //    if ((*it).first.first==pGlobalDetectorNumber) {
      if ((*it).first[0] == 1) {

         result = kTRUE;

         x = (*((*it).second))[0];
         if (xmap.count(x))
            xmap[x] = xmap[x] + 1;
         else
            xmap[x] = 1;

         y = (*((*it).second))[1];
         if (ymap.count(y))
            ymap[y] = ymap[y] + 1;
         else
            ymap[y] = 1;

         z = (*((*it).second))[2];
         if (zmap.count(z))
            zmap[z] = zmap[z] + 1;
         else
            zmap[z] = 1;
      }
   }

   nx = xmap.size();
   ny = ymap.size();
   nz = zmap.size();

   return result;
}

// -------------------------------------------------------------------------

void MpdZdcDigiScheme::GetZdcDimensions(Int_t &nx, Int_t &ny, Int_t &nz)
{
   nx = Nx;
   ny = Ny;
   nz = Nz;
}

// -------------------------------------------------------------------------

void MpdZdcDigiScheme::GetZdcPsdDimensions(Int_t &nx_psd, Int_t &ny_psd, Int_t &nz_psd)
{
   nx_psd = Nx_psd;
   ny_psd = Ny_psd;
   nz_psd = Nz_psd;
}

// -------------------------------------------------------------------------

MpdZdcDigiId_t MpdZdcDigiScheme::GetDigiIdFromCoords(Double_t x, Double_t y, Double_t z)
{
   Int_t         content[] = {-1, -1, -1, -1};
   MpdZdcVolId_t resultmc(content, content + sizeof(content) / sizeof(Int_t));

   MpdZdcDigiId_t result(-1, -1);
   if (gGeoManager) {
      TGeoNode *tgn = gGeoManager->FindNode(x, y, z);
      if (tgn) {
         resultmc[0] = 1 * (z > 0) + 2 * (z < 0);
         resultmc[1] = tgn->GetMotherVolume()->GetNumber();
         resultmc[2] = tgn->GetVolume()->GetNumber();
         resultmc[3] = tgn->GetNumber();
         if (fVolToDigiIdMap.find(resultmc) != fVolToDigiIdMap.end()) result = fVolToDigiIdMap[resultmc];
      }
   }
   return result;
}

// -------------------------------------------------------------------------

MpdZdcDigiId_t MpdZdcDigiScheme::GetDigiIdFromVolumeData(Int_t pMcVolumeNumber, Int_t pMcCopyNumber,
                                                         Int_t pMotherCopyNumber, Int_t pMotherMotherCopyNumber)
{
   Int_t          content[] = {pMotherMotherCopyNumber, pMotherCopyNumber, pMcVolumeNumber, pMcCopyNumber};
   MpdZdcVolId_t  pVolId(content, content + sizeof(content) / sizeof(Int_t));
   MpdZdcDigiId_t digiID = GetDigiId(&pVolId);
   return digiID;
}

// -------------------------------------------------------------------------

void MpdZdcDigiScheme::SplitDigiID(MpdZdcDigiId_t digiID, Int_t &detID, Int_t &modID, Int_t &chanID)
{
   detID  = digiID[0];
   modID  = digiID[1];
   chanID = digiID[2];
}

// -------------------------------------------------------------------------

void MpdZdcDigiScheme::Print()
{
   cout << "*********************************************" << endl;

   cout << "***  MpdZdcDigiScheme:" << endl;
   cout << " MpdZdc  Nx,Ny,Nz:  " << Nx << "," << Ny << "," << Nz << "       YesPsd: " << fYesPsd;

   if (fYesPsd) cout << "  MpdZdcPsd  Nx_psd,Ny_psd,Nz_psd:  " << Nx_psd << "," << Ny_psd << "," << Nz_psd;
   cout << endl;

   std::map<MpdZdcVolId_t, MpdZdcDigiId_t>::iterator it;

   for (it = fVolToDigiIdMap.begin(); it != fVolToDigiIdMap.end(); ++it)
      PrintVolume((*it).first[2], (*it).first[3], (*it).first[1], (*it).first[0]);

   cout << "*********************************************" << endl;
}

// -------------------------------------------------------------------------

Bool_t MpdZdcDigiScheme::GetDetIdModIdChanId(Int_t pMcVolumeNumber, Int_t pMcCopyNumber, Int_t pMotherCopyNumber,
                                             Int_t pMotherMotherCopyNumber, Int_t &pDetId, Int_t &pChanId,
                                             Int_t &pModId)
{
   Int_t          content[] = {pMotherMotherCopyNumber, pMotherCopyNumber, pMcVolumeNumber, pMcCopyNumber};
   MpdZdcVolId_t  pVolId(content, content + sizeof(content) / sizeof(Int_t));
   MpdZdcDigiId_t digiID = GetDigiId(&pVolId);
   pDetId                = digiID[0];
   pModId                = digiID[1];
   pChanId               = digiID[2];

   return kTRUE;
}

// -------------------------------------------------------------------------
ClassImp(MpdZdcDigiScheme);
