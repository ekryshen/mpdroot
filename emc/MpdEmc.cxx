#include "MpdEmc.h"

#include "MpdEmcPoint.h"
#include "MpdEmcGeo.h"
#include "MpdEmcGeoPar.h"

#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "MpdDetectorList.h"
#include "MpdStack.h"

#include "TClonesArray.h"
#include "TVirtualMC.h"

#include <iostream>
using std::cout;
using std::endl;

//
// To start ECal MC production (add to geometry_stage1.C)
//
//  FairDetector *Emc = new MpdEmc("EMC", kTRUE);
//  Emc->SetGeometryFileName("emc_v4.root");
//  fRun->AddModule(Emc);
//

MpdEmc::MpdEmc() : FairDetector("EMC", kTRUE, kECAL)
{
   /** create your collection for data points */
   fMpdEmcPointCollection = new TClonesArray("MpdEmcPoint");
}

MpdEmc::MpdEmc(const char *name, Bool_t active) : FairDetector(name, active, kECAL)
{
   fMpdEmcPointCollection = new TClonesArray("MpdEmcPoint");
   fVerboseLevel          = 1;
}

MpdEmc::~MpdEmc()
{
   if (fMpdEmcPointCollection) {
      fMpdEmcPointCollection->Delete();
      delete fMpdEmcPointCollection;
   }
}

/*
void MpdEmc::Initialize()
{
  FairDetector::Initialize();
  FairRuntimeDb *rtdb= FairRun::Instance()->GetRuntimeDb();
  //MpdEmcGeoPar* par=(MpdEmcGeoPar*)(rtdb->getContainer("MpdEmcGeoPar"));
}
*/

Bool_t MpdEmc::ProcessHits(FairVolume *vol)
{

   /** This method is called from the MC stepping */

   //#define EDEBUG
#ifdef EDEBUG
   static int my_counter = 0;
   if (my_counter < 10) {
      cout << "-I- MpdEmc::ProcessHits IsTrackEntering=" << (gMC->IsTrackEntering()) << endl;
      my_counter++;
   }
#endif
#undef EDEBUG

   // Set parameters at entrance of volume. Reset ELoss.
   if (gMC->IsTrackEntering()) {
      fELoss  = 0.;
      fTime   = gMC->TrackTime() * 1.0e09;
      fLength = gMC->TrackLength();
      gMC->TrackPosition(fPos);
      gMC->TrackMomentum(fMom);
   }

   // Sum energy loss for all steps in the active volume

   fELoss += gMC->Edep();

   // Create MpdEmcPoint at exit of active volume

   // Process Hit program

   if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) {
      fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
      ///////    fVolumeID = vol->getMCid();

      fVolumeID = GetDetectorEmcID();
      if (fELoss == 0.) return kFALSE;

      AddHit(fTrackID, fVolumeID, TVector3(fPos.X(), fPos.Y(), fPos.Z()), TVector3(fMom.Px(), fMom.Py(), fMom.Pz()),
             fTime, fLength, fELoss);

      // Increment number of tutorial det points in TParticle
      //    printf("Name x, y, z %s %f %f %f \n", vol->GetName(), fPos.X(), fPos.Y(), fPos.Z());

      MpdStack *stack = (MpdStack *)gMC->GetStack();
      stack->AddPoint(kECAL);
   }

   return kTRUE;
}

void MpdEmc::EndOfEvent()
{

   Print();

   fMpdEmcPointCollection->Clear();
}

Int_t MpdEmc::GetDetectorEmcID()
{

   Int_t       towerNumber = -1, towerXYNum, towerZNum;
   const Int_t totalTowerZ = 128, totalTowerXY = 300, totalTowerSectorXY = 12;
   Int_t       chamberNum, sectorNum, crateNum, moduleNum, towerType, towerNum;
   if (sscanf(gMC->CurrentVolPath(),
              "/cave_1/emcChamber_0/emcChH_%d/emcSector_%d/emcCrate_%d/emcModule%d_0/emc_box%d_%d/", &chamberNum,
              &sectorNum, &crateNum, &moduleNum, &towerType, &towerNum) > 0) {
      towerXYNum = sectorNum * totalTowerSectorXY + 2 * crateNum + towerNum % 2;
      towerZNum  = 0.5 * totalTowerZ + towerType - 1;
      if (chamberNum == 1) towerZNum = 0.5 * totalTowerZ - towerType;
      towerNumber = towerXYNum * 1000 + towerZNum;
   }

   return towerNumber;
}

void MpdEmc::Register()
{

   /** This will create a branch in the output tree called
       MpdEmcPoint, setting the last parameter to kFALSE means:
       this collection will not be written to the file, it will exist
       only during the simulation.
   */

   FairRootManager::Instance()->Register("EmcPoint", "MpdEmc", fMpdEmcPointCollection, kTRUE);
}

TClonesArray *MpdEmc::GetCollection(Int_t iColl) const
{
   if (iColl == 0)
      return fMpdEmcPointCollection;
   else
      return NULL;
}

void MpdEmc::Reset()
{
   fMpdEmcPointCollection->Clear();
}

void MpdEmc::Print() const
{
   Int_t nHits = fMpdEmcPointCollection->GetEntriesFast();
   cout << "-I- MpdEmc: " << nHits << " points registered in this event." << endl;

   if (fVerboseLevel > 1)
      for (Int_t i = 0; i < nHits; i++) (*fMpdEmcPointCollection)[i]->Print();
}

void MpdEmc::ConstructGeometry()
{

   TString fileName = GetGeometryFileName();
   if (fileName.EndsWith(".root")) {
      LOG(INFO) << "Constructing EMC geometry from ROOT file " << fileName.Data();
      ConstructRootGeometry();
   } else if (fileName.EndsWith(".geo")) {
      LOG(INFO) << "Constructing EMC geometry from ASCII file " << fileName.Data();
      ConstructAsciiGeometry();
   } else {
      LOG(FATAL) << "Geometry format of EMC file " << fileName.Data() << " not supported.";
   }
}

void MpdEmc::ConstructAsciiGeometry()
{

   FairGeoLoader    *geoLoad = FairGeoLoader::Instance();
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   MpdEmcGeo        *Geo     = new MpdEmcGeo();
   Geo->setGeomFile(GetGeometryFileName());
   geoFace->addGeoModule(Geo);

   Bool_t rc = geoFace->readSet(Geo);
   if (rc) Geo->create(geoLoad->getGeoBuilder());
   TList *volList = Geo->getListOfVolumes();

   // store geo parameter
   FairRun       *fRun       = FairRun::Instance();
   FairRuntimeDb *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdEmcGeoPar  *par        = (MpdEmcGeoPar *)(rtdb->getContainer("MpdEmcGeoPar"));
   TObjArray     *fSensNodes = par->GetGeoSensitiveNodes();
   TObjArray     *fPassNodes = par->GetGeoPassiveNodes();

   TListIter      iter(volList);
   FairGeoNode   *node = NULL;
   FairGeoVolume *aVol = NULL;

   while ((node = (FairGeoNode *)iter.Next())) {
      aVol = dynamic_cast<FairGeoVolume *>(node);
      if (node->isSensitive()) {
         fSensNodes->AddLast(aVol);
      } else {
         fPassNodes->AddLast(aVol);
      }
   }
   par->setChanged();
   par->setInputVersion(fRun->GetRunId(), 1);

   ProcessNodes(volList);
}

// Check sensitivity

Bool_t MpdEmc::CheckIfSensitive(std::string name)
{
   TString tsname = name;
   if (tsname.Contains("cl_sc")) {
      return kTRUE;
   }
   return kFALSE;
}

MpdEmcPoint *MpdEmc::AddHit(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                            Double_t ELoss)
{
   TClonesArray &clref = *fMpdEmcPointCollection;
   Int_t         size  = clref.GetEntriesFast();
   return new (clref[size]) MpdEmcPoint(trackID, detID, pos, mom, time, length, ELoss);
}

ClassImp(MpdEmc)
