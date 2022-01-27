/*
 * MpdMcord.cxx
 *
 *  Created on: 20 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdMcord.h"

#include "MpdMcordPoint.h"
#include "MpdDetectorList.h"  // for DetectorId::kTutDet
#include "FairGeoInterface.h" // for FairGeoInterface
#include "FairGeoLoader.h"    // for FairGeoLoader
#include "FairGeoNode.h"      // for FairGeoNode
#include "FairGeoVolume.h"    // for FairGeoVolume
#include "FairRootManager.h"  // for FairRootManager
#include "FairRun.h"          // for FairRun
#include "FairRuntimeDb.h"    // for FairRuntimeDb
#include "MpdMcordGeo.h"      // for FairTutorialDet1Geo
#include "MpdMcordGeoPar.h"   // for FairTutorialDet1GeoPar
#include "MpdMcordPoint.h"    // for FairTutorialDet1Point
#include "FairVolume.h"       // for FairVolume
#include "FairLogger.h"       // for logging
#include "TVirtualMC.h"
#include "MpdStack.h"
#include "MpdMcordGeo.h"
#include "MpdMcordGeoPar.h"
#include "MpdMcordPoint.h"
#include "TLorentzVector.h"

MpdMcordGeo *MpdMcord::fgGeo = nullptr;

MpdMcord::MpdMcord() : MpdMcord("mcord", kTRUE) {}

MpdMcord::MpdMcord(const char *name, Bool_t active)
   : FairDetector(name, active), fTrackID(0), fVolumeID(0), fPos(0, 0, 0, 0), fMom(0, 0, 0, 0), fTime(0), fLength(0),
     fELoss(0)
{
   fPointCollection = new TClonesArray("MpdMcordPoint");
}

MpdMcord::~MpdMcord()
{
   delete fPointCollection;
}

Bool_t MpdMcord::ProcessHits(FairVolume *vol)
{
   LOG(debug) << "In Mcord ::ProcessHits";
   // Set parameters at entrance of volume. Reset ELoss.
   if (TVirtualMC::GetMC()->IsTrackEntering()) {
      fELoss  = 0.;
      fTime   = TVirtualMC::GetMC()->TrackTime() * 1.0e09;
      fLength = TVirtualMC::GetMC()->TrackLength();
      TVirtualMC::GetMC()->TrackPosition(fPos);
      TVirtualMC::GetMC()->TrackMomentum(fMom);
   }
   // Sum energy loss for all steps in the active volume
   fELoss += TVirtualMC::GetMC()->Edep();

   // Create FairTutorialDet1Point at exit of active volume
   if (TVirtualMC::GetMC()->IsTrackExiting() || TVirtualMC::GetMC()->IsTrackStop() ||
       TVirtualMC::GetMC()->IsTrackDisappeared()) {
      fTrackID  = TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();
      fVolumeID = vol->getMCid();
      if (fELoss == 0.) {
         return kFALSE;
      }
      AddHit(fTrackID, fVolumeID, fPos.Vect(), fMom.Vect(), fTime, fLength, fELoss);
   }

   return kTRUE;
}

void MpdMcord::EndOfEvent()
{
   fPointCollection->Delete();
}

void MpdMcord::Register()
{
   if (!gMC->IsMT()) {
      FairRootManager::Instance()->Register("McordPoint", "Mcord", fPointCollection, kTRUE);
   } else {
      FairRootManager::Instance()->RegisterAny("McordPoint", fPointCollection, kTRUE);
   }
}

TClonesArray *MpdMcord::GetCollection(Int_t iColl) const
{
   if (iColl == 0) {
      return fPointCollection;
   } else {
      return nullptr;
   }
}

void MpdMcord::Print() const
{
   Int_t nHits = fPointCollection->GetEntriesFast();
   std::cout << "-I- MpdTof: " << nHits << " points registered in this event." << std::endl;
}

void MpdMcord::Reset()
{
   fPointCollection->Clear();
}

void MpdMcord::CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset) {}

void MpdMcord::ConstructGeometry()
{
   TString fileName = GetGeometryFileName();

   if (fileName.EndsWith(".root")) {
      LOG(info) << "Constructing MCORD geometry from ROOT file" << fileName.Data();
      ConstructRootGeometry();
      return;
   }

   LOG(info) << "Constructing MCORD geometry from ASCII file  " << fileName;

   FairGeoLoader    *loader       = FairGeoLoader::Instance();
   FairGeoInterface *GeoInterface = loader->getGeoInterface();
   MpdMcordGeo      *MGeo         = new MpdMcordGeo();
   MGeo->setGeomFile(GetGeometryFileName());
   GeoInterface->addGeoModule(MGeo);
   Bool_t rc = GeoInterface->readSet(MGeo);
   if (rc) {
      MGeo->create(loader->getGeoBuilder());
   }

   TList *volList = MGeo->getListOfVolumes();
   // store geo parameter
   // FairRun* fRun = FairRun::Instance();
   FairRuntimeDb  *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdMcordGeoPar *par        = static_cast<MpdMcordGeoPar *>(rtdb->getContainer("MpdMcordGeoPar"));
   TObjArray      *fSensNodes = par->GetGeoSensitiveNodes();
   TObjArray      *fPassNodes = par->GetGeoPassiveNodes();

   TListIter iter(volList);

   FairGeoNode   *node = nullptr;
   FairGeoVolume *aVol = nullptr;

   while ((node = (FairGeoNode *)iter.Next())) {
      aVol = dynamic_cast<FairGeoVolume *>(node);
      if (node->isSensitive())
         fSensNodes->AddLast(aVol);
      else
         fPassNodes->AddLast(aVol);
   }

   par->setChanged();
   par->setInputVersion(FairRun::Instance()->GetRunId(), 1);
   ProcessNodes(volList);
}

Bool_t MpdMcord::CheckIfSensitive(std::string name)
{
   TString tsname = name;
   if (tsname.BeginsWith("md01scintVol")) {
      return kTRUE;
   } else {
      return kFALSE;
   }
}

MpdMcordPoint *MpdMcord::AddHit(Int_t trackID, Int_t detId, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                                Double_t eloss)
{
   TClonesArray &clref = *fPointCollection;
   Int_t         size  = clref.GetEntriesFast();
   MpdStack     *stack = static_cast<MpdStack *>(TVirtualMC::GetMC()->GetStack());
   stack->AddPoint(kMCORD);
   return new (clref[size]) MpdMcordPoint(trackID, detId, pos, mom, time, length, eloss);
}
