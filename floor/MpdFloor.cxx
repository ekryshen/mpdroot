/*
 * MpdFloor.cxx
 *
 *  Created on: 20 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdFloor.h"

#include "FairGeoInterface.h"  // for FairGeoInterface
#include "FairGeoLoader.h"     // for FairGeoLoader
#include "FairGeoNode.h"       // for FairGeoNode
#include "FairGeoVolume.h"     // for FairGeoVolume
#include "FairLogger.h"        // for logging
#include "FairRootManager.h"   // for FairRootManager
#include "FairRun.h"           // for FairRun
#include "FairRuntimeDb.h"     // for FairRuntimeDb
#include "FairVolume.h"        // for FairVolume
#include "MpdDetectorList.h"   // for DetectorId::kTutDet
#include "MpdFloorGeo.h"       // for FairTutorialDet1Geo
#include "MpdFloorGeo.h"
#include "MpdFloorGeoPar.h"  // for FairTutorialDet1GeoPar
#include "MpdFloorGeoPar.h"
#include "MpdFloorPoint.h"
#include "MpdFloorPoint.h"  // for FairTutorialDet1Point
#include "MpdFloorPoint.h"
#include "MpdStack.h"
#include "TLorentzVector.h"
#include "TVirtualMC.h"

MpdFloorGeo* MpdFloor::fgGeo = nullptr;


MpdFloor::MpdFloor() : MpdFloor("Floor", kTRUE) {}

MpdFloor::MpdFloor(const char* name, Bool_t active) :
  FairDetector(name, active), fTrackID(0), fVolumeID(0), fPos(0, 0, 0, 0), fMom(0, 0, 0, 0), fTime(0), fLength(0), fELoss(0) {
  fPointCollection = new TClonesArray("MpdFloorPoint");
}

MpdFloor::~MpdFloor() { delete fPointCollection; }

Bool_t MpdFloor::ProcessHits(FairVolume* vol) {
  LOG(debug) << "In Floor ::ProcessHits";
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
  if (TVirtualMC::GetMC()->IsTrackExiting() || TVirtualMC::GetMC()->IsTrackStop() || TVirtualMC::GetMC()->IsTrackDisappeared()) {
    fTrackID  = TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();
    fVolumeID = vol->getMCid();
    // if (fELoss == 0.) { return kFALSE; }
    AddHit(fTrackID, fVolumeID, fPos.Vect(), fMom.Vect(), fTime, fLength, fELoss);
  }

  return kTRUE;
}

void MpdFloor::EndOfEvent() { fPointCollection->Delete(); }

void MpdFloor::Register() {
  if (!gMC->IsMT()) {
    FairRootManager::Instance()->Register("FloorPoint", "Floor", fPointCollection, kTRUE);
  } else {
    FairRootManager::Instance()->RegisterAny("FloorPoint", fPointCollection, kTRUE);
  }
}

TClonesArray* MpdFloor::GetCollection(Int_t iColl) const {
  if (iColl == 0) {
    return fPointCollection;
  } else {
    return nullptr;
  }
}

void MpdFloor::Print() const {
  Int_t nHits = fPointCollection->GetEntriesFast();
  std::cout << "-I- MpdTof: " << nHits << " points registered in this event." << std::endl;
}

void MpdFloor::Reset() { fPointCollection->Clear(); }

void MpdFloor::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset) {}

void MpdFloor::ConstructGeometry() {
  TString fileName = GetGeometryFileName();

  if (fileName.EndsWith(".root")) {
    LOG(info) << "Constructing Floor geometry from ROOT file" << fileName.Data();
    ConstructRootGeometry();
    return;
  }

  LOG(info) << "Constructing Floor geometry from ASCII file  " << fileName;

  FairGeoLoader* loader          = FairGeoLoader::Instance();
  FairGeoInterface* GeoInterface = loader->getGeoInterface();
  MpdFloorGeo* MGeo              = new MpdFloorGeo();
  MGeo->setGeomFile(GetGeometryFileName());
  GeoInterface->addGeoModule(MGeo);
  Bool_t rc = GeoInterface->readSet(MGeo);
  if (rc) { MGeo->create(loader->getGeoBuilder()); }

  TList* volList = MGeo->getListOfVolumes();
  // store geo parameter
  // FairRun* fRun = FairRun::Instance();
  FairRuntimeDb* rtdb   = FairRun::Instance()->GetRuntimeDb();
  MpdFloorGeoPar* par   = static_cast<MpdFloorGeoPar*>(rtdb->getContainer("MpdFloorGeoPar"));
  TObjArray* fSensNodes = par->GetGeoSensitiveNodes();
  TObjArray* fPassNodes = par->GetGeoPassiveNodes();

  TListIter iter(volList);

  FairGeoNode* node   = nullptr;
  FairGeoVolume* aVol = nullptr;

  while ((node = (FairGeoNode*) iter.Next())) {
    aVol = dynamic_cast<FairGeoVolume*>(node);
    if (node->isSensitive())
      fSensNodes->AddLast(aVol);
    else
      fPassNodes->AddLast(aVol);
  }

  par->setChanged();
  par->setInputVersion(FairRun::Instance()->GetRunId(), 1);
  ProcessNodes(volList);
}

Bool_t MpdFloor::CheckIfSensitive(std::string name) {
  // all detector active
  TString tsname = name;
  if (tsname.BeginsWith("fl")) {
    return kTRUE;
  } else {
    return kFALSE;
  }
}

MpdFloorPoint*
MpdFloor::AddHit(Int_t trackID, Int_t detId, TVector3 pos, TVector3 mom, Double_t time, Double_t length, Double_t eloss) {
  TClonesArray& clref = *fPointCollection;
  Int_t size          = clref.GetEntriesFast();
  MpdStack* stack     = static_cast<MpdStack*>(TVirtualMC::GetMC()->GetStack());
  stack->AddPoint(kFLOOR);
  return new (clref[size]) MpdFloorPoint(trackID, detId, pos, mom, time, length, eloss);
}
