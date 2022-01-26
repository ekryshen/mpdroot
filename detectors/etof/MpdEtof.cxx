#include <iostream>
#include <algorithm>
#include <limits>

#include <TGeoManager.h>
#include <TClonesArray.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "MpdStack.h"
#include "FairRuntimeDb.h"
#include "FairRunAna.h"
#include "FairVolume.h"

#include "MpdEtofGeo.h"
#include "MpdTofPoint.h"
#include "MpdEtofGeoPar.h"

#include "MpdEtof.h"

ClassImp(MpdEtof) using namespace std;
//------------------------------------------------------------------------------------------------------------------------
MpdEtof::MpdEtof(const char *name, Bool_t active)
   : FairDetector(name, active), nan(std::numeric_limits<double>::quiet_NaN())
{
   aTofHits      = new TClonesArray("MpdTofPoint");
   fPosIndex     = 0;
   fVerboseLevel = 1;
}
//------------------------------------------------------------------------------------------------------------------------
MpdEtof::~MpdEtof()
{
   if (aTofHits) {
      aTofHits->Delete();
      delete aTofHits;
   }
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdEtof::ProcessHits(FairVolume *vol)
{
   Int_t strip, box, sector; // detector,

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

   // Create MpdTofPoint at exit of active volume
   if (fELoss > 0 && (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())) {
      fTrackID = gMC->GetStack()->GetCurrentTrackNumber();

      TString padname = vol->GetName();
      padname.Remove(0, 7);
      strip = padname.Atoi(); // gMC->CurrentVolOffID(1, pad);
      // detector = 1;
      gMC->CurrentVolOffID(2, box);
      gMC->CurrentVolOffID(3, sector);

      fVolumeID = MpdTofPoint::GetSuid72(sector, box, strip); // stub!!!

      AddPoint(fTrackID, fVolumeID, fPos.Vect(), fMom.Vect(), fTime, fLength, fELoss);

      ((MpdStack *)gMC->GetStack())->AddPoint(kETOF);

      ResetParameters();
   }

   return kTRUE;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtof::EndOfEvent()
{
   if (fVerboseLevel) Print();
   aTofHits->Delete();
   fPosIndex = 0;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtof::Register()
{
   FairRootManager::Instance()->Register("ETOFPoint", "Etof", aTofHits, kTRUE);
}
//------------------------------------------------------------------------------------------------------------------------
TClonesArray *MpdEtof::GetCollection(Int_t iColl) const
{
   if (iColl == 0) return aTofHits;

   return nullptr;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtof::Print() const
{
   Int_t nHits = aTofHits->GetEntriesFast();
   cout << "-I- MpdEtof: " << nHits << " points registered in this event." << endl;

   if (fVerboseLevel > 1)
      for (Int_t i = 0; i < nHits; i++) (*aTofHits)[i]->Print();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtof::Reset()
{
   aTofHits->Delete();
   ResetParameters();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtof::CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset)
{
   Int_t nEntries = cl1->GetEntriesFast();
   cout << "-I- MpdEtof: " << nEntries << " entries to add." << endl;
   TClonesArray &clref    = *cl2;
   MpdTofPoint * oldpoint = nullptr;

   for (Int_t i = 0; i < nEntries; i++) {
      oldpoint    = (MpdTofPoint *)cl1->At(i);
      Int_t index = oldpoint->GetTrackID() + offset;
      oldpoint->SetTrackID(index);
      new (clref[fPosIndex]) MpdTofPoint(*oldpoint);
      fPosIndex++;
   }

   cout << "-I- MpdEtof: " << cl2->GetEntriesFast() << " merged entries." << endl;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtof::ConstructGeometry()
{
   FairGeoLoader *   geoLoad = FairGeoLoader::Instance();
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   MpdEtofGeo *      tofGeo  = new MpdEtofGeo();
   tofGeo->setGeomFile(GetGeometryFileName());
   geoFace->addGeoModule(tofGeo);

   Bool_t rc = geoFace->readSet(tofGeo);
   if (rc) tofGeo->create(geoLoad->getGeoBuilder());
   TList *volList = tofGeo->getListOfVolumes();

   // store geo parameter
   FairRun *      fRun       = FairRun::Instance();
   FairRuntimeDb *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdEtofGeoPar *par        = (MpdEtofGeoPar *)(rtdb->getContainer("MpdEtofGeoPar"));
   TObjArray *    fSensNodes = par->GetGeoSensitiveNodes();
   TObjArray *    fPassNodes = par->GetGeoPassiveNodes();

   FairGeoNode *  node = nullptr;
   FairGeoVolume *aVol = nullptr;
   TListIter      iter(volList);

   while ((node = (FairGeoNode *)iter.Next())) {
      aVol = dynamic_cast<FairGeoVolume *>(node);
      if (node->isSensitive())
         fSensNodes->AddLast(aVol);
      else
         fPassNodes->AddLast(aVol);
   }

   par->setChanged();
   par->setInputVersion(fRun->GetRunId(), 1);
   ProcessNodes(volList);
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofPoint *MpdEtof::AddPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                               Double_t eLoss)
{
   return new ((*aTofHits)[aTofHits->GetEntriesFast()]) MpdTofPoint(trackID, detID, pos, mom, time, length, eLoss);
}
//------------------------------------------------------------------------------------------------------------------------