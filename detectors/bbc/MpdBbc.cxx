//------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------
// -----                       MpdBbc source file                      -----
// -------------------------------------------------------------------------

#include <iostream>

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TVirtualMC.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairStack.h"
#include "FairRuntimeDb.h"
#include "TObjArray.h"
#include "FairRun.h"
#include "FairVolume.h"

#include "MpdBbc.h"
#include "MpdBbcGeo.h"
#include "MpdBbcPoint.h"
#include "MpdBbcGeoPar.h"

class FairVolume;

//------------------------------------------------------------------------------------------------------------------------
MpdBbc::MpdBbc() : FairDetector("BBC", kTRUE)
{
   fBbcCollection = new TClonesArray("MpdBbcPoint");
   fPosIndex      = 0;
   fVerboseLevel  = 1;
}
//------------------------------------------------------------------------------------------------------------------------
MpdBbc::MpdBbc(const char *name, Bool_t active) : FairDetector(name, active)
{
   fBbcCollection = new TClonesArray("MpdBbcPoint");
   fPosIndex      = 0;
   fVerboseLevel  = 1;
}
//------------------------------------------------------------------------------------------------------------------------
MpdBbc::~MpdBbc()
{
   if (fBbcCollection) {
      fBbcCollection->Delete();
      delete fBbcCollection;
   }
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdBbc::ProcessHits(FairVolume *vol)
{
   Int_t   gap, cell, module, region;
   TString Volname;

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

   // Create MpdBbcPoint at exit of active volume
   //	if(gMC->TrackCharge() != 0 &&  (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) )

   // Create MpdBbcPoint at ENTER of active volume; fELoss INVALID!!!
   if (gMC->IsTrackEntering()) {
      fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
      Volname  = vol->getRealName();
      region   = Volname[5] - '0'; //?????????????????????????
      gMC->CurrentVolID(gap);
      gMC->CurrentVolOffID(1, cell);
      gMC->CurrentVolOffID(2, module);

      fVolumeID = ((region - 1) << 24); ////////////// + ((module-1)<<14) + ((cell-1)<<4) + (gap-1);

      AddHit(fTrackID, fVolumeID, TVector3(fPos.X(), fPos.Y(), fPos.Z()), TVector3(fMom.Px(), fMom.Py(), fMom.Pz()),
             fTime, fLength, fELoss);

      ((FairStack *)gMC->GetStack())->AddPoint(kBBC);

      ResetParameters();
   }

   return kTRUE;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBbc::EndOfEvent()
{
   if (fVerboseLevel) Print();
   fBbcCollection->Clear();
   fPosIndex = 0;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBbc::Register()
{
   FairRootManager::Instance()->Register("BBCPoint", "Bbc", fBbcCollection, kTRUE);
}
//------------------------------------------------------------------------------------------------------------------------
TClonesArray *MpdBbc::GetCollection(Int_t iColl) const
{
   if (iColl == 0) return fBbcCollection;

   return NULL;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBbc::Print() const
{
   Int_t nHits = fBbcCollection->GetEntriesFast();
   cout << "-I- MpdBbc: " << nHits << " points registered in this event." << endl;

   if (fVerboseLevel > 1)
      for (Int_t i = 0; i < nHits; i++) (*fBbcCollection)[i]->Print();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBbc::Reset()
{
   fBbcCollection->Clear();
   ResetParameters();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBbc::CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset)
{
   Int_t nEntries = cl1->GetEntriesFast();
   cout << "-I- MpdBbc: " << nEntries << " entries to add." << endl;
   TClonesArray &clref    = *cl2;
   MpdBbcPoint  *oldpoint = NULL;

   for (Int_t i = 0; i < nEntries; i++) {
      oldpoint    = (MpdBbcPoint *)cl1->At(i);
      Int_t index = oldpoint->GetTrackID() + offset;
      oldpoint->SetTrackID(index);
      new (clref[fPosIndex]) MpdBbcPoint(*oldpoint);
      fPosIndex++;
   }

   cout << "-I- MpdBbc: " << cl2->GetEntriesFast() << " merged entries." << endl;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBbc::ConstructGeometry()
{

   Int_t count     = 0;
   Int_t count_tot = 0;

   FairGeoLoader    *geoLoad = FairGeoLoader::Instance();
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   MpdBbcGeo        *bbcGeo  = new MpdBbcGeo();
   bbcGeo->setGeomFile(GetGeometryFileName());
   geoFace->addGeoModule(bbcGeo);

   Bool_t rc = geoFace->readSet(bbcGeo);
   if (rc) bbcGeo->create(geoLoad->getGeoBuilder());
   TList *volList = bbcGeo->getListOfVolumes();

   // store geo parameter
   FairRun       *fRun       = FairRun::Instance();
   FairRuntimeDb *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdBbcGeoPar  *par        = (MpdBbcGeoPar *)(rtdb->getContainer("MpdBbcGeoPar"));
   TObjArray     *fSensNodes = par->GetGeoSensitiveNodes();
   TObjArray     *fPassNodes = par->GetGeoPassiveNodes();

   FairGeoNode   *node = NULL;
   FairGeoVolume *aVol = NULL;
   TListIter      iter(volList);

   while ((node = (FairGeoNode *)iter.Next())) {
      aVol = dynamic_cast<FairGeoVolume *>(node);
      if (node->isSensitive()) {
         fSensNodes->AddLast(aVol);
         count++;
      } else
         fPassNodes->AddLast(aVol);
      count_tot++;
   }

   par->setChanged();
   par->setInputVersion(fRun->GetRunId(), 1);
   ProcessNodes(volList);
}
//------------------------------------------------------------------------------------------------------------------------
MpdBbcPoint *MpdBbc::AddHit(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                            Double_t eLoss)
{
   TClonesArray &clref = *fBbcCollection;
   Int_t         size  = clref.GetEntriesFast();

   return new (clref[size]) MpdBbcPoint(trackID, detID, pos, mom, time, length, eLoss);
}
//------------------------------------------------------------------------------------------------------------------------

ClassImp(MpdBbc)
