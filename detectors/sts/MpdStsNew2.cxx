// -------------------------------------------------------------------------
// -----                       MpdSts source file                      -----
// -----                        18.04.2012   VK                        -----
// -------------------------------------------------------------------------
//  Independant strip readout for each sector
//  Geometry: its_sec_cables.geo ( Three types of sectors consist of one,
//            two or three sensors)

#include "MpdStsNew2.h"
#include "MpdStsGeo.h"
#include "MpdStsGeoPar.h"
#include "MpdStsHit.h"
#include "MpdStsPoint.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "MpdStack.h"
#include "FairRuntimeDb.h"
#include "TObjArray.h"
#include "FairRun.h"
#include "FairVolume.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TVirtualMC.h"

#include <iostream>

//--------------------------------------------------------------------------
MpdStsNew2::MpdStsNew2() : FairDetector("STS", kTRUE)
{
   fStsCollection = new TClonesArray("MpdStsPoint");
   fPosIndex      = 0;
   fVerboseLevel  = 1;
}
//--------------------------------------------------------------------------
MpdStsNew2::MpdStsNew2(const char *name, Bool_t active) : FairDetector(name, active)
{
   fStsCollection = new TClonesArray("MpdStsPoint");
   fPosIndex      = 0;
   fVerboseLevel  = 1;
}
//--------------------------------------------------------------------------
MpdStsNew2::~MpdStsNew2()
{
   if (fStsCollection) {
      fStsCollection->Delete();
      delete fStsCollection;
   }
}
//--------------------------------------------------------------------------
Bool_t MpdStsNew2::ProcessHits(FairVolume *vol)
{
   Int_t   sector, sectorType, layer, ladder, side = 0;
   TString Volname;

   // Set parameters at entrance of volume. Reset ELoss.
   if (gMC->IsTrackEntering()) {
      fELoss  = 0.;
      fTime   = gMC->TrackTime() * 1.0e09;
      fLength = gMC->TrackLength();
      gMC->TrackPosition(fPos);
      gMC->TrackMomentum(fMom);
      // cout << gMC->GetStack()->GetCurrentTrackNumber() << " " << vol->getRealName() << endl;
   }

   // Sum energy loss for all steps in the active volume
   try {
      fELoss += gMC->Edep();
   } catch (...) {
      cout << "-E- MpdStsNew2::ProcessHits gMC->Edep() exception   !!!" << ((UInt_t *)gMC) << endl;
   }

   // Create MpdStsPoint at exit of active volume
   if ((gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) && fELoss > 0) {

      // Exiting volume or stopped track -> save point
      fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
      Volname  = vol->getRealName(); // EL
      // cout << fTrackID << " " << Volname << " " << gMC->CurrentVolPath() << endl;

      if (Volname.Contains("sector")) {

         // Sector geometry
         layer      = Volname[12] - 48; // Number of Layer
         sectorType = Volname[11] - 48; // Type of sector

         gMC->CurrentVolOffID(0, sector); // Copy Nr of sector
         gMC->CurrentVolOffID(1, ladder); // Copy Nr of ladder

         gMC->CurrentVolPath(); // path from cave to sts01sector

         fVolumeID = side; // # of side      on bit 0
         // fVolumeID |= (sector << 1);      // # of sector    on bit 1 - 5
         // fVolumeID |= (ladder << 6);      // # of ladder    on bit 6 - 10
         // fVolumeID |= (layer << 11);      // # of layer     on bit 11 - 13
         // fVolumeID |= (sectorType << 14); // type of sector (1,2 or 3)    on bit 14 - 16
         fVolumeID |= (sector << MpdStsHit::kSensorS);      // # of sector    on bit 1 - 5
         fVolumeID |= (ladder << MpdStsHit::kLadderS);      // # of ladder    on bit 6 - 10
         fVolumeID |= (layer << MpdStsHit::kLayerS);        // # of layer     on bit 11 - 13
         fVolumeID |= (sectorType << MpdStsHit::kSecTypeS); // type of sector (1,2 or 3)    on bit 14 - 16

         //**************************************************************************************
         /*
          cout << fTrackID << " " << Volname << " " << gMC->CurrentVolPath() << endl;
          cout << layer << " " << ladder  << " " << sectorType << " " << sector << " " << side << " " <<endl;
          cout << fVolumeID << endl;
          cout << ((fVolumeID >> 10) & 7) << " " << ((fVolumeID >> 5) & 31) << " " << ((fVolumeID >> 13) & 3)
          << " " << ((fVolumeID >> 1) & 15) << " " << (fVolumeID & 1) << endl;
          */
         //*************************************************************************************
      }

      TVector3 vposIn = TVector3(fPos.X(), fPos.Y(), fPos.Z());
      TVector3 vmomIn = TVector3(fMom.Px(), fMom.Py(), fMom.Pz());
      gMC->TrackPosition(fPos);
      TVector3 vposOut = TVector3(fPos.X(), fPos.Y(), fPos.Z());

      AddHit(fTrackID, fVolumeID, vposIn, vmomIn, vposOut, fTime, fLength, fELoss);

      ((MpdStack *)gMC->GetStack())->AddPoint(kSTS);

      ResetParameters();
   }
   return kTRUE;
}

//-------------------------------------------------------------------------------------------------------
void MpdStsNew2::EndOfEvent()
{
   if (fVerboseLevel) Print();
   fStsCollection->Delete();
   fPosIndex = 0;
}
//------------------------------------------------------------------------------------------------------
void MpdStsNew2::Register()
{
   FairRootManager::Instance()->Register("StsPoint", "Sts", fStsCollection, kTRUE);
}
//-------------------------------------------------------------------------------------------------------
TClonesArray *MpdStsNew2::GetCollection(Int_t iColl) const
{
   if (iColl == 0) return fStsCollection;

   return NULL;
}
//-------------------------------------------------------------------------------------------------------
void MpdStsNew2::Print() const
{
   Int_t nHits = fStsCollection->GetEntriesFast();
   cout << "-I- MpdStsNew2: " << nHits << " points registered in this event." << endl;

   if (fVerboseLevel > 1)
      for (Int_t i = 0; i < nHits; i++) (*fStsCollection)[i]->Print();
}
//-------------------------------------------------------------------------------------------------------
void MpdStsNew2::Reset()
{
   fStsCollection->Delete();
   ResetParameters();
}
//-------------------------------------------------------------------------------------------------------
void MpdStsNew2::CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset)
{
   Int_t nEntries = cl1->GetEntriesFast();
   cout << "-I- MpdStsNew2: " << nEntries << " entries to add." << endl;
   TClonesArray &clref    = *cl2;
   MpdStsPoint * oldpoint = NULL;

   for (Int_t i = 0; i < nEntries; i++) {
      oldpoint    = (MpdStsPoint *)cl1->At(i);
      Int_t index = oldpoint->GetTrackID() + offset;
      oldpoint->SetTrackID(index);
      new (clref[fPosIndex]) MpdStsPoint(*oldpoint);
      fPosIndex++;
   }

   cout << "-I- MpdStsNew2: " << cl2->GetEntriesFast() << " merged entries." << endl;
}
//-------------------------------------------------------------------------------------------------------
void MpdStsNew2::ConstructGeometry()
{
   Int_t count     = 0;
   Int_t count_tot = 0;

   FairGeoLoader *   geoLoad = FairGeoLoader::Instance();
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   MpdStsGeo *       stsGeo  = new MpdStsGeo();
   stsGeo->setGeomFile(GetGeometryFileName());
   geoFace->addGeoModule(stsGeo);

   Bool_t rc = geoFace->readSet(stsGeo);
   if (rc) stsGeo->create(geoLoad->getGeoBuilder());
   TList *volList = stsGeo->getListOfVolumes();

   // store geo parameter
   FairRun *      fRun       = FairRun::Instance();
   FairRuntimeDb *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdStsGeoPar * par        = (MpdStsGeoPar *)(rtdb->getContainer("MpdStsGeoPar"));
   TObjArray *    fSensNodes = par->GetGeoSensitiveNodes();
   TObjArray *    fPassNodes = par->GetGeoPassiveNodes();

   FairGeoNode *  node = NULL;
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
//--------------------------------------------------------------------------------------------------------
MpdStsPoint *MpdStsNew2::AddHit(Int_t trackID, Int_t detID, TVector3 posIn, TVector3 momIn, TVector3 posOut,
                                Double_t time, Double_t length, Double_t eLoss)
{
   TClonesArray &clref = *fStsCollection;
   Int_t         size  = clref.GetEntriesFast();

   return new (clref[size]) MpdStsPoint(trackID, detID, posIn, momIn, posOut, time, length, eLoss);
}
//--------------------------------------------------------------------------------------------------------

ClassImp(MpdStsNew2);
