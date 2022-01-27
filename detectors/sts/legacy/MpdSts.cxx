//--------------------------------------------------------------------------
// -------------------------------------------------------------------------
// -----                       MpdSts source file                      -----
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
#include "MpdStack.h"
#include "FairRuntimeDb.h"
#include "TObjArray.h"
#include "FairRun.h"
#include "FairVolume.h"

#include "MpdSts.h"
#include "MpdStsGeo.h"
#include "MpdStsPoint.h"
#include "MpdStsGeoPar.h"

class FairVolume;

//--------------------------------------------------------------------------
MpdSts::MpdSts() : FairDetector("STS", kTRUE)
{
   fStsCollection = new TClonesArray("MpdStsPoint");
   fPosIndex      = 0;
   fVerboseLevel  = 1;
}
//--------------------------------------------------------------------------
MpdSts::MpdSts(const char *name, Bool_t active) : FairDetector(name, active)
{
   fStsCollection = new TClonesArray("MpdStsPoint");
   fPosIndex      = 0;
   fVerboseLevel  = 1;
}
//--------------------------------------------------------------------------
MpdSts::~MpdSts()
{
   if (fStsCollection) {
      fStsCollection->Delete();
      delete fStsCollection;
   }
}
//--------------------------------------------------------------------------
Bool_t MpdSts::ProcessHits(FairVolume *vol)
{
   Int_t   sensor, detector, layer, ladder, side = 0;
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
   try {
      fELoss += gMC->Edep();
   } catch (...) {
      cout << "-E- MpdSts::ProcessHits gMC->Edep() exception   !!!" << ((UInt_t *)gMC) << endl;
   }

   // Create MpdStsPoint at exit of active volume
   // if (gMC->TrackCharge() != 0 &&  (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) )
   if ((gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) && fELoss > 0) {
      // Exiting volume or stopped track - save point

      fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
      Volname  = vol->getRealName(); // EL

      if (Volname.Contains("sensor")) {
         // Modular geometry
         layer = Volname[11] - 48; // Layer

         gMC->CurrentVolID(sensor);         // CopyNo of sensor
         gMC->CurrentVolOffID(1, detector); // CopyNo of detector
         gMC->CurrentVolOffID(2, ladder);   // CopyNo of ladder

         gMC->CurrentVolOffName(2); // name of ladder
         gMC->CurrentVolOffName(1); // name of detector
         gMC->CurrentVolOffName(0); // name of sensor
         gMC->CurrentVolPath();     // path from cave to sts01sensor

         fVolumeID = side; // side     on bit 0 - 1   because number  of sides=1 (0,1)
         fVolumeID |=
            (layer << 1); // layer    on bit 1 - 4   because number of layers=4, 4(dec)->4(bin)=100 => num of bytes = 3
         fVolumeID |= (ladder << 5); // ladder   on bit 5 - 9   because max number of ladders=18, 18(dec)->18(bin)=10010
                                     // => num of bytes = 5
         fVolumeID |= (detector << 10); // detector on bit 10 - 14  because max number of detectors=24,
                                        // 24(dec)->24(bin)=11000 => num of bytes = 5

         //**************************************************************************************

         // cout << fTrackID << " " << Volname << " " << gMC->CurrentVolOffName(1) << " "
         //	 << gMC->CurrentVolOffName(0) << " " << gMC->CurrentVolPath() << " " << cell << " "
         //	 << gMC->CurrentVolOffID(1, ladder) << " " << endl;
         /*
         cout << fTrackID << " " << Volname << " " << gMC->CurrentVolPath() << endl;
         cout << layer << " " << ladder  << " " << detector << " " << side << " " <<endl;
         cout << fVolumeID << endl;
         cout << (fVolumeID >> 10) << " " << ((fVolumeID >> 5) & 31) << " " << ((fVolumeID >> 1) & 15) << " " <<
         (fVolumeID & 1) << endl;
         */
         //*************************************************************************************
      } else {
         // Simple "cylinder" geometry
         Int_t region = Volname[5] - 96;
         fVolumeID    = (region << 1);
      }

      TVector3 vposIn = TVector3(fPos.X(), fPos.Y(), fPos.Z());
      TVector3 vmomIn = TVector3(fMom.Px(), fMom.Py(), fMom.Pz());
      gMC->TrackPosition(fPos);
      TVector3 vposOut = TVector3(fPos.X(), fPos.Y(), fPos.Z());

      //		AddHit(fTrackID, fVolumeID, TVector3(fPos.X(),  fPos.Y(),  fPos.Z()),
      //	   		TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime, fLength, fELoss);
      AddHit(fTrackID, fVolumeID, vposIn, vmomIn, vposOut, fTime, fLength, fELoss);

      ((MpdStack *)gMC->GetStack())->AddPoint(kSTS);

      ResetParameters();
   }
   return kTRUE;
}

//-------------------------------------------------------------------------------------------------------
void MpdSts::EndOfEvent()
{
   if (fVerboseLevel) Print();
   fStsCollection->Delete();
   fPosIndex = 0;
}
//------------------------------------------------------------------------------------------------------
void MpdSts::Register()
{
   FairRootManager::Instance()->Register("StsPoint", "Sts", fStsCollection, kTRUE);
}
//-------------------------------------------------------------------------------------------------------
TClonesArray *MpdSts::GetCollection(Int_t iColl) const
{
   if (iColl == 0) return fStsCollection;

   return NULL;
}
//-------------------------------------------------------------------------------------------------------
void MpdSts::Print() const
{
   Int_t nHits = fStsCollection->GetEntriesFast();
   cout << "-I- MpdSts: " << nHits << " points registered in this event." << endl;

   if (fVerboseLevel > 1)
      for (Int_t i = 0; i < nHits; i++) (*fStsCollection)[i]->Print();
}
//-------------------------------------------------------------------------------------------------------
void MpdSts::Reset()
{
   fStsCollection->Delete();
   ResetParameters();
}
//-------------------------------------------------------------------------------------------------------
void MpdSts::CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset)
{
   Int_t nEntries = cl1->GetEntriesFast();
   cout << "-I- MpdSts: " << nEntries << " entries to add." << endl;
   TClonesArray &clref    = *cl2;
   MpdStsPoint  *oldpoint = NULL;

   for (Int_t i = 0; i < nEntries; i++) {
      oldpoint    = (MpdStsPoint *)cl1->At(i);
      Int_t index = oldpoint->GetTrackID() + offset;
      oldpoint->SetTrackID(index);
      new (clref[fPosIndex]) MpdStsPoint(*oldpoint);
      fPosIndex++;
   }

   cout << "-I- MpdSts: " << cl2->GetEntriesFast() << " merged entries." << endl;
}
//-------------------------------------------------------------------------------------------------------
void MpdSts::ConstructGeometry()
{

   Int_t count     = 0;
   Int_t count_tot = 0;

   FairGeoLoader    *geoLoad = FairGeoLoader::Instance();
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   MpdStsGeo        *stsGeo  = new MpdStsGeo();
   stsGeo->setGeomFile(GetGeometryFileName());
   geoFace->addGeoModule(stsGeo);

   Bool_t rc = geoFace->readSet(stsGeo);
   if (rc) stsGeo->create(geoLoad->getGeoBuilder());
   TList *volList = stsGeo->getListOfVolumes();

   // store geo parameter
   FairRun       *fRun       = FairRun::Instance();
   FairRuntimeDb *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdStsGeoPar  *par        = (MpdStsGeoPar *)(rtdb->getContainer("MpdStsGeoPar"));
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
//--------------------------------------------------------------------------------------------------------
MpdStsPoint *MpdSts::AddHit(Int_t trackID, Int_t detID, TVector3 posIn, TVector3 momIn, TVector3 posOut, Double_t time,
                            Double_t length, Double_t eLoss)
{
   TClonesArray &clref = *fStsCollection;
   Int_t         size  = clref.GetEntriesFast();

   return new (clref[size]) MpdStsPoint(trackID, detID, posIn, momIn, posOut, time, length, eLoss);
}
//--------------------------------------------------------------------------------------------------------

ClassImp(MpdSts);
