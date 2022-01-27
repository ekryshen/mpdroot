
//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTof, LRectangle, LStrip
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <algorithm>
#include <cmath>

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
#include "MpdMCTrack.h"

#include "MpdTofGeo.h"
#include "MpdTofGeoPar.h"
#include "MpdTofPoint.h"
#include "MpdTofHit.h"

#include "MpdTof.h"
using namespace std;

ClassImp(MpdTof);
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
MpdTof::MpdTof(const char *name, Bool_t active) : FairDetector(name, active)
{
   aTofHits      = new TClonesArray("MpdTofPoint");
   fVerboseLevel = 1;
}
//------------------------------------------------------------------------------------------------------------------------
MpdTof::~MpdTof()
{
   aTofHits->Delete();
   delete aTofHits;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::ResetParameters()
{
   fTrackID = fVolumeID = 0;
   fPos.SetXYZM(NAN, NAN, NAN, NAN);
   fMom  = fPos;
   fTime = fLength = fELoss = NAN;
   fPosIndex                = 0;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdTof::ProcessHits(FairVolume *vol)
{
   Int_t sector, detector, strip;

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

      gMC->CurrentVolOffID(0, strip);    // [1,72]
      gMC->CurrentVolOffID(1, detector); // [1,20]
      gMC->CurrentVolOffID(2, sector);   // [1,14]
      fVolumeID = MpdTofPoint::GetSuid72(sector, detector, strip);

      AddPoint(fTrackID, fVolumeID, fPos.Vect(), fMom.Vect(), fTime, fLength, fELoss);

      ((MpdStack *)gMC->GetStack())->AddPoint(kTOF);

      ResetParameters();
   }

   return kTRUE;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::EndOfEvent()
{
   if (fVerboseLevel) Print();

   aTofHits->Delete();
   fPosIndex = 0;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::Register()
{
   FairRootManager::Instance()->Register("TOFPoint", "Tof", aTofHits, kTRUE);
}
//------------------------------------------------------------------------------------------------------------------------
TClonesArray *MpdTof::GetCollection(Int_t iColl) const
{
   if (iColl == 0) return aTofHits;

   return nullptr;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::Print() const
{
   Int_t nHits = aTofHits->GetEntriesFast();
   cout << "-I- MpdTof: " << nHits << " points registered in this event." << endl;

   if (fVerboseLevel > 1)
      for (Int_t i = 0; i < nHits; i++) (*aTofHits)[i]->Print();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::Dump(TClonesArray *aPoints, TClonesArray *aHits, TClonesArray *aTracks, const char *comment, ostream &os)
{
   if (comment != nullptr) os << comment;

   // Fill & sort MpdTofPoint&MpdTofHit by mctid key.
   Int_t nHits   = aHits->GetEntriesFast();
   Int_t nPoints = aPoints->GetEntriesFast();
   Int_t nTracks = aTracks->GetEntriesFast();

   multimap<Int_t, Int_t> mmPoints, mmHits, mmTracks; // < mctid, poinIndex or hitIndex>  < mctid parent, mctid>
   set<Int_t>             sMcIndex;

   vector<Int_t> vec;
   for (Int_t index = 0; index < nHits; index++) {
      auto *pHit = (MpdTofHit *)aHits->At(index);
      vec.clear();
      pHit->getLinks(MpdTofUtils::mcTrackIndex, vec);

      for (auto it = vec.begin(), itEnd = vec.end(); it != itEnd; it++) {
         Int_t mctid = *it;
         mmHits.insert({mctid, index});
         sMcIndex.insert(mctid);
      }
   }
   for (Int_t index = 0; index < nPoints; index++) {
      auto *pPoint = (MpdTofPoint *)aPoints->At(index);
      Int_t mctid  = pPoint->GetTrackID();
      mmPoints.insert({mctid, index});
      sMcIndex.insert(mctid);
   }
   for (Int_t index = 0; index < nTracks; index++) {
      auto pMCtrack = (MpdMCTrack *)aTracks->UncheckedAt(index);
      auto ptid     = pMCtrack->GetMotherId();
      mmTracks.insert({ptid, index});
   }

   // Print to ostream os.
   os << "\n[MpdTof::Dump] points: " << nPoints << ", hits: " << nHits
      << " --------------------------------------------------------->>>";
   for (auto mctid : sMcIndex) {
      os << "\ntid=" << mctid << " points: (" << mmPoints.count(mctid) << ") ";

      auto range = mmPoints.equal_range(mctid);
      for (auto it = range.first; it != range.second; ++it) os << " pid=" << it->second;

      os << " hits: (" << mmHits.count(mctid) << ") ";

      auto range2 = mmHits.equal_range(mctid);
      for (auto it = range2.first; it != range2.second; ++it) {
         auto  hid  = it->second;
         auto  pHit = (MpdTofHit *)aHits->At(hid);
         Int_t suid = pHit->GetDetectorID();

         Int_t sector, detector, gap, strip;
         MpdTofPoint::ParseSuid(suid, sector, detector, gap, strip);

         os << " hid=" << hid << "[" << suid << "]{" << sector << "," << detector << "," << gap << "," << strip << "}";
      }
   }
   os << "\nTrack tree "
         "---------------------------------------------------------------------------------------------------------";
   for (Int_t index = 0; index < nTracks; index++) {
      if (mmTracks.find(index) != mmTracks.end()) continue; // pass only ended track line

      os << "\n tid=" << index;

      auto pMCtrack = (MpdMCTrack *)aTracks->At(index);
      do {
         Int_t ptid = pMCtrack->GetMotherId();
         if (ptid != -1) {
            os << "<tid=" << ptid;
            pMCtrack = (MpdMCTrack *)aTracks->At(ptid);
         } else
            pMCtrack = nullptr;
      } while (pMCtrack);
   }
   os << "\n[MpdTof::Dump] "
         "--------------------------------------------------------------------------------------------------<<<";
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::Reset()
{
   aTofHits->Delete();
   ResetParameters();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset)
{
   Int_t nEntries = cl1->GetEntriesFast();
   cout << "-I- MpdTof: " << nEntries << " entries to add." << endl;
   TClonesArray &clref    = *cl2;
   MpdTofPoint  *oldpoint = nullptr;

   for (Int_t i = 0; i < nEntries; i++) {
      oldpoint    = (MpdTofPoint *)cl1->At(i);
      Int_t index = oldpoint->GetTrackID() + offset;
      oldpoint->SetTrackID(index);
      new (clref[fPosIndex]) MpdTofPoint(*oldpoint);
      fPosIndex++;
   }

   cout << "-I- MpdTof: " << cl2->GetEntriesFast() << " merged entries." << endl;
}
//--------------------------------------------------------------------------------------------------------------------------------------
void MpdTof::ConstructGeometry()
{
   TString fileName = GetGeometryFileName();
   if (fileName.EndsWith(".root")) {
      LOG(INFO) << "Constructing TOF geometry from ROOT file " << fileName.Data();
      ConstructRootGeometry();
   } else if (fileName.EndsWith(".geo")) {
      LOG(INFO) << "Constructing TOF geometry from ASCII file " << fileName.Data();
      ConstructAsciiGeometry();
   } else
      LOG(FATAL) << "Geometry format of TOF file " << fileName.Data() << " not supported.";
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::ConstructAsciiGeometry()
{
   FairGeoLoader    *geoLoad = FairGeoLoader::Instance();
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   MpdTofGeo        *tofGeo  = new MpdTofGeo();
   tofGeo->setGeomFile(GetGeometryFileName());
   geoFace->addGeoModule(tofGeo);

   Bool_t rc = geoFace->readSet(tofGeo);
   if (rc) tofGeo->create(geoLoad->getGeoBuilder());
   TList *volList = tofGeo->getListOfVolumes();

   // store geo parameter

   FairRuntimeDb *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdTofGeoPar  *par        = (MpdTofGeoPar *)(rtdb->getContainer("MpdTofGeoPar"));
   TObjArray     *fSensNodes = par->GetGeoSensitiveNodes();
   TObjArray     *fPassNodes = par->GetGeoPassiveNodes();

   FairGeoNode   *node = nullptr;
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
   par->setInputVersion(FairRun::Instance()->GetRunId(), 1);
   ProcessNodes(volList);
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofPoint *MpdTof::AddPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                              Double_t eLoss)
{
   return new ((*aTofHits)[aTofHits->GetEntriesFast()]) MpdTofPoint(trackID, detID, pos, mom, time, length, eLoss);
}
//--------------------------------------------------------------------------------------------------------------------------------------
Bool_t MpdTof::CheckIfSensitive(string name)
{
   if (name.find("Active") != string::npos) return kTRUE;

   return kFALSE;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::Print(const TVector3 &v, const char *comment, ostream &os)
{
   if (comment != nullptr) os << comment;
   os << "(" << v.X() << "," << v.Y() << "," << v.Z() << "; " << v.Perp() << "," << v.Mag() << ")";
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTof::GetDelta(const TVector3 &mcPos, const TVector3 &estPos, double &dev, double &devZ, double &devR,
                      double &devPhi)
{
   dev    = (mcPos - estPos).Mag();
   devZ   = mcPos.Z() - estPos.Z();
   devR   = mcPos.Perp() - estPos.Perp();
   devPhi = sqrt(dev * dev - devZ * devZ - devR * devR);
   devPhi = (mcPos.Phi() > estPos.Phi()) ? devPhi : -1. * devPhi;
}
//------------------------------------------------------------------------------------------------------------------------
