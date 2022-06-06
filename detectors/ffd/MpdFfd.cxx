//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFfd
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include "iostream"

#include "MpdMCTrack.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TVirtualMC.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "MpdStack.h"
#include "TObjArray.h"

#include "MpdFfd.h"
#include "MpdFfdGeo.h"
#include "MpdFfdGeoPar.h"
#include "MpdFfdPoint.h"
#include "MpdFfdHitProducer.h"

class FairVolume;
using namespace std;

ClassImp(MpdFfd);
//------------------------------------------------------------------------------------------------------------------------
MpdFfd::MpdFfd() : FairDetector("FFD", kTRUE)
{
   aFfdPoints    = new TClonesArray("MpdFfdPoint");
   fPosIndex     = 0;
   fVerboseLevel = 1;
}
//------------------------------------------------------------------------------------------------------------------------
MpdFfd::MpdFfd(const char *name, Bool_t active, Int_t verbose) : FairDetector(name, active)
{
   aFfdPoints    = new TClonesArray("MpdFfdPoint");
   fPosIndex     = 0;
   fVerboseLevel = verbose;
}
//------------------------------------------------------------------------------------------------------------------------
MpdFfd::~MpdFfd()
{
   aFfdPoints->Delete();
   delete aFfdPoints;
}
//------------------------------------------------------------------------------------------------------------------------
MpdFfdPoint *MpdFfd::FindPoint(Int_t ptid, Int_t suid)
{
   if (cPtid == ptid && cSuid == suid) return cPoint; // cached result

   // looking for new pair <tid, suid>
   for (int i = 0, N = aFfdPoints->GetEntriesFast(); i < N; i++) // cycle already created MpdFfdPoint
   {
      auto entry = (MpdFfdPoint *)aFfdPoints->UncheckedAt(i);

      if (entry->IsSame(ptid, suid)) {
         cPoint = entry; // update current point
         cPtid  = ptid;
         cSuid  = suid;
         return cPoint;
      }
   }
   return nullptr;
}
//------------------------------------------------------------------------------------------------------------------------
MpdFfdPoint *MpdFfd::CreatePoint(Int_t tid, Int_t suid)
{
   return new ((*aFfdPoints)[aFfdPoints->GetEntriesFast()]) MpdFfdPoint(tid, suid);
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdFfd::ProcessHits(FairVolume *vol)
{
   Int_t suid;
   gMC->CurrentVolID(suid);
   Int_t tid = gMC->GetStack()->GetCurrentTrackNumber();
   Int_t pid = gMC->TrackPid(); // Cherenkov 50000050  FeedbackPhoton 50000051

   if (gMC->IsTrackEntering()) {
      if (pid != 50000050) //  event = non-op track ENTER to quartz volume
      {
         gMC->TrackPosition(fPos);
         gMC->TrackMomentum(fMom);

         auto ret = FindTrackParams(tid, suid);
         if (ret.second == false) // don't found
         {
            TrackParam param(suid, fPos.Vect(), fMom, gMC->TrackTime() * 1.e+9, gMC->TrackLength());
            params.insert(make_pair(tid, param)); // insert new parent track parameters
         }
      } else // event = op create at quartz volume
      {
         TLorentzVector mom;
         gMC->TrackMomentum(mom);

         double energy = mom.E() * 1.e+9; // op energy [eV]

         bool addEntry = true;
         if (MpdFfdPoint::fCurrentMode == MpdFfdPoint::kPhotoElectron) // current mode = save pe only
         {
            if (!MpdFfdHitProducer::IsPeCreated(energy)) addEntry = false; // failed create pe from op
         }

         if (addEntry) {
            Int_t ptid = gMC->GetStack()->GetCurrentParentTrackNumber();

            auto point = FindPoint(ptid, suid);
            if (point == nullptr) point = CreatePoint(ptid, suid);

            if (!point->AddOp(energy,
                              gMC->TrackTime() *
                                 1.e+9)) // [ns]; if MpdFfdPoint not closed, need to call SaveParentTrackParams ONCE!
            {
               auto ret = FindTrackParams(ptid, suid);
               if (ret.second) {
                  auto param = ret.first->second;
                  point->SaveParentTrackParams(param.posIn, param.posOut, param.mom, param.time, param.length);

                  ((MpdStack *)gMC->GetStack())->AddPoint(kFFD, ptid); // add marker for parent track
               }
            }

         } // addEntry

      } // op

   } // track entering

   // Update MpdFfdPoint at non-op track EXIT from quartz
   if (pid != 50000050 && (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())) {
      TLorentzVector posOut;
      gMC->TrackPosition(posOut);

      auto ret = FindTrackParams(tid, suid);
      if (ret.second) ret.first->second.posOut = posOut.Vect(); // update parent track output position

      ResetParameters();

   } // track exiting

   return kTRUE;
}
//------------------------------------------------------------------------------------------------------------------------
pair<MpdFfd::Tparams::iterator, bool> MpdFfd::FindTrackParams(Int_t tid, Int_t suid)
{
   auto ret = params.equal_range(tid);
   for (auto it = ret.first; it != ret.second; ++it) {
      if (it->second.suid == suid) {
         return make_pair(it, true);
      }
   }

   return make_pair(params.end(), false);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfd::EndOfEvent()
{
   if (fVerboseLevel) Print();

   aFfdPoints->Delete();
   fPosIndex = 0;

   // event data cleanup
   params.clear();
   cPoint = nullptr;
   cPtid = cSuid = -1;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfd::Register()
{
   FairRootManager::Instance()->Register("FFDPoint", "Ffd", aFfdPoints, kTRUE);
}
//------------------------------------------------------------------------------------------------------------------------
TClonesArray *MpdFfd::GetCollection(Int_t iColl) const
{
   if (iColl == 0) return aFfdPoints;

   return nullptr;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfd::Print() const
{
   Int_t nHits = aFfdPoints->GetEntriesFast();
   cout << "-I- MpdFfd: " << nHits << " points registered in this event.\n";

   if (fVerboseLevel > 1)
      for (Int_t i = 0; i < nHits; i++) (*aFfdPoints)[i]->Print();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfd::Reset()
{
   aFfdPoints->Delete();
   ResetParameters();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfd::CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset)
{
   Int_t nEntries = cl1->GetEntriesFast();
   cout << "-I- MpdFfd: " << nEntries << " entries to add." << endl;
   TClonesArray &clref    = *cl2;
   MpdFfdPoint  *oldpoint = nullptr;

   for (Int_t i = 0; i < nEntries; i++) {
      oldpoint    = (MpdFfdPoint *)cl1->At(i);
      Int_t index = oldpoint->GetTrackID() + offset;
      oldpoint->SetTrackID(index);
      new (clref[fPosIndex]) MpdFfdPoint(*oldpoint);
      fPosIndex++;
   }

   cout << "-I- MpdFfd: " << cl2->GetEntriesFast() << " merged entries." << endl;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfd::ConstructGeometry()
{
   TString fileName = GetGeometryFileName();

   if (fileName.EndsWith(".root")) {
      LOG(info) << "Constructing FFD geometry from ROOT file " << fileName.Data() << endl;
      ConstructRootGeometry();
   } else if (fileName.EndsWith(".geo")) {
      LOG(info) << "Constructing FFD geometry from ASCII file " << fileName.Data() << endl;
      ConstructAsciiGeometry();
   }
   /*else if ( fileName.EndsWith(".gdml") )
   {
           LOG(info) << "Constructing CPC geometry from GDML file " <<
   fileName.Data() << endl; ConstructGDMLGeometry();
   }*/
   else {
      LOG(fatal) << "Geometry format of FFD file " << fileName.Data() << " not supported." << endl;
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfd::ConstructAsciiGeometry()
{

   int count     = 0;
   int count_tot = 0;

   FairGeoLoader    *geoLoad = FairGeoLoader::Instance();
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   MpdFfdGeo        *Geo     = new MpdFfdGeo();
   Geo->setGeomFile(GetGeometryFileName());
   geoFace->addGeoModule(Geo);

   Bool_t rc = geoFace->readSet(Geo);
   if (rc)
      Geo->create(geoLoad->getGeoBuilder());
   else
      std::cerr << "FfdDetector:: geometry could not be read!" << std::endl;

   TList *volList = Geo->getListOfVolumes();

   // store geo parameter
   FairRun       *fRun       = FairRun::Instance();
   FairRuntimeDb *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdFfdGeoPar  *par        = (MpdFfdGeoPar *)(rtdb->getContainer("MpdFfdGeoPar"));
   TObjArray     *fSensNodes = par->GetGeoSensitiveNodes();
   TObjArray     *fPassNodes = par->GetGeoPassiveNodes();

   TListIter      iter(volList);
   FairGeoNode   *node = nullptr;
   FairGeoVolume *aVol = nullptr;

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
Bool_t MpdFfd::CheckIfSensitive(std::string name)
{
   if (name.find("Active") != string::npos) return kTRUE;

   return kFALSE;
}
//------------------------------------------------------------------------------------------------------------------------
