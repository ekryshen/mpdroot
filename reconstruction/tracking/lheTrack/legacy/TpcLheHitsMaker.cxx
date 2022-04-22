#include "TpcLheHitsMaker.h"

#include "TpcLheHit.h"
#include "TpcLheTrack.h"

#include "FairRunAna.h"
#include "TpcPoint.h"
#include "MpdMCTrack.h"
#include "FairMCApplication.h"
#include "FairRootManager.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "FairVolume.h"
#include "FairTask.h"
#include "TParticle.h"
#include "TRandom.h"
#include "Riostream.h"
#include "TH1F.h"

#include <iomanip>

using namespace std;

//#define TST
//#define TPC
#define TRACKS
#define GEANT
#define MASSPROD

//_________________________________________________________________
TpcLheHitsMaker::TpcLheHitsMaker()
{
   //---

   fLheHits     = new TClonesArray("TpcLheHit");
   fGeantTracks = new TClonesArray("TpcLheTrack");
}

//_________________________________________________________________
TpcLheHitsMaker::TpcLheHitsMaker(const char *name, const char *title) : FairTask(name)
{
   //---
   fLheHits     = new TClonesArray("TpcLheHit");
   fGeantTracks = new TClonesArray("TpcLheTrack");
}

//_________________________________________________________________
TpcLheHitsMaker::~TpcLheHitsMaker()
{
   //
   FairRootManager *fManager = FairRootManager::Instance();
   fManager->Write();

   delete fGeantTracks;
}

//_________________________________________________________________
InitStatus TpcLheHitsMaker::ReInit()
{
   // --- Update the pointers after reinitialisation

   //  cout << "TpcLheHitsMaker::ReInit() \n";

#if 0
  FairRunAna* ana = FairRunAna::Instance();
  FairRuntimeDb* rtdb=ana->GetRuntimeDb();
  fGeoPar=(FairGeoStsPar*)(rtdb->getContainer("FairGeoStsPar"));
#endif

   return kSUCCESS;
}

//_________________________________________________________________
InitStatus TpcLheHitsMaker::Init()
{

   //  cout << "InitStatus TpcLheHitsMaker::Init\n\n";

   fCuts = TpcLheTrackCuts::Instance();

   FairRootManager *fManager = FairRootManager::Instance();

   fTpcPoints = (TClonesArray *)fManager->GetObject("TpcPoint");
   if (!fTpcPoints) {
      cout << "-I- " << GetName() << "::Init: No TpcPoint array!" << endl;
      return kERROR;
   }

#if 0
  fTstPoints = (TClonesArray*) fManager->GetObject("TSTPoint");
  if ( ! fTstPoints ) {
    cout << "-I- "<< GetName() << "::Init: No TSTPoint array!" << endl;
    return kERROR;
  }
#endif

#ifdef GEANT
   fListMCtracks = (TClonesArray *)fManager->GetObject("MCTrack");
#endif

   fNofEvents = 0;

   Register();

   return kSUCCESS;
}

//_________________________________________________________________
void TpcLheHitsMaker::GetTstHits()
{
   //---

#ifdef TST
   cout << " TpcLheHitsMaker::GetTstHits() : TST Hit entries " << fTstPoints->GetEntriesFast() << endl;
#endif

   TpcLheHit *hit = NULL;

#if 0
  for (int j=0; j < fTstPoints->GetEntries(); j++ ) {
    FairTstPoint *point = (FairTstPoint*) fTstPoints->At(j);

#ifdef KTST
    cout << " hit: " << j << " detID " <<  point->GetStationNr() <<
    " event: " << info->GetEventNumber() << " fake: " <<
      info->IsFake() << "\n ";
#endif

    hit = AddHit();
    hit->SetHitNumber(fNHit++);

    hit->SetX(point->GetX());
    hit->SetY(point->GetY());
    hit->SetZ(point->GetZ());

    hit->SetTrackID(point->GetTrackID());

    hit->SetRefIndex(j);

#if 0
    hit->SetXerr(point->GetDx());
    hit->SetYerr(point->GetDy());
    hit->SetZerr(point->GetDz());

    hit->SetStation(point->GetStationNr()); //

    if (info->GetEventNumber() == 0 ) {
      //&& !info->IsFake()) {
      //FairMCPoint *mcpoint = (FairMCPoint *) fTpcPoints->At(point->GetRefIndex());
      //	hit->SetTrackID(mcpoint->GetTrackID());

	hit->SetTrackID(info->GetTrackId());
	hit->SetRefIndex(point->GetRefIndex());
      }
      else {
	hit->SetTrackID(-1);
      }
#endif

#ifdef TST
    hit->Print();
#endif

  }
#endif
}

//_________________________________________________________________
void TpcLheHitsMaker::GetTpcHits()
{

#ifdef TPC
   cout << " TpcLheHitsMaker::GetTpcHits(): Tpc points entries " << fTpcPoints->GetEntriesFast() << endl;
#endif

   Int_t nTpcPoints = fTpcPoints->GetEntriesFast();
   for (int j = 0; j < nTpcPoints; ++j) {

      TpcPoint *point = (TpcPoint *)fTpcPoints->UncheckedAt(j);

#ifdef TPC
      point->Print(" ");
      PR(point->GetTrackID());
#endif

      TpcLheHit *hit = AddHit();
      hit->SetHitNumber(fNHit++);

      hit->SetX(point->GetX());
      hit->SetY(point->GetY());
      hit->SetZ(point->GetZ());

      hit->SetXerr(0.03);
      hit->SetYerr(0.03);
      hit->SetZerr(0.03);

      hit->SetTrackID(point->GetTrackID());
      hit->SetRefIndex(j);

      // AZ
      Double_t dRPhi = 0, dZ = 0;
      // gRandom->Rannor(dRPhi,dZ); //AZ
      hit->SetZerr(0.1);                            // 1 mm error in Z
      hit->SetZ(hit->GetZ() + dZ * hit->GetZerr()); // add error
      TVector3 pos;
      point->Position(pos);
      hit->SetR(pos.Pt());
      hit->SetRphi(pos.Phi() * hit->GetR() + dRPhi * hit->GetRphiErr()); // add error
      Double_t phi = hit->GetRphi() / hit->GetR();
      hit->SetX(hit->GetR() * TMath::Cos(phi));
      hit->SetY(hit->GetR() * TMath::Sin(phi));
      hit->SetEdep(point->GetEnergyLoss());
      hit->SetStep(point->GetStep());
      // AZ

#ifdef TPC
      hit->Print();
#endif

   } // for (int j=0; j <  fStsHits->GetEntries(); j++ )
}

//_________________________________________________________________
void TpcLheHitsMaker::Exec(Option_t *option)
{

   cout << "\n\n  ***  Event # " << ++fNofEvents << endl;

   cout << " =====   TpcLheHitsMaker   =====\n";

   Reset();

   //  GetTstHits();

   GetTpcHits();

   cout << "  Total number of hits for tracking: " << setw(5) << fLheHits->GetEntries() << endl;

#ifdef GEANT
   CheckTracks();
#ifdef TRACK
   PrintTracks(0);
#endif
#endif

#ifdef MASSPROD
   //  fTpcPoints->Clear("C");
   //  fTstPoints->Clear("C");
#endif
}

//_________________________________________________________________
TpcLheHit *TpcLheHitsMaker::AddHit()
{
   // Creates a new hit in the TClonesArray.

   TClonesArray &hitRef = *fLheHits;
   Int_t         size   = hitRef.GetEntriesFast();
   return new (hitRef[size]) TpcLheHit();
}

//_________________________________________________________________
void TpcLheHitsMaker::Register()
{
   //---

   FairRootManager::Instance()->Register("LheGeantTrack", "Lhe", fGeantTracks, kFALSE);
   //    Register("LheGeantTrack","Lhe", fGeantTracks, kTRUE);

   FairRootManager::Instance()->Register("LheHit", "Lhe", fLheHits, kFALSE);
   // Register("LheHit","Lhe", fLheHits, kTRUE);
}

//_________________________________________________________________
void TpcLheHitsMaker::Finish()
{

   //   FairRootManager *fManager =FairRootManager::Instance();
   //   fManager->Fill();

   cout << " TpcLheHitsMaker::Finish : Geant tracks " << fGeantTracks->GetEntriesFast() << endl;

   fGeantTracks->Delete();
   fLheHits->Delete();
}

//_________________________________________________________________
void TpcLheHitsMaker::Reset()
{
   //---
   //  cout << " TpcLheHitsMaker::Reset  " << endl;

   fNHit   = 0;
   fNTrack = 0;

   if (fLheHits->GetEntries() != 0) fLheHits->Clear("C");
   if (fGeantTracks->GetEntries() != 0) fGeantTracks->Clear("C");
}

#ifdef GEANT
//_______________________________________________________________
void TpcLheHitsMaker::CheckTracks()
{

   using namespace std;

   Int_t       parent, tr_num, tr_pdg;
   MpdMCTrack *gtrack = 0;

   // AZ Int_t nMCtracks = fListMCtracks->GetEntries();
   Int_t nMCtracks = fListMCtracks->GetEntriesFast();

   Int_t selected_tracks = 0;

   // AZ for (Int_t ih = 0; ih < fLheHits->GetEntries(); ih++) {
   for (Int_t ih = 0; ih < fLheHits->GetEntriesFast(); ih++) {
      TpcLheHit *hit = (TpcLheHit *)fLheHits->UncheckedAt(ih);
      SetTrack(hit);
      //    hit->Print();
   }

   //#ifdef TRACK
   cout << " MC tracks in event: " << fListMCtracks->GetEntries();
   cout << " tracks in TPC " << fNTrack << endl;
   //#endif

   for (Int_t itrack = 0; itrack < fNTrack; itrack++) {

      TpcLheTrack *track = (TpcLheTrack *)fGeantTracks->UncheckedAt(itrack);
      tr_num             = track->GetTrackNumber();

      if (tr_num >= 0 && tr_num < nMCtracks) {
         tr_pdg = 0;
         parent = 99999;

         gtrack = (MpdMCTrack *)fListMCtracks->At(tr_num);

         tr_pdg = gtrack->GetPdgCode();
         track->SetPid(tr_pdg);
         track->SetCharge(TMath::Sign(1, tr_pdg));
         parent = gtrack->GetMotherId();
         if (parent == -1) {
            track->ComesFromMainVertex(kTRUE);
            track->SetGood(kTRUE);
         }

         TVector3 vertex;
         gtrack->GetStartVertex(vertex);
         track->SetVertex(vertex.X(), vertex.Y(), vertex.Z());

         TVector3 moment;
         gtrack->GetMomentum(moment);
         track->SetPx(moment.X());
         track->SetPy(moment.Y());
         track->SetPz(moment.Z());

         // Int_t nhit = track->GetNumberOfHits();
         // TObjArray  *chits = track->GetRHits();
         //  print some info about track

         /*
         cout << "\n *** track *** " << track->GetTrackNumber() <<
      " pdg " << track->GetPid() <<
      " nhits " <<nhit << " p_t " << track->GetPt() << endl;

         track->Print();
         TpcLheHit *hit = (TpcLheHit *)chits->First();
         hit->Print();
         */

         if (fCuts->IsGoodGeantTrack(track)) {
            selected_tracks++;
         }

         if (parent == 99999) cout << "  track parent unknown";
      }
   }

#if 0
  for (Int_t ih = 0; ih < fLheHits->GetEntries(); ih++) {
    TpcLheHit *hit = (TpcLheHit *)fLheHits->UncheckedAt(ih);
    Int_t mTrackNumber = hit->GetTrackID();
    for (Int_t itrack = 0; itrack < fNTrack; itrack++) {
      TpcLheTrack *track = (TpcLheTrack*)fGeantTracks->UncheckedAt(itrack);
      Int_t tr_num = track->GetTrackNumber();
      if (mTrackNumber ==  tr_num && !track->IsGood()) {
	hit->SetUsage(kTRUE);
      }
    }
  }
#endif

   cout << "Total number of tracks in TPC: " << setw(5) << fGeantTracks->GetEntries() << endl;
   cout << "           Good tracks in TPC: " << setw(5) << selected_tracks << endl;
}

//________________________________________________________________
void TpcLheHitsMaker::SetTrack(TpcLheHit *hit)
{

   Int_t mTrackNumber = hit->GetTrackID();
   Int_t mTN;

#ifdef HITS
   cout << " fNHit = " << hit->GetHitNumber();
   cout << " fNTrack = " << fNTrack;
   cout << " mTrackNumber = " << mTrackNumber;
   cout << endl;
#endif

   TpcLheTrack *track = NULL;

   Bool_t newtrack = kTRUE;

   if (fNTrack > 0) {
      for (Int_t i = 0; i < fGeantTracks->GetEntriesFast(); i++) {
         TpcLheTrack *currtrk = ((TpcLheTrack *)fGeantTracks->At(i));
         mTN                  = currtrk->GetTrackNumber();
         //      PR( currtrk->GetTrackNumber());
         if (mTN == mTrackNumber) {
            track    = currtrk;
            newtrack = kFALSE;
            break;
         }
      }
   }
   if (newtrack) track = AddTrack(mTrackNumber);

   track->AddHit(hit);
}

//________________________________________________________________
TpcLheTrack *TpcLheHitsMaker::AddTrack(Int_t mTrackNumber)
{
   // Add a new track to the list of tracks for this event.

   TClonesArray &tracks = *fGeantTracks;
   TpcLheTrack  *track  = new (tracks[fNTrack++]) TpcLheTrack(mTrackNumber);
   return track;
}

//_________________________________________________________________
void TpcLheHitsMaker::PrintTracks(Int_t ntr)
{

   //  cout << " " << fNTrack << " tracks  \n\n ";
   //  cout << "Station X Y Z Px Py Pz \n\n";

   Int_t ptracks;

   if (ntr != 0)
      ptracks = ntr;
   else
      ptracks = fGeantTracks->GetEntries();

   for (Int_t itrack = 0; itrack < ptracks; itrack++) {

      TpcLheTrack *track = (TpcLheTrack *)fGeantTracks->UncheckedAt(itrack);

      track->Print(); //
   }
}

#endif

ClassImp(TpcLheHitsMaker);
