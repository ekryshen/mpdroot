#include "TpcLheTrackFinder.h"

#include "TpcLheCMTrack.h"
#include "TpcLheCMPoint.h"
#include "lhe.h"

#include "FairMCApplication.h"
#include "FairTask.h"
#include "FairRunAna.h"
#include "TpcGeoPar.h"
#include "FairGeoNode.h"
#include "FairGeoVector.h"
#include "FairGeoMedium.h"
#include "FairRootManager.h"

#include "TObjectTable.h"
#include "TClonesArray.h"

ClassImp(TpcLheTrackFinder);

//________________________________________________________________
TpcLheTrackFinder::TpcLheTrackFinder()
{
   //---

   fFoundTracks = new TClonesArray("TpcLheTrack");
   fCMHits      = new TClonesArray("TpcLheCMPoint");
   fBench       = new TBenchmark();
   fVertex      = NULL;
}

//_________________________________________________________________
TpcLheTrackFinder::TpcLheTrackFinder(const char *name, const char *title) : FairTask(name)
{
   //---

   fFoundTracks = new TClonesArray("TpcLheTrack");
   fCMHits      = new TClonesArray("TpcLheCMPoint");
   fBench       = new TBenchmark();
   fVertex      = NULL;
}

//_________________________________________________________________
TpcLheTrackFinder::~TpcLheTrackFinder()
{

   if (fCMTracks) {
      fCMTracks->Delete();
      delete fCMTracks;
   }
   if (fCMHits) {
      fCMHits->Delete();
      delete fCMHits;
   }

   if (fSegments) delete fSegments;

   FairRootManager *fManager = FairRootManager::Instance();
   fManager->Write();
}

//_________________________________________________________________
void TpcLheTrackFinder::Register()
{
   //---
   FairRootManager::Instance()->Register("TpcLheTrack", "Lhe", fFoundTracks, kTRUE);

   FairRootManager::Instance()->Register("TpcLheCMPoint", "Lhe", fCMHits, kTRUE);
}

//________________________________________________________________
InitStatus TpcLheTrackFinder::Init()
{
   // ---

   FairRootManager *fManager = FairRootManager::Instance();
   fTrackCuts                = TpcLheTrackCuts::Instance();

   if (fOption.Contains("geant")) {
      fGeantTracks = (TClonesArray *)fManager->GetObject("LheGeantTrack");
   } else {
      fLheHits = (TClonesArray *)fManager->GetObject("LheHit");
   }

   Register();

   // create TObjArrays
   fCMTracks = new TObjArray(64);

   fVertex = new TpcLhePoint(0.0, 0.0, 0.0);

   fSegment = new TpcLheSegments();
   fSegment->Init();

   fSegments = fSegment->GetSegments();

   return kSUCCESS;
}

//_________________________________________________________________
void TpcLheTrackFinder::Exec(Option_t *option)
{

   if (fOption.Contains("geant")) {
      cout << "  Copy Geant tracks  \n" << endl;
      CopyClones(fGeantTracks, fFoundTracks);
      return;
   }

   if (fBench) {
      fBench->Start("finder");
   }

   Reset();

   Int_t n_hits = fLheHits->GetEntriesFast(); // number of hits

   TpcLheHit *ghit = NULL;

   Int_t good_hits = 0;

   for (Int_t ih = 0; ih < n_hits; ih++) {
      ghit                  = (TpcLheHit *)fLheHits->UncheckedAt(ih);
      TClonesArray  &cmhits = *fCMHits;
      TpcLheCMPoint *cmhit  = new (cmhits[good_hits++]) TpcLheCMPoint(ghit);
      cmhit->SetHitNumber(ghit->GetHitNumber());
      cmhit->Setup(fVertex);
      cmhit->SetUsage(kFALSE);
   }

   cout << " Working with " << good_hits << " hits\n";

   fSegment->FillSegments(fCMHits);

   DoTracking();

   if (fBench) {
      cout << endl;
      fBench->Show("finder");
   }
}

//________________________________________________________________
void TpcLheTrackFinder::AddTrackForFit(TpcLheCMTrack *track_in)
{
   //--- Add a new track to the list of tracks for this event.

   TClonesArray &tracks = *fFoundTracks;
   Int_t         size   = tracks.GetEntriesFast();

   TObjArray *rhits = (TObjArray *)track_in->GetCMHits();
   cout << "TpcLheTrackFinder::AddTrackForFit " << rhits << endl;
   if (rhits) {
      Int_t NHits = rhits->GetEntriesFast();

      TpcLheTrack *track = new (tracks[size]) TpcLheTrack(size);
      track->SetCharge(track_in->GetCharge());

      for (Int_t lh = 0; lh < NHits; lh++) {
         TpcLheHit *hit = dynamic_cast<TpcLheHit *>(rhits->At(lh));
         track->AddHit(hit);
      }
   }
}

//_________________________________________________________________
void TpcLheTrackFinder::Finish()
{

   cout << " found  " << fFoundTracks->GetEntries() << " tracks\n";

   FairRootManager *fManager = FairRootManager::Instance();
   fManager->Fill();
}

//________________________________________________________________
void TpcLheTrackFinder::Reset()
{
   //---

   if (fCMTracks->GetEntriesFast() != 0) fCMTracks->Clear("C");
   if (fCMHits->GetEntriesFast() != 0) fCMHits->Clear("C");
   if (fFoundTracks->GetEntriesFast() != 0) fFoundTracks->Clear("C");
}

//________________________________________________________________
void TpcLheTrackFinder::DoTracking()
{
   //--- Tracking

   for (Int_t im = HIGH; im <= HIGH; im++) {

      switch (im) {

      case HIGH: fTrackCuts->SetHighMomTrackCuts(); break;

      case LOW: fTrackCuts->SetLowMomTrackCuts(); break;

      default: fTrackCuts->SetHighMomTrackCuts();
      }

      LoopOverHits();

      char info[2] = {' ', '\0'};
      TrackingInfo(&info[0]);
   }
}

//________________________________________________________________
void TpcLheTrackFinder::CheckClones()
{
   //---

   TIter next(fCMTracks);
   next.Reset();
   TpcLheCMTrack *trk;
   Int_t          ntracks = fCMTracks->GetEntriesFast();

   next.Reset();
   while ((trk = (TpcLheCMTrack *)next())) {
      trk->SetGood(kTRUE);
   }

   for (Int_t it = 0; it < ntracks; it++) {
      trk = (TpcLheCMTrack *)fCMTracks->At(it);

      TpcLhePoint hel_ex = trk->GetCircle();

      Double_t xc_ex  = hel_ex.GetX();
      Double_t yc_ex  = hel_ex.GetY();
      Double_t rad_ex = hel_ex.GetZ();

      for (Int_t jt = it + 1; jt < ntracks; jt++) {
         TpcLheCMTrack *trk_in = (TpcLheCMTrack *)fCMTracks->At(jt);
         TpcLhePoint    hel_in = trk_in->GetCircle();

         Double_t xc_in  = hel_in.GetX();
         Double_t yc_in  = hel_in.GetY();
         Double_t rad_in = hel_in.GetZ();
         if (TMath::Abs(xc_ex - xc_in) < 1. && TMath::Abs(yc_ex - yc_in) < 1. && TMath::Abs(rad_ex - rad_in) < 1.)
            trk_in->SetGood(kFALSE);
         //      cout << "\n Clones\n"; trk->Print(); trk_in->Print();
      }
   }

   next.Reset();
   while ((trk = (TpcLheCMTrack *)next())) {
      if (trk->IsGood()) AddTrackForFit(trk);
   }
}

//________________________________________________________________
Bool_t TpcLheTrackFinder::TrackExtension(TpcLheCMTrack *track)
{

   TpcLheCMPoint *hit = (TpcLheCMPoint *)(track->GetRHits())->First();

   Bool_t back;

   if (hit->GetZ() < 0.) {
      back = kTRUE;
   } else {
      back = kFALSE;
   }

   GetClosestHit(track, 0, back);

   if (!fTrackCuts->IsGoodFoundTrack(track)) {
      // cleanup
      track->Clear();
      return kFALSE;
   } else {
      return kTRUE;
   }
}

//________________________________________________________________
void TpcLheTrackFinder::TrackingInfo(char *info)
{
   // Information about the tracking process.

   cout << " found ";

   cout.width(5);
   cout << fFoundTracks->GetEntriesFast() << " tracks ";
}

//________________________________________________________________
void TpcLheTrackFinder::LoopOverHits()
{
   //--- loops over all hits

   Int_t entries = fCMHits->GetEntriesFast();

   for (Int_t nhit = 0; nhit < entries; nhit++) {
      TpcLheCMPoint *hit = (TpcLheCMPoint *)fCMHits->At(nhit);
      hit->Print();
      // start hit was used before
      if (!hit->GetUsage()) {
         CreateTrack(hit);
      }
   }

   CheckClones();
}

//_______________________________________________________________________
void TpcLheTrackFinder::CreateTrack(TpcLheCMPoint *seed_hit)
{
   //---

   Double_t theta_cut  = 3. * .0002;
   Int_t    theta_low  = fSegment->GetThetaSegm(seed_hit->GetTheta() - theta_cut);
   Int_t    theta_high = fSegment->GetThetaSegm(seed_hit->GetTheta() + theta_cut);

   Double_t phi_cut  = 3. * .004;
   Int_t    phi_low  = fSegment->GetPhiSegm(seed_hit->GetPhi() - phi_cut);
   Int_t    phi_high = fSegment->GetPhiSegm(seed_hit->GetPhi() + phi_cut);

   for (Int_t theta_segm = theta_low; theta_segm <= theta_high; theta_segm++) {

      //     phi
      for (Int_t phi_segm = phi_low; phi_segm <= phi_high; phi_segm++) {

         // loop over entries in one sub segment
         TObjArray *cell = (TObjArray *)fSegments->At
                           //	(fSegment->GetSegm(theta_segm, phi_segm));
                           (fSegment->GetSegm(theta_segm, phi_segm));

         Int_t cell_entries = cell->GetEntriesFast();
         if (cell_entries) {

            for (Int_t hit_num = 0; hit_num < cell_entries; hit_num++) {
               TpcLheCMPoint *hit = (TpcLheCMPoint *)cell->At(hit_num);

               if (TMath::Abs(seed_hit->GetTheta() - hit->GetTheta()) > theta_cut) continue;
               if (TMath::Abs(seed_hit->GetPhi() - hit->GetPhi()) > phi_cut) continue;
               if (!hit->GetUsage() && hit->GetHitNumber() != seed_hit->GetHitNumber()) {

                  Int_t ntracks = fCMTracks->GetEntries();
                  fCMTracks->AddAtAndExpand(new TpcLheCMTrack(ntracks), ntracks);
                  TpcLheCMTrack *new_track = (TpcLheCMTrack *)fCMTracks->At(ntracks);
                  //	    fNTracks++;

                  new_track->AddPoint(seed_hit, kFALSE);
                  new_track->AddPoint(hit, kFALSE);

                  seed_hit->SetUsage(kTRUE);
                  hit->SetUsage(kTRUE);
                  TrackExtension(new_track);

               } //  if (!hit->GetUsage())
            }    //  over entries in segment
         }       //  over segments
      }          //  for (phi)
   }             //  for (theta)
}

//________________________________________________________________
void TpcLheTrackFinder::GetClosestHit(TpcLheCMTrack *cmtrack, Int_t cur_station, Bool_t back)
{
   //--- Returns the nearest hit in the next station

   TpcLheCMPoint *hit = NULL;

   TArrayI *philist   = new TArrayI();
   TArrayI *thetalist = new TArrayI();

   if (cmtrack->GetNumberOfHits() < 3) {
      GetPhiRange(philist, cmtrack);
      GetThetaRange(thetalist, cmtrack);
   }

   Int_t *jphi   = philist->GetArray();
   Int_t *jtheta = thetalist->GetArray();

   for (Int_t it = 0; it < thetalist->GetSize(); it++) {
      Int_t theta_segm = jtheta[it];

      for (Int_t ip = 0; ip < philist->GetSize(); ip++) {
         Int_t phi_segm = jphi[ip];

         // loop over entries in one sub segment
         TObjArray *cell = (TObjArray *)fSegments->At(fSegment->GetSegm(theta_segm, phi_segm));

         Int_t entries = cell->GetEntriesFast();
         if (entries) {
            for (Int_t sub_hit_num = 0; sub_hit_num < entries; sub_hit_num++) {
               hit = (TpcLheCMPoint *)cell->At(sub_hit_num);

               if (!hit->GetUsage()) { // hit was not used before
                  if (fTrackCuts->VerifyTrack(cmtrack, hit, back)) {
                     cmtrack->AddPoint(hit, kFALSE);
                     hit->SetUsage(kTRUE);
                     GetClosestHit(cmtrack, 0, back);
                  } else {
                  }
               } //  if (!hit->GetUsage())
            }    //  over entries
         }       //       //  over segments
         else {
         }
      } //  for (phi)
   }    //  for (theta)

   delete philist;
   delete thetalist;
}

//______________________________________________________________________
void TpcLheTrackFinder::GetThetaRange(TArrayI *list, TpcLheCMTrack *cmtrack)
{
   //---

   TObjArray     *trackpoint = cmtrack->GetCMHits();
   TpcLheCMPoint *hit        = (TpcLheCMPoint *)trackpoint->First();
   Double_t       thet       = hit->GetTheta();

   Int_t theta_high, theta_low;

   if (fSegment->GetThetaSegm(thet) < fSegment->GetNumThetaSegments() / 2) {
      theta_high = fSegment->GetThetaSegm(thet);
      theta_low  = 0;
   } else {
      theta_high = fSegment->GetNumThetaSegments() - 1;
      theta_low  = fSegment->GetThetaSegm(thet);
   }

   Int_t nTheta = theta_high - theta_low;

   list->Set(nTheta);

   for (Int_t s = 0; s < nTheta; s++) {
      list->AddAt(theta_low + s, s);
   }
}

//______________________________________________________________________
void TpcLheTrackFinder::GetPhiRange(TArrayI *cl1, TpcLheCMTrack *cmtrack)
{
   //---

   TObjArray     *trackpoint = cmtrack->GetCMHits();
   TpcLheCMPoint *hit        = (TpcLheCMPoint *)trackpoint->First();
   TpcLheCMPoint *hit1       = (TpcLheCMPoint *)trackpoint->Last();

   Double_t phi = hit->GetPhi();

   Double_t phi_c = fTrackCuts->GetPhiPrediction(cmtrack);

   Double_t phi_m = phi_c - TMath::Pi() / 2.;
   Double_t phi_p = phi_c + TMath::Pi() / 2.;

   Double_t OutRadTpc = .42;                         // m
   Double_t Rad       = cmtrack->GetRadius() / 100.; // m
   Double_t MagField  = 2;                           //    Tesla
   Double_t pt        = .2998 * Rad * MagField;      // PR(pt);

   TpcLhePoint vert  = cmtrack->GetVertex();
   Double_t    xtr_c = vert.GetX() / 100.; // m
   Double_t    ytr_c = vert.GetY() / 100.; // m

   Double_t xcross_1 = 0.;
   Double_t ycross_1 = 0.;
   Double_t xcross_2 = 0.;
   Double_t ycross_2 = 0.;
   Int_t crossing = fTrackCuts->Circle_Circle_Intersection(0., 0., OutRadTpc, xtr_c, ytr_c, Rad, &xcross_1, &ycross_1,
                                                           &xcross_2, &ycross_2);

   Int_t range, phi_low, phi_high;

   Double_t phi_cr;

   // distance from 1st point to line (vertex - 2nd point)
   Double_t dist = hit->GetX() * hit1->GetY() - hit->GetY() * hit1->GetX();

   if (dist > 0.) {
      // negative track -- counterclockwize
      cmtrack->SetCharge(-1);

      Double_t phi_r = phi_c - TMath::Pi() / 2.;
      if (phi_r < 0.) phi_r = TMath::TwoPi() + phi_r;
      phi_low = fSegment->GetPhiSegm(phi_r);

      if (crossing) {
         TVector2 v(xcross_2, ycross_2);
         phi_cr   = (v.Phi() >= 0.) ? phi_cr = v.Phi() : phi_cr = v.Phi() + 2. * TMath::Pi();
         phi_high = fSegment->GetPhiSegm(phi_cr);
      } else {
         phi_high = fSegment->GetPhiSegm(phi_c);
      }

      range = phi_high - phi_low + 1;

      if (phi_high < phi_low) {
         range = fSegment->GetNumPhiSegments() - phi_low + phi_high + 1;
      }
      cl1->Set(range);

      for (Int_t s = 0; s < range; s++) {
         //      cout << "range " << range << " cont " << phi_low + s << endl;
         Int_t segm = phi_low + s;
         if (segm >= fSegment->GetNumPhiSegments()) segm = segm - fSegment->GetNumPhiSegments();
         cl1->AddAt(segm, s);
      }

   } else {

      cmtrack->SetCharge(1);
      TVector2 v(xcross_1, ycross_1);
      phi_cr = (v.Phi() >= 0.) ? phi_cr = v.Phi() : phi_cr = v.Phi() + 2. * TMath::Pi();

      phi_low  = fSegment->GetPhiSegm(phi_cr);
      phi_high = fSegment->GetPhiSegm(phi_c + TMath::Pi() / 2.);

      if (crossing) {
         phi_low = fSegment->GetPhiSegm(phi_cr);
      } else {
         phi_low = fSegment->GetPhiSegm(phi_c);
      }

      range = phi_high - phi_low + 1;

      if (phi_low > phi_high) {
         range = fSegment->GetNumPhiSegments() - phi_low + phi_high + 1;
      }

      cl1->Set(range);

      for (Int_t s = 0; s < range; s++) {
         Int_t segm = phi_high - s;
         if (segm < 0) segm = fSegment->GetNumPhiSegments() + segm;
         cl1->AddAt(segm, s);
      }
   }
}

//______________________________________________________________________
void TpcLheTrackFinder::CopyClones(TClonesArray *cl1, TClonesArray *cl2)
{
   //---

   Int_t nEntries = cl1->GetEntriesFast();

   TClonesArray &clref = *cl2;
   TpcLheTrack  *old   = NULL;

   for (Int_t i = 0; i < nEntries; i++) {
      old = (TpcLheTrack *)cl1->At(i);
      new (clref[i]) TpcLheTrack(*old);
   }

   cout << " - " << cl2->GetEntriesFast() << " Geant tracks are copied.\n";

} //&:)=

/* ============================================================================
===============================================================================
===============================================================================
===============================================================================
===============================================================================
*/
