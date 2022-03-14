// -------------------------------------------------------------------------
// -----                    MpdTpcDedxTask source file                 -----
// -----                 Created 9/02/10  by A. Zinchenko              -----
// -------------------------------------------------------------------------

/**  MpdTpcDedxTask.cxx
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** dE/dx determination in MPD TPC
 **/

#include "MpdTpcDedxTask.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdTpcHit.h"

#include "FairRootManager.h"
#include <TFile.h>
#include <TMath.h>

#include <iostream>
#include <set>

using std::cout;
using std::endl;
using std::set;

//__________________________________________________________________________
MpdTpcDedxTask::MpdTpcDedxTask(const char *name, Int_t iVerbose) : FairTask(name, iVerbose)
{
   fHistoDir = nullptr;
   fTracks = nullptr, fMCTracks = nullptr;
}

//__________________________________________________________________________
MpdTpcDedxTask::~MpdTpcDedxTask()
{
   // delete fKHits;
   // delete fTracks;
   // delete fTrackCand;
}

//__________________________________________________________________________
InitStatus MpdTpcDedxTask::Init()
{
   fTracks = 0x0;
   fTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("ItsTrack");
   if (fTracks == 0x0) fTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("TpcKalmanTrack");
   fHits0 = (TClonesArray *)FairRootManager::Instance()->GetObject("TpcHit");
   fHits  = (TClonesArray *)FairRootManager::Instance()->GetObject("TpcRecPoint");

   fMCTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("MCTrack");

   return kSUCCESS;
}

//__________________________________________________________________________
void MpdTpcDedxTask::Reset()
{
   ///
}

//__________________________________________________________________________
void MpdTpcDedxTask::SetParContainers() {}

//__________________________________________________________________________
void MpdTpcDedxTask::Finish()
{
   // Write();
}

//__________________________________________________________________________
void MpdTpcDedxTask::Exec(Option_t *option)
{

   const Double_t trunc        = 0.7; // truncation parameter
   static Int_t   eventCounter = 0;
   cout << " MpdTpcDedxTask event " << ++eventCounter << endl;

   Reset();

   Int_t nTracks = fTracks->GetEntriesFast();
   for (Int_t itr = 0; itr < nTracks; ++itr) {
      MpdTpcKalmanTrack *track = (MpdTpcKalmanTrack *)fTracks->UncheckedAt(itr);
      TClonesArray      *hits  = track->GetTrHits();
      Int_t              nHits = hits->GetEntriesFast(), nOK = 0;
      Double_t          *dedx = new Double_t[nHits];
      for (Int_t j = 0; j < nHits; ++j) {
         MpdKalmanHit *hit = (MpdKalmanHit *)hits->UncheckedAt(j);
         // cout << j << " " << hit->GetUniqueID() << " " << hit->GetDedx() << endl;
         if (hit->GetUniqueID() == 1) continue; // ITS hit
         if (hit->GetSignal() <= 0.) continue;  // strange
         // Apply flag cuts
         if (hit->GetFlag() & MpdKalmanHit::kSinglePix) continue;
         if (hit->GetFlag() & MpdKalmanHit::kSinglePad) continue;
         // if (hit->GetFlag() & MpdKalmanHit::kVirtual) continue;
         //
         Double_t sig = hit->GetSignal();
         if (!fHits) {
            // For hit producer
            MpdTpcHit *thit = (MpdTpcHit *)fHits0->UncheckedAt(hit->GetIndex());
            Double_t   x    = thit->GetStep();
            x               = TMath::Max(x, 0.4);
            x               = TMath::Min(x, 10.0);
            sig /= (7.57541 * TMath::Power(x - 0.292613, 0.0160641) - 6.64079);
         } else {
            // For clusters: correct for track angles
            while (1) {
               MpdKalmanHit *hit1 = nullptr;
               if (j == 0)
                  hit1 = (MpdKalmanHit *)hits->UncheckedAt(1);
               else
                  hit1 = (MpdKalmanHit *)hits->UncheckedAt(j - 1);
               Int_t isec  = hit->GetDetectorID() % 1000000;
               Int_t isec1 = hit1->GetDetectorID() % 1000000;
               if (isec != isec1 && TMath::Abs(isec - isec1) != 12) {
                  if (j == 0 || j == nHits - 1) break; // no correction
                  hit1  = (MpdKalmanHit *)hits->UncheckedAt(j + 1);
                  isec1 = hit1->GetDetectorID() % 1000000;
                  if (isec != isec1 && TMath::Abs(isec - isec1) != 12) break; // no correction
               }
               MpdTpcHit *thit  = (MpdTpcHit *)fHits->UncheckedAt(hit->GetIndex());
               MpdTpcHit *thit1 = (MpdTpcHit *)fHits->UncheckedAt(hit1->GetIndex());
               // Double_t padh = thit->GetStep();
               Double_t dx = thit->GetLocalX() - thit1->GetLocalX();
               Double_t dz = thit->GetLocalZ() - thit1->GetLocalZ();
               Double_t dy = thit->GetLocalY() - thit1->GetLocalY();
               // Double_t tancros = dx / padh;
               // Double_t tandip = dz / padh;
               Double_t tancros = TMath::Abs(dx / dy);
               Double_t tandip  = TMath::Abs(dz / dy);
               sig /= TMath::Sqrt(1 + tancros * tancros + tandip * tandip);
               // sig /= TMath::Power (1 + tancros*tancros + tandip*tandip, 0.2);
               // sig /= TMath::Sqrt (1 + tancros * tancros);
               break;
            }
         }
         if (fHits && sig > -800)
            dedx[nOK++] = sig; // threshold
         else if (!fHits)
            dedx[nOK++] = sig;
         hit->SetSignal(sig);
      }
      if (nOK == 0) continue;
      Int_t *indx = new Int_t[nOK];
      TMath::Sort(nOK, dedx, indx, kFALSE);
      Double_t sum    = 0.;
      Int_t    nTrunc = TMath::Nint(nOK * trunc);
      if (nTrunc == 0) nTrunc = 1;
      for (Int_t j = 0; j < nTrunc; ++j) sum += dedx[indx[j]];
      track->SetPartID(sum / nTrunc);
      delete[] dedx;
      delete[] indx;

      // Apply signal correction for the drift length
      if (fHits)
         // for (Int_t iter = 0; iter < 3; ++iter) DriftCorrection(track);
         for (Int_t iter = 0; iter < 1; ++iter) DriftCorrection(track);
   }
}

//__________________________________________________________________________

void MpdTpcDedxTask::DriftCorrection(MpdTpcKalmanTrack *track)
{
   // Signal correction for the drift length

   Double_t dedx = TMath::Log10(track->GetDedx(6.036e-3)); // in kev
   // cout << dedx << " " << track->GetPartID() << endl;
   dedx           = TMath::Min(dedx, 2.5);
   dedx           = TMath::Max(dedx, 1.0);
   Double_t dedx2 = dedx * dedx;
   Double_t coef  = -3.51299e-05 + 8.99592e-05 * dedx - 6.55681e-05 * dedx2 + 1.25936e-05 * dedx2 * dedx;
   // coef = TMath::Min (coef, 0.0);

   TClonesArray  *hits  = track->GetTrHits();
   Int_t          nhits = hits->GetEntriesFast();
   const Double_t trunc = 0.7; // truncation parameter
   set<Double_t>  signals;

   for (Int_t ih = 0; ih < nhits; ++ih) {
      MpdKalmanHit *hit = (MpdKalmanHit *)hits->UncheckedAt(ih);
      if (hit->GetSignal() < -800) continue;                    // threshold
      Double_t sig = TMath::Log10(hit->GetSignal() * 6.036e-3); // in keV
      // cout << ih << " " << " " << hit->GetMeas(1) << " " << sig << " ";
      // sig /= (1 + p01 * TMath::Abs(hit->GetMeas(1)) + p02 * hit->GetMeas(1) * hit->GetMeas(1));
      sig /= (1 + coef * hit->GetMeas(1) * hit->GetMeas(1));
      // cout << sig << endl;
      signals.insert(sig);
   }
   Int_t nTrunc = TMath::Nint(signals.size() * trunc);
   if (nTrunc == 0) nTrunc = 1;
   Int_t    nok = 0;
   Double_t sum = 0;

   for (set<Double_t>::iterator sit = signals.begin(); sit != signals.end(); ++sit) {
      sum += TMath::Power(10, *sit) / 6.036e-3;
      ++nok;
      if (nok == nTrunc) break;
   }
   track->SetPartID(sum / nTrunc);
   // cout << track->GetPartID() << endl;
}

//__________________________________________________________________________

void MpdTpcDedxTask::Write()
{
   /// Write

   TFile histoFile("Vertex.root", "RECREATE");
   Writedir2current(fHistoDir);
   histoFile.Close();
}

//__________________________________________________________________________
void MpdTpcDedxTask::Writedir2current(TObject *obj)
{
   /// Write

   if (!obj->IsFolder())
      obj->Write();
   else {
      TDirectory *cur = gDirectory;
      TDirectory *sub = cur->mkdir(obj->GetName());
      sub->cd();
      TList *listSub = ((TDirectory *)obj)->GetList();
      TIter  it(listSub);
      while (TObject *obj1 = it()) Writedir2current(obj1);
      cur->cd();
   }
}
//__________________________________________________________________________

ClassImp(MpdTpcDedxTask);
