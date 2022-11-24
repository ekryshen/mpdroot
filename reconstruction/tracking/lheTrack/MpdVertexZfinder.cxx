// -------------------------------------------------------------------------
// -----                   MpdVertexZfinder source file                -----
// -----                 Created 19/08/09  by A. Zinchenko             -----
// -------------------------------------------------------------------------

/**  MpdVertexZfinder.h
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Primary vertex Z-position evaluator in MPD using hits from TPC
 **/

#include "MpdVertexZfinder.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
//#include "TpcPadPlane.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <TPad.h>
#include <TVector3.h>

#include <iostream>
#include <set>

using std::cout;
using std::endl;
using std::multiset;

// const Double_t MpdTrackFinderIts::fgkChi2Cut = 20; //20; //100;

//__________________________________________________________________________
MpdVertexZfinder::MpdVertexZfinder(BaseTpcSectorGeo &fSecGeo, const char *name, Int_t iVerbose)
   : FairTask(name, iVerbose)
{
   secGeo = dynamic_cast<TpcSectorGeoAZ *>(&fSecGeo);
   if (!secGeo) Fatal("MpdVertexZfinder::MpdVertexZfinder", " !!! Wrong geometry type !!! ");

   fKHits = nullptr, fhLays = nullptr;
   fhZ = nullptr, fUnc = nullptr, fhPhLay = nullptr;
}

//__________________________________________________________________________
MpdVertexZfinder::~MpdVertexZfinder()
{
   delete fhZ;
   delete fUnc;
   delete fhPhLay;
}

//__________________________________________________________________________
InitStatus MpdVertexZfinder::Init()
{
   // AZ fhZ = new TH1F("hZfinder","Zv",500,-100.2,99.8);
   fhZ     = new TH1F("hZfinder", "Zv", 300, -150., 150.);
   fhZDip  = new TH2F("hZDip", "Z-Dip", 300, -150., 150., 20, -TMath::Pi() / 3, TMath::Pi() / 3);
   fUnc    = new TF1("fUnc", "[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])+[3]", -10, 10);
   fhPhLay = new TH2F("hPhLay", "Phi-Layer", 40, -TMath::Pi(), TMath::Pi(), 60, 0, 60);

   // return ReInit();
   ReInit();

   return kSUCCESS;
}

//__________________________________________________________________________
InitStatus MpdVertexZfinder::ReInit()
{
   return kSUCCESS;
}

//__________________________________________________________________________
void MpdVertexZfinder::Reset()
{
   ///
   cout << " MpdVertexZfinder::Reset  " << endl;

   fhZ->Reset();
   fhZ->GetXaxis()->SetRange(10, 0); // reset axis range
   fhZDip->Reset();
   fhPhLay->Reset();
}

//__________________________________________________________________________
void MpdVertexZfinder::SetParContainers()
{
   // Get run and runtime database
   /*
   FairRunAna* run = FairRunAna::Instance();
   if ( ! run ) Fatal("SetParContainers", "No analysis run");

   FairRuntimeDb* db = run->GetRuntimeDb();
   if ( ! db ) Fatal("SetParContainers", "No runtime database");

   // Get STS geometry parameter container
   db->getContainer("MpdStsGeoPar");
   */
}

//__________________________________________________________________________
void MpdVertexZfinder::Finish()
{
   // Write();
}

//__________________________________________________________________________
void MpdVertexZfinder::Exec(Option_t *option)
{

   static int eventCounter = 0;
   cout << " - - - - \n Vertex Z-finder event " << ++eventCounter << endl;

   Reset();
}

//__________________________________________________________________________
void MpdVertexZfinder::SetHits(const TClonesArray *hits, const TH1F *hLays)
{
   /// Set hits container and histogram

   if (fKHits) return;
   fKHits = hits;
   fhLays = hLays;
}

//__________________________________________________________________________
Double_t MpdVertexZfinder::FindZ(const Int_t *layPointers, Int_t &flag)
{
   /// Evaluate vertex Z

   // Int_t layMax0 = ((MpdKalmanHit*)fKHits->First())->GetLayer();
   Int_t modular = 0;
   if (((MpdKalmanHit *)fKHits->First())->GetType() == MpdKalmanHit::kFixedP) modular = 1;
   // const TpcPadPlane *padPlane = TpcPadPlane::Instance();
   // MpdTpcSectorGeo* secGeo = MpdTpcSectorGeo::Instance();

   // Estimate Z-position of vertex
   // Loop over layers
   // Int_t layBeg = layMax0, layEnd = layBeg - 2, iDir = -1, dLays = 6; //10;
   Int_t iDir = 1, dLays = 5, isec = 0, isec1 = 0;
   // AZ-140921 Int_t layBeg = ((MpdKalmanHit*)fKHits->Last())->GetLayer(), layEnd = layBeg + 2;
   Int_t layBeg = ((MpdKalmanHit *)fKHits->Last())->GetLayer() + 6, layEnd = layBeg + 2;
   // if (fhLays->GetBinContent(layBeg+1,0) < 100) layEnd = layBeg + 6; //AZ-140921
   layEnd = ((MpdKalmanHit *)fKHits->First())->GetLayer() - dLays * iDir; // AZ-160921
   if (layBeg >= layEnd) {                                                // AZ-251021
      layBeg = ((MpdKalmanHit *)fKHits->Last())->GetLayer();
      layEnd = ((MpdKalmanHit *)fKHits->First())->GetLayer();
   }

   // Exclude some hit phase space regions (originating from loopers)
   // to regularize the task for low-multiplicity events
   Int_t nKhits = fKHits->GetEntriesFast();

   for (Int_t i = 0; i < nKhits; ++i) {
      MpdKalmanHit *hit = (MpdKalmanHit *)fKHits->UncheckedAt(i);
      fhPhLay->Fill(hit->GetPhi(), hit->GetLayer());
   }

   multiset<Int_t> contSet;

   for (Int_t i = 1; i <= fhPhLay->GetNbinsX(); ++i) {
      for (Int_t j = 1; j <= fhPhLay->GetNbinsY(); ++j) {
         Double_t cont = fhPhLay->GetBinContent(i, j);
         if (cont < 0.1) continue;
         contSet.insert(cont + 0.1);
      }
   }
   Int_t nok = 0, cutoff = 0, ntrunc = contSet.size() * 0.30; // 30% truncation
   for (multiset<Int_t>::iterator sit = contSet.begin(); sit != contSet.end(); ++sit) {
      ++nok;
      if (nok >= ntrunc) {
         cutoff = *sit;
         break;
      }
   }
   cutoff = cutoff + 3 * TMath::Sqrt(cutoff);
   cout << fhPhLay->GetMaximum() << " " << cutoff << endl;

   // Int_t layMax = layMax0, dLays = 6; //10;
   TVector3 pos1, pos2, posLoc;
   Double_t rad1, phi1, rad2, phi2;

   for (Int_t lay = layBeg; lay != layEnd; lay += iDir) {
      // for (Int_t lay = layMax; lay > layMax-2; --lay) {
      Int_t nHits1 = (Int_t)fhLays->GetBinContent(lay + 1, 0);
      Int_t nHits2 = (Int_t)fhLays->GetBinContent(lay + dLays * iDir + 1, 0);
      // cout << "Hits: " << nHits1 << " " << nHits2 << endl;
      //  Loop over hits in first layer
      for (Int_t ihit1 = 0; ihit1 < nHits1; ++ihit1) {
         MpdKalmanHit *hit1 = (MpdKalmanHit *)fKHits->UncheckedAt(layPointers[lay] + ihit1);
         if (hit1->GetFlag() & MpdKalmanHit::kUsed) continue;
         if (fhPhLay->GetBinContent(fhPhLay->FindBin(hit1->GetPhi(), hit1->GetLayer())) > cutoff) continue;

         if (modular) {
            posLoc.SetXYZ(hit1->GetMeas(0), hit1->GetDist(), hit1->GetMeas(1));
            isec = hit1->GetDetectorID() % 1000000;
            // pos1 = padPlane->fromSectorReferenceFrame(posLoc, isec);
            secGeo->Local2Global(isec, posLoc, pos1);
            rad1 = pos1.Pt();
            phi1 = pos1.Phi();
         } else {
            rad1 = hit1->GetPos();
            phi1 = hit1->GetMeas(0) / rad1;
         }

         // Loop over hits in second (closer to the beam) layer
         for (Int_t ihit2 = 0; ihit2 < nHits2; ++ihit2) {
            MpdKalmanHit *hit2 = (MpdKalmanHit *)fKHits->UncheckedAt(layPointers[lay + dLays * iDir] + ihit2);
            // if (h2->GetTrackID() != h1->GetTrackID()) continue; //!!! for test - exact ID matching

            if (hit2->GetFlag() & MpdKalmanHit::kUsed) continue;
            if (fhPhLay->GetBinContent(fhPhLay->FindBin(hit2->GetPhi(), hit2->GetLayer())) > cutoff) continue;

            if (modular) {
               posLoc.SetXYZ(hit2->GetMeas(0), hit2->GetDist(), hit2->GetMeas(1));
               isec1 = hit2->GetDetectorID() % 1000000;
               // if (TMath::Abs(isec-isec1) > 1 && TMath::Abs(isec-isec1) < padPlane->nSectors()-1) continue;
               // pos2 = padPlane->fromSectorReferenceFrame(posLoc, isec1);
               if (TMath::Abs(isec - isec1) > 1 && TMath::Abs(isec - isec1) < secGeo->NofSectors() - 1) continue;
               secGeo->Local2Global(isec1, posLoc, pos2);
               rad2 = pos2.Pt();
               // if ((rad2 - rad1) * iDir < 0.5) continue;
               phi2 = pos2.Phi();
            } else {
               rad2 = hit2->GetPos();
               // if (TMath::Abs(rad1 - rad2) < 0.5) continue;
               phi2 = hit2->GetMeas(0) / rad2;
            }
            Double_t dPhi = MpdKalmanFilter::Instance()->Proxim(phi1, phi2) - phi1;
            // AZ-150921 if (TMath::Abs(dPhi) > TMath::PiOver4()) continue;
            if (TMath::Abs(dPhi) > TMath::Pi() / 6) continue;
            // Double_t zvert = pos2.Z() - (pos1.Z()-pos2.Z()) / (rad1-rad2) * rad2;
            // AZ-150921 Double_t zvert = hit2->GetMeas(1) -
            //(hit1->GetMeas(1)-hit2->GetMeas(1)) / (rad1-rad2) * rad2;
            Double_t zvert = hit1->GetMeas(1) - (hit1->GetMeas(1) - hit2->GetMeas(1)) / (rad1 - rad2) * rad1;
            fhZ->Fill(zvert);
            fhZDip->Fill(zvert, TMath::ATan2(hit1->GetMeas(1) - hit2->GetMeas(1), TMath::Abs(rad1 - rad2)));
         }
      }
      if (fhZ->GetMaximum() > 50) break; // AZ-160921
   }                                     // for (Int_t lay = layBeg; lay != layEnd;

   // AZ-170921
   fhZDip->Divide(fhZDip);
   TH1D *projx = fhZDip->ProjectionX("px");
   projx->Scale(1 / projx->GetMaximum());
   fhZ->Multiply(projx);
   //

   Double_t z   = fhZ->GetXaxis()->GetBinCenter(fhZ->GetMaximumBin());
   Double_t max = fhZ->GetMaximum();
   fUnc->SetParameters(max, z, 5.0, 1.0);
   fUnc->SetParLimits(0, max * 0.2, max * 9.0);
   fUnc->SetParLimits(1, z - 5., z + 5.);
   fUnc->SetParLimits(2, 0.2, 10.);
   // fhZ->Fit(fUnc,"ww","",z-10.,z+10.);
   fhZ->Fit(fUnc, "wwNQ", "", z - 10., z + 10.);
   Double_t zFit  = fUnc->GetParameter(1);
   Double_t zVert = (TMath::Abs(z - zFit) > 2.0) ? z : zFit;
   // Set quality flag
   // AZ 14-09-21 if (max < 6 || (max < 10 && TMath::Abs(z-zFit) > 1)) { flag = -1; zVert = 0.0; }
   if (max < 6 || (max < 10 && TMath::Abs(z - zFit) > 1)) {
      flag  = -1;
      zVert = (z + zFit) / 2;
   } // AZ 14-09-21
   else
      flag = 0;
   cout << " MpdVertexZfinder::FindZ: " << z << " " << zFit << " " << zVert << " " << max << " " << flag << endl;
   /*
   fhZ->Draw("hist");
   //fhZDip->Draw("lego20");
   gPad->Update();
   Char_t go1[1];
   gets(go1);
   */
   return zVert;
}

ClassImp(MpdVertexZfinder);
