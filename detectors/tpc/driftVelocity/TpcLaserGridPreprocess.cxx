#include "TpcLaserGridPreprocess.h"

#include "FairRootManager.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include <Math/Point3D.h>
#include <Math/Rotation3D.h>
#include <Math/RotationZ.h>
#include <Math/Translation3D.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TFile.h>
#include <TRandom2.h>
#include <TSystem.h>
#include <TStyle.h>

#include "TpcPoint.h"
#include "MpdMCTrack.h"
#include "TaskHelpers.h"

#include <iostream>
#include <vector>

TpcLaserGridPreprocess::TpcLaserGridPreprocess()
   : FairTask("TPC laser grid preprocessor"), clrSecondaries(kTRUE), makeQA(kFALSE)
{
   inputBranchName = "TpcPoint";
}

//---------------------------------------------------------------------------
TpcLaserGridPreprocess::~TpcLaserGridPreprocess() {}

//---------------------------------------------------------------------------
InitStatus TpcLaserGridPreprocess::Init()
{
   FairRootManager *ioman = FairRootManager::Instance();

   if (!ioman) {
      std::cout << "\n-E- [TpcLaserGridPreprocess::Init]: RootManager not instantiated!" << std::endl;
      return kFATAL;
   }
   mcPointArray = (TClonesArray *)ioman->GetObject(inputBranchName);
   mcTrackArray = (TClonesArray *)ioman->GetObject("MCTrack");

   std::cout << "-I- TpcLaserGridPreprocess: Initialization successful." << std::endl;
   return kSUCCESS;
}

//---------------------------------------------------------------------------

void TpcLaserGridPreprocess::Exec(Option_t *opt)
{
   std::cout << "TpcLaserGridPreprocess::Exec started" << std::endl;

   Int_t nPoints = mcPointArray->GetEntriesFast();
   Int_t nTracks = mcTrackArray->GetEntriesFast();
   if (nPoints < 2) {
      Warning("TpcLaserGridPreprocess::Exec", "Not enough Hits (<2)");
      return;
   }

   std::cout << "Number of MC points is " << nPoints << std::endl;
   Int_t rmPrimPtsCount = 0, rmSecPtsCount = 0;

   for (Int_t ip = 0; ip < nPoints; ++ip) {
      TpcPoint *point   = (TpcPoint *)mcPointArray->UncheckedAt(ip);
      Int_t     ptTrack = point->GetTrackID();

      TVector3 l0, mom;
      if (ptTrack < beamsCount) {
         MpdMCTrack *tr = (MpdMCTrack *)mcTrackArray->UncheckedAt(ptTrack);
         tr->GetStartVertex(l0);
         tr->GetMomentum(mom);

         TVector3 xy0(0., 0., l0.Z());
         TVector3 l1(l0.X() + mom.X() * centerBeamMaxLength, l0.Y() + mom.Y() * centerBeamMaxLength, l0.Z());

         if (distance_to_line(xy0, l0, l1) < tubeInR) {
            if (point->GetLength() > centerBeamMaxLength) {
               mcPointArray->RemoveAt(ip);
               rmPrimPtsCount++;
               continue;
            }
         }
      } else if (clrSecondaries) {
         MpdMCTrack *tr = (MpdMCTrack *)mcTrackArray->UncheckedAt(point->GetTrackID());
         if (tr->GetMotherId() != -1 && tr->GetMotherId() < beamsCount) {
            mcPointArray->RemoveAt(ip);
            rmSecPtsCount++;
            continue;
         }
      }

      if (mcPointsSmearing && ptTrack < beamsCount) {
         ROOT::Math::XYZVector z0(0., 0., 1.);
         Double_t              cosPhi = mom.X() / (sqrt(mom.X() * mom.X() + mom.Y() * mom.Y()));
         Double_t              sinPhi = sqrt(1 - cosPhi * cosPhi);
         if (mom.Y() < 0) sinPhi = -sinPhi;
         Double_t cosTheta = mom.Z() / (sqrt(mom.X() * mom.X() + mom.Y() * mom.Y() + mom.Z() * mom.Z()));
         Double_t sinTheta = sqrt(1 - cosTheta * cosTheta);

         TRandom2 rnd;
         rnd.SetSeed();
         Double_t rndval = rnd.Gaus();
         while (abs(rndval) > 3.0 /*sigma*/) {
            rndval = rnd.Gaus();
         }
         Double_t              r = rndval * beamRadius;
         Double_t              a = rnd.Uniform(TMath::TwoPi());
         ROOT::Math::XYZVector p0(0., r, 0.);
         ROOT::Math::RotationZ rot(a);
         p0 = rot * p0;
         ROOT::Math::Rotation3D rotglo(cosTheta * cosPhi, -sinPhi, sinTheta * cosPhi, sinPhi * cosTheta, cosPhi,
                                       sinTheta * sinPhi, -sinTheta, 0., cosTheta);
         p0 = rotglo * p0;
         ROOT::Math::XYZVector trglo(point->GetX(), point->GetY(), point->GetZ());
         p0 = trglo + p0;

         point->SetXYZ(p0.X(), p0.Y(), p0.Z());
      }
   }

   std::cout << "Number of removed MC points: " << rmPrimPtsCount << " primary and " << rmSecPtsCount << " secondary"
             << std::endl;
   std::cout << "TpcLaserGridPreprocess::Exec finished" << std::endl;
}

void TpcLaserGridPreprocess::Finish()
{
   // fMCPointArray->Compess();

   if (makeQA) {
      toDirectory("QA/TPC/Laser");

      TGraph2D *gr = new TGraph2D();
      gr->SetName("LaserPreprocessed");
      gr->SetMarkerSize(0.25);
      gr->GetXaxis()->SetTitle("Z");
      gr->GetYaxis()->SetTitle("X");
      gr->GetZaxis()->SetTitle("Y");
      int grPts = 0;

      TGraph *beam_layer = new TGraph();
      beam_layer->SetName("LaserPreprocessedBeamLayer");
      beam_layer->SetTitle(";X [cm];Y [cm]");
      beam_layer->SetMarkerStyle(kFullCircle);
      beam_layer->SetMarkerColor(kRed);
      beam_layer->SetMarkerSize(0.6);
      int layPts = 0;

      const Int_t           num_layers              = 8;
      const Double_t        beam_layer_z_pos_offset = 30.;
      std::vector<Double_t> z_layers(num_layers);
      Int_t                 l = 0;
      for (Int_t j = -num_layers / 2; j < num_layers / 2; ++j) {
         z_layers[l] = j < 0 ? j * beam_layer_z_pos_offset : (j + 1) * beam_layer_z_pos_offset;
         l++;
      }
      const Double_t beamsLayer = 2;
      const Double_t eps_z      = 1.; // beam layer thickness

      Int_t nPoints = mcPointArray->GetEntriesFast();
      for (Int_t ip = 0; ip < nPoints; ++ip) {
         TpcPoint *point = (TpcPoint *)mcPointArray->UncheckedAt(ip);
         if (point) {
            Double_t gX = point->GetX();
            Double_t gY = point->GetY();
            Double_t gZ = point->GetZ();

            gr->SetPoint(grPts++, gZ, gX, gY);

            if (TMath::Abs(z_layers[beamsLayer] - gZ) < eps_z) beam_layer->SetPoint(layPts++, gX, gY);
         }
      }

      gr->Write();
      beam_layer->Write();
      gFile->cd();
   }
}

Double_t TpcLaserGridPreprocess::distance_to_line(TVector3 pt, TVector3 l1, TVector3 l2) const
{
   TVector3 AB   = l2 - l1;
   TVector3 AC   = pt - l1;
   Double_t area = (AB.Cross(AC)).Mag();
   Double_t CD   = area / AB.Mag();
   return CD;
}

ClassImp(TpcLaserGridPreprocess)
