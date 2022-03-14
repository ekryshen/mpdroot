// -------------------------------------------------------------------------
// -----              MpdKfPrimaryVertexFinder source file             -----
// -----                 Created 4/07/09  by A. Zinchenko              -----
// -------------------------------------------------------------------------

/**  MpdKfPrimaryVertexFinder.h
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Primary vertex finder in MPD based on the Kalman filter
 **/

#include "MpdKfPrimaryVertexFinder.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanTrack.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdVertex.h"

#include "FairMCPoint.h"
#include "MpdMCTrack.h"
#include "FairRootManager.h"
#include <TAxis.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixFSym.h>
#include <TPad.h>
#include <TRandom.h>
#include <TFile.h>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

const Double_t MpdKfPrimaryVertexFinder::fgkChi2Cut = 3.5 * 3.5;
FILE          *lunVtx                               = 0x0; // fopen("vtx.dat","w");

//__________________________________________________________________________
MpdKfPrimaryVertexFinder::MpdKfPrimaryVertexFinder(const char *name, Int_t iVerbose)
   : FairTask(name, iVerbose), fConstrFlag(0), fSmoothSame(0), fVertTracks(nullptr)
{
   fHistoDir = nullptr;
   fHist[0] = fHist[1] = fHist[2] = nullptr;

   fCovar.ResizeTo(3, 3);
   for (Int_t i = 0; i < 3; ++i) fXYZ[i] = 0;
   fUnc = nullptr;
}

//__________________________________________________________________________
MpdKfPrimaryVertexFinder::~MpdKfPrimaryVertexFinder()
{
   // delete fKHits;
   // delete fTracks;
   // delete fTrackCand;
   for (Int_t i = 0; i < 3; ++i)
      if (fHist[i] != nullptr) fHist[i]->Delete();
   delete fUnc;
}

//__________________________________________________________________________
InitStatus MpdKfPrimaryVertexFinder::Init()
{

   fVertexCont = new TClonesArray("MpdVertex", 5);
   if (fConstrFlag) fVertTracks = new TClonesArray("MpdTpcKalmanTrack");

   fUnc = new TF1("fUnc", "[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])", -5, 5);

   fHist[0] = new TH1D("hXv", "Xv", 40, -2.05, 1.95);
   fHist[1] = new TH1D("hYv", "Yv", 40, -2.05, 1.95);
   // AZ fHist[2] = new TH1D("hZv","Zv",1000,-100.1,99.9);
   fHist[2] = new TH1D("hZv", "Zv", 1000, -150.0, 150.0);

   fTracks = 0x0;
   fTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("TpcTrackRefit");
   if (fTracks == 0x0) fTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("ItsTrack");
   if (fTracks == 0x0) fTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("ItsTrackKF");
   if (fTracks == 0x0) fTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("TpcKalmanTrack");

   fMCTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("MCTrack");

   FairRootManager::Instance()->Register("Vertex", "MPDVertex", fVertexCont, kTRUE);
   if (fConstrFlag) FairRootManager::Instance()->Register("PrimTracks", "ConstrTrs", fVertTracks, kTRUE);

   fNPass = 3;
   return kSUCCESS;
}

//__________________________________________________________________________
InitStatus MpdKfPrimaryVertexFinder::ReInit()
{
   return kSUCCESS;
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::Reset()
{
   ///
   cout << " MpdKfPrimaryVertexFinder::Reset  " << endl;

   if (fVertexCont) fVertexCont->Delete();
   if (fConstrFlag) fVertTracks->Delete();
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::SetParContainers() {}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::Finish()
{
   // Write();
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::Exec(Option_t *option)
{

   static int eventCounter = 0;
   cout << " VertexFinder event " << ++eventCounter << endl;

   Reset();

   EvalVertex();
   FindVertex();
   if (fConstrFlag) Smooth();
   Chi2Vertex();

   FillVertex(); // fill vertex info
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::EvalVertex()
{
   /// Primary vertex position evaluator

   Int_t    nTracks = fTracks->GetEntriesFast();
   Double_t dMax = 0., dMin = 0.09, xyzM[3], xyzF[3] = {0, 0, 0}, xyz[3] = {0};
   TAxis   *axis[3];
   for (Int_t i = 0; i < 3; ++i) {
      axis[i] = fHist[i]->GetXaxis();
      axis[i]->SetRange(10, 0); // reset axes ranges
      // AZ xyz[i] = fXYZ[i];
      fXYZ[i] = 0; // AZ-140921
   }
   fCovar = 0.;

   // Repeat until shift in transverse position becomes small
   for (Int_t iter = 0; iter < 5; ++iter) {

      for (Int_t i = 0; i < 3; ++i) fHist[i]->Reset();
      for (Int_t itr = 0; itr < nTracks; ++itr) {
         MpdKalmanTrack *track = (MpdKalmanTrack *)fTracks->UncheckedAt(itr);
         if (track->GetUniqueID() == 0 && track->GetNofHits() < 10) continue; // AZ-150921
         MpdKalmanTrack tr = *track;
         // if (iter > 0) {
         //  Find new PCA
         tr.SetParamNew(*tr.GetParam());
         tr.SetPos(tr.GetPosNew());
         // track1.ReSetWeight();
         MpdKalmanFilter::Instance()->FindPca(&tr, xyz);
         tr.SetParam(*tr.GetParamNew());
         //}
         Double_t r = tr.GetPosNew();
         // cout << " pos " << iter << " " << itr << " " << r << endl;
         // tr.GetParam()->Print();
         // Double_t phi = 0;
         Double_t phi = tr.GetParam(2);
         if (r > 1.e-7) phi = tr.GetParam(0) / r;
         fHist[0]->Fill(r * TMath::Cos(phi)); // x of PCA to (0,0,0)
         fHist[1]->Fill(r * TMath::Sin(phi));
         fHist[2]->Fill(tr.GetParam(1));
      }

      dMax = 0.;
      for (Int_t i = 0; i < 3; ++i) {
         // Take mean value
         xyzM[i] = fHist[i]->GetMean();
         // AZ fHist[i]->SetAxisRange(xyzM[i]-2.,xyzM[i]+2.);
         // AZ-170921 fHist[i]->SetAxisRange(xyzM[i]-5.,xyzM[i]+5.);
         fHist[i]->SetAxisRange(xyzM[i] - 10., xyzM[i] + 10.);
         xyzM[i] = fHist[i]->GetMean(); // restricted mean
         // AZ-150921 if (fHist[i]->GetMaximum() > 0) xyzM[i] = fHist[i]->GetBinCenter(fHist[i]->GetMaximumBin()); //
         // peak
         if (fHist[i]->GetMaximum() > 1) xyzM[i] = fHist[i]->GetBinCenter(fHist[i]->GetMaximumBin()); // peak
         // AZ Double_t rms = fHist[i]->Integral() > 1.1 ? TMath::Min (fHist[i]->GetRMS(), fHist[i]->GetMeanError() * 3)
         // : 5.0;
         Double_t rms = fHist[i]->Integral() > 1.1 ? TMath::Min(fHist[i]->GetRMS(), fHist[i]->GetMeanError() * 1) : 5.0;
         /*
         fUnc->SetParameters(fHist[i]->GetMaximum(),xyzM[i],rms/2);
         fUnc->SetParLimits(0,fHist[i]->GetMaximum()*0.8,fHist[i]->GetMaximum()*9);
         fUnc->SetParLimits(1,xyzM[i]-1.,xyzM[i]+1.);
         //fUnc->SetParLimits(2,0.0001,rms);
         fHist[i]->Fit(fUnc,"w","",xyzM[i]-2.,xyzM[i]+2.);
         xyzF[i] = fUnc->GetParameter(1);
         Double_t sigma = TMath::Abs (fUnc->GetParameter(2));
         */
         Double_t sigma = 999;
         if (sigma < rms) {
            xyz[i]       = xyzF[i];
            fCovar(i, i) = sigma * sigma;
         } else {
            xyz[i]       = xyzM[i];
            fCovar(i, i) = rms * rms;
         }
         cout << " Hist " << i << " " << fHist[i]->GetMaximum() << " " << xyzM[i] << " " << fXYZ[i] << " " << xyzF[i]
              << " " << rms << " " << sigma << endl;
         if (i == 2) continue;
         Double_t d = TMath::Abs(xyz[i] - fXYZ[i]);
         if (d > dMax) dMax = d;
      }
      for (Int_t i = 0; i < 3; ++i) fXYZ[i] = xyz[i];
      cout << iter << " " << dMax << endl;

      /*
      for (Int_t i = 0; i < 3; ++i) {
        fHist[i]->Draw();
        gPad->Update();
        Char_t go[1];
        gets(go);
      }
      */
      if (dMax < dMin) break;
   } // for (Int_t iter = 0; iter < 5;
   for (Int_t i = 0; i < 3; ++i) fCovar(i, i) = TMath::Max(fCovar(i, i), 0.02 * 0.02);

   // Use beam tolerances
   Double_t sigx = 9, sigy = 9, sigxy = 0.02; // sigma 200 um
   // AZ-150921 ??? while (sigx > 3*sigxy || sigy > 3*sigxy)
   while (TMath::Abs(sigx) > 3 || TMath::Abs(sigy) > 3) gRandom->Rannor(sigx, sigy);
   fXYZ[0]      = sigx * sigxy; // position
   fCovar(0, 0) = sigxy * sigxy;
   fXYZ[1]      = sigy * sigxy; // position
   fCovar(1, 1) = sigxy * sigxy;
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::FindVertex()
{
   /// Kalman primary vertex finder

   // AZ-040121 const Double_t chi2acc[3] = {35., 25., 15.};
   const Double_t chi2acc[3] = {10., 15., 20.};
   Int_t          nTracks = fTracks->GetEntriesFast(), nOK = 0, ok = 1;
   MpdVertex     *vtx = new ((*fVertexCont)[0]) MpdVertex();
   TMatrixD       c(3, 3), cov(3, 3), xk0(3, 1), xk(3, 1), ck0(5, 1);
   TMatrixD       a(5, 3), b(5, 3);
   MpdKalmanHit   hit;
   Double_t       chi2 = 0;

   // xk.Zero();
   // xk = 0.5;
   xk.SetMatrixArray(fXYZ);
   cout << " Tracks: " << nTracks << " " << xk(0, 0) << " " << xk(1, 0) << " " << xk(2, 0) << endl;
   // hit.SetR(35.); // 4 cm
   // xk(0,0) = 1;
   // cov(0,0) = cov(1,1) = 1;
   // cov(2,2) = 25 * 25;
   vector<Int_t> inds;

   for (Int_t ipass = 0; ipass < fNPass; ++ipass) {

      chi2             = 0.;
      Double_t cutChi2 = chi2acc[TMath::Min(ipass, 2)];
      // c.Zero();
      // c(0,0) = c(1,1) = c(2,2) = 100.;
      // c = fCovar;
      // cov = TMatrixD(TMatrixD::kInverted,c);
      if (ipass == 0) {
         c   = fCovar;
         cov = TMatrixD(TMatrixD::kInverted, c);
      }
      // AZ-040121 Int_t nPrim = 0, ibeg = 0, iend = nTracks, istep = 1;
      // if (ipass % 2 > 0) { ibeg = nTracks-1; iend = -1; istep = -1; }
      Int_t nPrim = 0, ibeg = nTracks - 1, iend = -1, istep = -1;
      if (ipass % 2 > 0) {
         ibeg  = 0;
         iend  = nTracks;
         istep = 1;
      }
      nOK = 0;

      for (Int_t itr = ibeg; itr != iend; itr += istep) {
         // for (Int_t itr = 0; itr < nTracks; ++itr) {
         // for (Int_t itr = nTracks-1; itr > -1; --itr) {
         if (ok) {
            xk0 = xk; // xk0 stands for x(k-1)
            cov = TMatrixD(TMatrixD::kInverted, c);
         }

         MpdKalmanTrack *track = (MpdKalmanTrack *)fTracks->UncheckedAt(itr);
         // if (track->GetNode() != "") continue; // exclude failed tracks
         //  Select primaries
         // MpdMCTrack *mcTr = (MpdMCTrack*) fMCTracks->UncheckedAt(track->GetTrackID());
         // if (mcTr->GetMotherId() >= 0) continue; // secondary
         // Double_t th = TMath::PiOver2() - track->GetParam(3);
         // if (TMath::Abs(TMath::Log(TMath::Tan(th/2))) > 1.) continue; // eta-cut
         // if (1./TMath::Abs(track->GetParam(4)) < 0.3) continue; // pt-cut
         if (track->GetUniqueID() == 0 && track->GetNofHits() < 10) continue; // AZ-150921: exclude short TPC tracks
         ++nPrim;

         MpdKalmanTrack track1 = *track;
         track1.SetParamNew(*track1.GetParam());
         track1.SetPos(track1.GetPosNew());
         track1.ReSetWeight();
         TMatrixD g = *track1.GetWeight(); // track weight matrix
         if (ipass == 0 && track->GetUniqueID())
            g *= 0.04; // AZ-050121 - downweight ITS track
         else if (ipass == 0)
            g *= 0.1; // AZ-140921 - downweight TPC track
         // track1.GetWeight()->Print();

         if (track->GetNode() == "") {
            hit.SetType(MpdKalmanHit::kFixedR);
            hit.SetPos(track->GetPos());
         } else {
            // 03-jan-2018 - skip track not pointing to the beam line
            ok = 0;
            continue;

            hit.SetType(MpdKalmanHit::kFixedP);
            TString detName = track->GetNode();
            if (track->GetUniqueID()) {
               // ITS
               detName = detName(16, detName.Length());
               detName += "#0";
            }
            MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();
            hit.SetDetectorID(geo->DetId(detName));
            // Find distance from the current track position to the last point (plane) -
            // to define direction (mainly for ITS)
            TVector3 pos   = geo->GlobalPos(&hit);
            TVector3 norm  = geo->Normal(&hit);
            Double_t v7[7] = {0.0};
            TString  node  = track1.GetNode();
            track1.SetNode("");
            MpdKalmanFilter::Instance()->SetGeantParamB(&track1, v7, 1);
            Double_t d = -(pos * norm); // Ax+By+Cz+D=0, A=nx, B=ny, C=nz
            TVector3 v3(v7[0], v7[1], v7[2]);
            d += v3 * norm;
            if (d < 0) track1.SetDirection(MpdKalmanTrack::kOutward);
         }

         MpdKalmanFilter::Instance()->PropagateToHit(&track1, &hit, kFALSE, kTRUE);
         // track1.GetParamNew()->Print();
         // cout << " Meas. radius " << track1.GetPosNew() << endl;

         // TMatrixD g = *track1.GetWeight(); // track weight matrix
         // track->GetCovariance()->Print();
         // track1.GetCovariance()->Print();
         // g.Print();
         // track1.Weight2Cov();
         // track1.GetCovariance()->Print();
         ComputeAandB(xk0, track, track1, a, b, ck0); // compute matrices of derivatives

         // W = (Bt*G*B)'
         TMatrixD tmp(g, TMatrixD::kMult, b);
         TMatrixD w(b, TMatrixD::kTransposeMult, tmp);
         w.Invert();

         // Gb = G - G*B*W*Bt*G
         TMatrixD tmp1(b, TMatrixD::kTransposeMult, g);
         TMatrixD tmp2(w, TMatrixD::kMult, tmp1);
         TMatrixD tmp3(b, TMatrixD::kMult, tmp2);
         TMatrixD tmp4(g, TMatrixD::kMult, tmp3);
         TMatrixD gb = g;
         gb -= tmp4;

         // Ck = ((Ck-1)' + At*Gb*A)'
         TMatrixD tmp5(gb, TMatrixD::kMult, a);
         c = TMatrixD(a, TMatrixD::kTransposeMult, tmp5);
         // cout << "gb, cov, c" << endl;
         // gb.Print();
         // cov.Print();
         // c.Print();
         c += cov;
         c.Invert();

         // xk = Ck*((Ck-1)'x(k-1)+At*Gb*(m-ck0))
         TMatrixD m = *track1.GetParamNew();
         m -= ck0; // m-ck0
         // cout << " m-ck0: " << endl;
         // ck0.Print();
         // m.Print();
         TMatrixD tmp11(gb, TMatrixD::kMult, m);
         TMatrixD tmp12(a, TMatrixD::kTransposeMult, tmp11);
         TMatrixD tmp13(cov, TMatrixD::kMult, xk0);
         tmp13 += tmp12;
         xk = TMatrixD(c, TMatrixD::kMult, tmp13);
         // cout << " Vertex: " << itr << " " << track->GetTrackID() << " " << xk(0,0) << " " << xk(1,0) << " " <<
         // xk(2,0) << endl;

         // qk = W*Bt*G*(m-ck0-A*xk)
         TMatrixD tmp21(a, TMatrixD::kMult, xk);
         tmp21 *= -1;
         tmp21 += m; // m-ck0-A*xk
         TMatrixD tmp22(g, TMatrixD::kMult, tmp21);
         TMatrixD tmp23(b, TMatrixD::kTransposeMult, tmp22);
         TMatrixD qk(w, TMatrixD::kMult, tmp23);

         // Residual m-ck0-A*xk-B*qk
         TMatrixD r = tmp21;
         TMatrixD tmp31(b, TMatrixD::kMult, qk);
         r -= tmp31;
         // cout << " Residual matrix: " << endl;
         // r.Print();
         // qk.Print();

         // Chi2 = rt*G*r + (xk-x(k-1))t*(Ck-1)'*(xk-x(k-1))
         TMatrixD tmp41(g, TMatrixD::kMult, r);
         TMatrixD chi21(r, TMatrixD::kTransposeMult, tmp41);
         // chi21.Print();
         TMatrixD dx = xk;
         dx -= xk0;
         // dx.Print();
         TMatrixD tmp42(cov, TMatrixD::kMult, dx);
         TMatrixD chi22(dx, TMatrixD::kTransposeMult, tmp42);
         // chi22.Print();
         if (chi21(0, 0) < 0 || chi22(0, 0) < 0) {
            cout << " FindVertex: chi2 < 0 " << chi21(0, 0) << " " << chi22(0, 0) << " " << ipass << " " << itr << endl;
            cout << track->GetNode() << " " << track->GetNodeNew() << " " << track->GetUniqueID() << " "
                 << track->Momentum() << " " << track->GetDirection() << endl;
            // exit(0);
            chi21(0, 0) = TMath::Max(chi21(0, 0), 0.);
            chi22(0, 0) = TMath::Max(chi22(0, 0), 0.);
         }
         chi21 += chi22;
         if (chi21(0, 0) > cutChi2) {
            ok = 0;
            continue;
         } // skip track
         chi2 += chi21(0, 0);
         ++nOK;
         ok = 1;
         if (ipass == fNPass - 1) inds.push_back(itr);
         // if (lunVtx) fprintf(lunVtx,"%10.3e \n",chi21(0,0));
         // chi21.Print();
         // c.Print();
         // TMatrixD www = c;
         // www.Invert();
         // www.Print();
      } // for (Int_t itr = 0; itr < nTracks;
      if (ok)
         cout << " Vertex position: " << ipass << " " << xk(0, 0) << " " << xk(1, 0) << " " << xk(2, 0) << " " << chi2
              << " " << nOK << " " << nPrim << endl;
      else
         cout << " Vertex position: " << ipass << " " << xk0(0, 0) << " " << xk0(1, 0) << " " << xk0(2, 0) << " "
              << chi2 << " " << nOK << " " << nPrim << endl;
      // fCovar = c;
      // xk.GetMatrix2Array(fXYZ);
   } // for (Int_t ipass = 0; ipass < fNPass;

   if (ok == 0) {
      // Take vertex after last accepted track
      c  = TMatrixD(TMatrixD::kInverted, cov);
      xk = xk0;
   }

   TMatrixFSym covar(3);
   for (Int_t i = 0; i < 3; ++i) {
      for (Int_t j = 0; j < 3; ++j) covar(i, j) = c(i, j);
   }
   vtx->SetVertex(xk(0, 0), xk(1, 0), xk(2, 0), chi2, 2 * nOK, nOK, covar);
   TArrayI *ind = vtx->GetIndices();
   ind->Set(nOK);
   for (Int_t i = 0; i < nOK; ++i) (*ind)[i] = inds[i];
   fCovar = c; // save for further use

   if (lunVtx) {
      for (Int_t i = 0; i < 3; ++i) {
         for (Int_t j = i; j < 3; ++j) fprintf(lunVtx, "%10.3e ", covar(i, j));
      }
      fprintf(lunVtx, "\n");
   } //*/
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::Smooth()
{
   /// Smooth tracks from primary vertex (update momentum and track length -
   /// covariance matrix is not updated !!!)

   MpdVertex   *vtx   = (MpdVertex *)fVertexCont->UncheckedAt(0);
   TArrayI     *ind   = vtx->GetIndices();
   Int_t        nPrim = ind->GetSize();
   MpdKalmanHit hit;
   TMatrixD     c(3, 3), xk(3, 1), ck0(5, 1);
   TMatrixD     a(5, 3), b(5, 3);
   TVector3     vert;
   vtx->Position(vert);
   xk(0, 0)     = vert.X();
   xk(1, 0)     = vert.Y();
   xk(2, 0)     = vert.Z();
   Double_t rad = vert.Pt();

   for (Int_t itr = 0; itr < nPrim; ++itr) {
      MpdKalmanTrack *track  = (MpdKalmanTrack *)fTracks->UncheckedAt((*ind)[itr]);
      MpdKalmanTrack *trVert = nullptr;
      if (fConstrFlag) trVert = new ((*fVertTracks)[itr]) MpdTpcKalmanTrack(*(MpdTpcKalmanTrack *)track);

      MpdKalmanTrack track1 = *track;
      track1.SetParamNew(*track1.GetParam());
      track1.SetPos(track1.GetPosNew());
      track1.ReSetWeight();
      track1.SetLength(0.);
      TMatrixD g = *track1.GetWeight(); // track weight matrix
      // track1.GetWeight()->Print();

      if (track->GetNode() == "") {
         hit.SetType(MpdKalmanHit::kFixedR);
         hit.SetPos(track->GetPos());
      } else {
         hit.SetType(MpdKalmanHit::kFixedP);
         TString detName = track->GetNode();
         if (track->GetUniqueID()) {
            // ITS
            detName = detName(16, detName.Length());
            detName += "#0";
         }
         MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();
         hit.SetDetectorID(geo->DetId(detName));
         // Find distance from the current track position to the last point (plane) -
         // to define direction (mainly for ITS)
         TVector3 pos   = geo->GlobalPos(&hit);
         TVector3 norm  = geo->Normal(&hit);
         Double_t v7[7] = {0.0};
         track1.SetNode("");
         MpdKalmanFilter::Instance()->SetGeantParamB(&track1, v7, 1);
         Double_t d = -(pos * norm); // Ax+By+Cz+D=0, A=nx, B=ny, C=nz
         TVector3 v3(v7[0], v7[1], v7[2]);
         d += v3 * norm;
         if (d < 0) track1.SetDirection(MpdKalmanTrack::kOutward);
      }

      MpdKalmanFilter::Instance()->PropagateToHit(&track1, &hit, kTRUE, kTRUE);
      // track1.GetParamNew()->Print();
      // cout << " Meas. radius " << track1.GetPosNew() << endl;

      ComputeAandB(xk, track, track1, a, b, ck0); // compute matrices of derivatives

      // W = (Bt*G*B)'
      TMatrixD tmp(g, TMatrixD::kMult, b);
      TMatrixD w(b, TMatrixD::kTransposeMult, tmp);
      w.Invert();

      TMatrixD m = *track1.GetParamNew();
      m -= ck0; // m-ck0

      // qk = W*Bt*G*(m-ck0-A*xk)
      TMatrixD tmp21(a, TMatrixD::kMult, xk);
      tmp21 *= -1;
      tmp21 += m; // m-ck0-A*xk
      TMatrixD tmp22(g, TMatrixD::kMult, tmp21);
      TMatrixD tmp23(b, TMatrixD::kTransposeMult, tmp22);
      TMatrixD qk(w, TMatrixD::kMult, tmp23);

      // Update momentum and last coordinates
      TMatrixD par = *track->GetParamNew();
      for (Int_t i = 0; i < 3; ++i) par(i + 2, 0) = qk(i, 0);
      par(0, 0) = rad * vert.Phi();
      par(1, 0) = vert.Z();
      if (trVert) {
         trVert->SetParam(par);
         trVert->SetPosNew(rad);
      } else if (fSmoothSame) {
         track->SetParam(par);
         track->SetPosNew(rad);
      }

      // Update track length
      Double_t dLeng = track1.GetLength(); // track length from DCA to last saved position
      track1.SetParam(par);
      track1.SetPos(rad);
      track1.SetLength(0.);
      if (track->GetNode() == "")
         MpdKalmanFilter::Instance()->PropagateParamR(&track1, &hit, kTRUE);
      else
         MpdKalmanFilter::Instance()->PropagateParamP(&track1, &hit, kTRUE, kTRUE);

      if (trVert)
         trVert->SetLength(track->GetLength() - dLeng + track1.GetLength());
      else if (fSmoothSame)
         track->SetLength(track->GetLength() - dLeng + track1.GetLength());
   } // for (Int_t itr = 0; itr < nPrim;
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::ComputeAandB(TMatrixD &xk0, const MpdKalmanTrack *track, const MpdKalmanTrack &trackM,
                                            TMatrixD &a, TMatrixD &b, TMatrixD &ck0)
{
   /// Compute matrices of derivatives w.r.t. vertex coordinates and track momentum

   Double_t vert0[3], *vert = xk0.GetMatrixArray(); // zero[3] = {0},
   for (Int_t i = 0; i < 3; ++i) vert0[i] = vert[i];

   MpdKalmanTrack trackk = *track;
   trackk.SetPos(trackk.GetPosNew());
   // trackk.GetParam()->Print();
   //  Propagate track to PCA w.r.t. point xk0
   MpdKalmanFilter::Instance()->FindPca(&trackk, vert0);
   // MpdKalmanFilter::Instance()->FindPca(&trackk,zero); // just for test
   // cout << trackk.GetPosNew() << endl;
   trackk.SetParam(*trackk.GetParamNew());
   // trackk.GetParam()->Print();

   // Put track at xk0
   Double_t r   = TMath::Sqrt(vert0[0] * vert0[0] + vert0[1] * vert0[1]);
   Double_t phi = trackk.GetParamNew(2); // track Phi
   if (r > 1.e-7) phi = TMath::ATan2(vert0[1], vert0[0]);
   trackk.SetPos(r);
   trackk.SetParam(0, r * phi);
   trackk.SetParam(1, vert0[2]);
   trackk.SetNode("");
   MpdKalmanTrack track0 = trackk;

   // Propagate track to chosen radius
   MpdKalmanHit hit;
   // hit.SetR(35.);
   // hit = *(MpdKalmanHitR*)track->GetTrHits()->Last();
   if (track->GetNode() == "") {
      hit.SetType(MpdKalmanHit::kFixedR);
      // hit.SetR(35.);
      // hit = *(MpdKalmanHitR*)track->GetTrHits()->Last();
      hit.SetPos(track->GetPos());
      MpdKalmanFilter::Instance()->PropagateParamR(&trackk, &hit, kFALSE);
      // trackk.GetParamNew()->Print();
      Proxim(trackM, trackk);
   } else {
      hit.SetType(MpdKalmanHit::kFixedP);
      TString detName = track->GetNode();
      if (track->GetUniqueID()) {
         // ITS
         detName = detName(16, detName.Length());
         detName += "#0";
      }
      MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();
      hit.SetDetectorID(geo->DetId(detName));
      // Find distance from the current track position to the last point (plane) -
      // to define direction (mainly for ITS)
      TVector3 pos   = geo->GlobalPos(&hit);
      TVector3 norm  = geo->Normal(&hit);
      Double_t v7[7] = {0.0};
      MpdKalmanFilter::Instance()->SetGeantParamB(&trackk, v7, 1);
      Double_t d = -(pos * norm); // Ax+By+Cz+D=0, A=nx, B=ny, C=nz
      TVector3 v3(v7[0], v7[1], v7[2]);
      d += v3 * norm;
      if (d < 0) trackk.SetDirection(MpdKalmanTrack::kOutward);
      MpdKalmanFilter::Instance()->PropagateParamP(&trackk, &hit, kFALSE, kTRUE);
      track0.SetDirection(trackk.GetDirection());
   }

   // Double_t shift = 0.01; // 100 um coordinate shift
   Double_t shift = 0.1; // 1 mm coordinate shift
   for (Int_t i = 0; i < 3; ++i) {
      MpdKalmanTrack track1 = track0;
      vert0[i] += shift;
      if (i > 0) vert0[i - 1] -= shift;
      r = TMath::Sqrt(vert0[0] * vert0[0] + vert0[1] * vert0[1]);
      if (r > 1.e-7)
         phi = TMath::ATan2(vert0[1], vert0[0]);
      else
         phi = track0.GetParamNew(2); // track Phi
      track1.SetPos(r);
      track1.SetParam(0, r * phi);
      track1.SetParam(1, vert0[2]);
      if (track->GetNode() == "") {
         MpdKalmanFilter::Instance()->PropagateParamR(&track1, &hit, kFALSE);
         Proxim(trackk, track1);
         // Proxim(track1,trackk);
      } else
         MpdKalmanFilter::Instance()->PropagateParamP(&track1, &hit, kFALSE, kTRUE);
      // Derivatives
      for (Int_t j = 0; j < 5; ++j) {
         a(j, i) = (track1.GetParamNew(j) - trackk.GetParamNew(j)) / shift;
      }
   }

   for (Int_t i = 0; i < 3; ++i) {
      MpdKalmanTrack track1 = track0;
      Int_t          j      = i + 2;
      shift                 = (*track->GetCovariance())(j, j);
      shift                 = TMath::Sqrt(shift);
      if (j == 4) shift *= TMath::Sign(1., -track0.GetParamNew(j)); // 1/p
      track1.SetParam(j, track0.GetParamNew(j) + shift);
      // if (j == 2 && track1.GetParamNew(j)*TMath::Sign(1.,track1.GetParamNew(j)) > TMath::Pi())
      // track1.SetParam(j,track0.GetParamNew(j)-shift);
      if (track->GetNode() == "") {
         MpdKalmanFilter::Instance()->PropagateParamR(&track1, &hit, kFALSE);
         Proxim(trackk, track1);
         // Proxim(track1,trackk);
      } else
         MpdKalmanFilter::Instance()->PropagateParamP(&track1, &hit, kFALSE, kTRUE);
      // Derivatives
      for (Int_t k = 0; k < 5; ++k) {
         b(k, i) = (track1.GetParamNew(k) - trackk.GetParamNew(k)) / shift;
      }
   }

   TMatrixD qk0(3, 1);
   for (Int_t i = 0; i < 3; ++i) qk0(i, 0) = track0.GetParamNew(i + 2);
   // qk0.Print();
   ck0 = *trackk.GetParamNew();
   ck0 -= TMatrixD(a, TMatrixD::kMult, xk0);
   ck0 -= TMatrixD(b, TMatrixD::kMult, qk0);
   // cout << " Derivatives: " << endl;
   // a.Print();
   // b.Print();
   // ck0.Print();
   // TMatrixD(b,TMatrixD::kMult,qk0).Print();
   // TMatrixD(a,TMatrixD::kMult,xk0).Print();
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::Proxim(const MpdKalmanTrack &track0, MpdKalmanTrack &track)
{
   /// Adjust track parameters

   if (track0.GetType() != MpdKalmanTrack::kBarrel) {
      cout << " !!! Implemented only for kBarrel tracks !!!" << endl;
      exit(0);
   }

   // Double_t tmp = track.GetParamNew(0);
   Double_t phi0 = track0.GetParamNew(0) / track0.GetPosNew();
   Double_t phi  = track.GetParamNew(0) / track.GetPosNew();
   phi           = MpdKalmanFilter::Instance()->Proxim(phi0, phi);
   TMatrixD *par = track.GetParamNew();
   (*par)(0, 0)  = phi * track.GetPosNew();
   phi0          = track0.GetParamNew(2);
   phi           = track.GetParamNew(2);
   phi           = MpdKalmanFilter::Instance()->Proxim(phi0, phi);
   (*par)(2, 0)  = phi;
   track.SetParamNew(*par);
   // cout << " Proxim: " << track0.GetParamNew(0) << " " << track.GetParamNew(0) << " " << tmp << endl;
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::Write()
{
   /// Write

   TFile histoFile("Vertex.root", "RECREATE");
   Writedir2current(fHistoDir);
   histoFile.Close();
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::Writedir2current(TObject *obj)
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
void MpdKfPrimaryVertexFinder::FillVertex()
{
   /// Fill vertex info
}

//__________________________________________________________________________
void MpdKfPrimaryVertexFinder::Chi2Vertex()
{
   /// Compute Chi2-distance of tracks from the primary vertex

   Int_t        nTracks = fTracks->GetEntriesFast(), nPrim = 0;
   MpdVertex   *vtx = (MpdVertex *)fVertexCont->UncheckedAt(0);
   TMatrixD     c(3, 3), cov(3, 3), xk0(3, 1), xk(3, 1), ck0(5, 1);
   TMatrixD     a(5, 3), b(5, 3);
   MpdKalmanHit hit;
   hit.SetType(MpdKalmanHit::kFixedR);
   TVector3 vert;
   vtx->Position(vert);
   xk0(0, 0) = vert.X();
   xk0(1, 0) = vert.Y();
   xk0(2, 0) = vert.Z();
   // xk(0,0) = xk(1,0) = xk(2,0) = 0.;
   // cov = fCovar;
   TMatrixFSym covM(3);
   vtx->CovMatrix(covM);
   cov = covM;
   // AZ cov.Invert();
   Int_t iii = 0;
   MpdKalmanFilter::Instance()->MnvertLocal(cov.GetMatrixArray(), 3, 3, 3, iii);

   // vvv.Print();
   // fCovar.Print();

   for (Int_t itr = 0; itr < nTracks; ++itr) {
      MpdKalmanTrack *track = (MpdKalmanTrack *)fTracks->UncheckedAt(itr);
      if (track->GetNode() != "") {
         track->SetChi2Vertex(999.);
         continue; // exclude failed tracks
      }

      // Select primaries
      // MpdMCTrack *mcTr = (MpdMCTrack*) fMCTracks->UncheckedAt(track->GetTrackID());
      // if (mcTr->GetMotherId() >= 0) continue; // secondary
      // Double_t th = TMath::PiOver2() - track->GetParam(3);
      // if (TMath::Abs(TMath::Log(TMath::Tan(th/2))) > 1.) continue; // eta-cut
      // if (1./TMath::Abs(track->GetParam(4)) < 0.2) continue; // pt-cut
      ++nPrim;

      // hit = *(MpdKalmanHitR*)track->GetTrHits()->Last();
      // MpdKalmanFilter::Instance()->GetGeo()->SetGlobalPos(&hit,TVector3(track->GetPos(),0.,0.),kTRUE);
      hit.SetPos(track->GetPos());
      MpdKalmanTrack track1 = *track;
      track1.SetParamNew(*track1.GetParam());
      track1.SetPos(track1.GetPosNew());
      track1.ReSetWeight();
      TMatrixD g = *track1.GetWeight(); // track weight matrix
      MpdKalmanFilter::Instance()->PropagateToHit(&track1, &hit, kFALSE);

      ComputeAandB(xk0, track, track1, a, b, ck0); // compute matrices of derivatives

      // W = (Bt*G*B)'
      TMatrixD tmp(g, TMatrixD::kMult, b);
      TMatrixD w(b, TMatrixD::kTransposeMult, tmp);
      w.Invert();

      // Gb = G - G*B*W*Bt*G
      TMatrixD tmp1(b, TMatrixD::kTransposeMult, g);
      TMatrixD tmp2(w, TMatrixD::kMult, tmp1);
      TMatrixD tmp3(b, TMatrixD::kMult, tmp2);
      TMatrixD tmp4(g, TMatrixD::kMult, tmp3);
      TMatrixD gb = g;
      gb -= tmp4;

      // Ck = ((Ck-1)' + At*Gb*A)'
      TMatrixD tmp5(gb, TMatrixD::kMult, a);
      c = TMatrixD(a, TMatrixD::kTransposeMult, tmp5);
      c += cov;
      c.Invert();

      // xk = Ck*((Ck-1)'x(k-1)+At*Gb*(m-ck0))
      TMatrixD m = *track1.GetParamNew();
      m -= ck0; // m-ck0
      TMatrixD tmp11(gb, TMatrixD::kMult, m);
      TMatrixD tmp12(a, TMatrixD::kTransposeMult, tmp11);
      TMatrixD tmp13(cov, TMatrixD::kMult, xk0);
      tmp13 += tmp12;
      xk = TMatrixD(c, TMatrixD::kMult, tmp13);

      // qk = W*Bt*G*(m-ck0-A*xk)
      TMatrixD tmp21(a, TMatrixD::kMult, xk);
      tmp21 *= -1;
      tmp21 += m; // m-ck0-A*xk
      TMatrixD tmp22(g, TMatrixD::kMult, tmp21);
      TMatrixD tmp23(b, TMatrixD::kTransposeMult, tmp22);
      TMatrixD qk(w, TMatrixD::kMult, tmp23);

      // Residual m-ck0-A*xk-B*qk
      TMatrixD r = tmp21;
      TMatrixD tmp31(b, TMatrixD::kMult, qk);
      r -= tmp31;

      // Chi2 = rt*G*r + (xk-x(k-1))t*(Ck-1)'*(xk-x(k-1))
      TMatrixD tmp41(g, TMatrixD::kMult, r);
      TMatrixD chi21(r, TMatrixD::kTransposeMult, tmp41);
      TMatrixD dx = xk;
      dx -= xk0;
      TMatrixD tmp42(cov, TMatrixD::kMult, dx);
      TMatrixD chi22(dx, TMatrixD::kTransposeMult, tmp42);
      // cout << " Chi2Vertex: " << chi21(0,0) << " " << chi22(0,0) << endl;
      chi21 += chi22;
      // chi21.Print();

      track->SetChi2Vertex(chi21(0, 0));
   } // for (Int_t itr = 0; itr < nTracks;
}

ClassImp(MpdKfPrimaryVertexFinder);
