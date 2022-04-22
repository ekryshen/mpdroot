/// \class MpdKfV0Fitter
///
/// Kalman filter V0 fitter for the MPD detector
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdKfV0Fitter.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdHelix.h"
#include "MpdTpcKalmanTrack.h"

#include "FairField.h"
#include "MpdMCTrack.h"
#include "FairRunAna.h"
#include "FairTask.h"

#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>

#include <iostream>
using std::cout;
using std::endl;

MpdKfV0Fitter *MpdKfV0Fitter::fgV0 = 0x0;

//__________________________________________________________________________
MpdKfV0Fitter::MpdKfV0Fitter() : FairTask()
{
   /// Default constructor
}

//__________________________________________________________________________
MpdKfV0Fitter::MpdKfV0Fitter(const char *name, const char *title) : FairTask(name)
{
   /// Constructor
   fgV0 = this;
}

//__________________________________________________________________________
MpdKfV0Fitter *MpdKfV0Fitter::Instance()
{
   /// Get pointer to the V0 fitter singleton object
   if (!fgV0) {
      fgV0 = new MpdKfV0Fitter;
      fgV0->Init();
      // Automatic destruction is forced
      std::atexit(DestroyInstance);
   }
   return fgV0;
}

//__________________________________________________________________________
MpdKfV0Fitter *MpdKfV0Fitter::Instance(const char *name, const char *title)
{
   /// Get pointer to the V0 fitter singleton object
   if (!fgV0) {
      fgV0 = new MpdKfV0Fitter(name, title);
      fgV0->Init();
      // Automatic destruction is forced
      std::atexit(DestroyInstance);
   }
   return fgV0;
}

//__________________________________________________________________________
MpdKfV0Fitter::~MpdKfV0Fitter()
{
   /// Destructor
   // FairRootManager *manager = FairRootManager::Instance();
   // manager->Write();
}

//__________________________________________________________________________
InitStatus MpdKfV0Fitter::Init()
{

   cout << "InitStatus MpdKfV0Fitter::Init\n\n";
   Double_t   bZ       = 5.;
   FairField *magField = FairRunAna::Instance()->GetField();
   if (!magField || TMath::Abs(magField->GetBz(0, 0, 0)) < 0.01) {
      cout << " -I- MpdKfV0Fitter::Init - Using the default constant magnetic field Bz = 5 kG " << endl;
   } else
      bZ = TMath::Abs(magField->GetBz(0, 0, 0));
   fieldConst = 0.3 * 0.01 * bZ / 10;
   return kSUCCESS;
}

//__________________________________________________________________________
InitStatus MpdKfV0Fitter::ReInit()
{
   // fMagField = FairRunAna::Instance()->GetField(); // !!! interim solution !!!
   return kSUCCESS;
}

//__________________________________________________________________________
void MpdKfV0Fitter::Reset()
{
   ///
   // cout << " MpdKfV0Fitter::Reset  " << endl;
}

//__________________________________________________________________________
void MpdKfV0Fitter::Register()
{
   ///
   // FairRootManager::Instance()->
   // Register("TpcLheKalmanTrack","TpcLheKalmanTrack", fTracks, kFALSE);
}

//__________________________________________________________________________
void MpdKfV0Fitter::Finish()
{
   ///
}

//__________________________________________________________________________
void MpdKfV0Fitter::Exec(Option_t *option)
{
   ///
}

//__________________________________________________________________________
void MpdKfV0Fitter::EvalVertex(MpdKalmanTrack *tr[2])
{
   /// Evaluate V0 position as PCA of 2 helices

   // Create MpdHelix's
   MpdHelix *helix[2];
   TVector3  sum;
   for (Int_t i = 0; i < 2; ++i) {
      Double_t r   = tr[i]->GetPosNew();
      Double_t phi = tr[i]->GetParam(0) / r;
      Double_t x   = r * TMath::Cos(phi);
      Double_t y   = r * TMath::Sin(phi);

      Double_t dip = tr[i]->GetParam(3);
      Double_t cur = TMath::Abs(tr[i]->GetParam(4)) * fieldConst;
      TVector3 o(x, y, tr[i]->GetParam(1));

      Int_t h = (Int_t)TMath::Sign(1.1, tr[i]->GetParam(4));
      // helix[itr++] = MpdHelix(cur, dip, phi, o, h);
      helix[i] = new MpdHelix(cur, dip, tr[i]->GetParam(2) - TMath::PiOver2() * h, o, h);
      // cout << helix[i]->xcenter() << " " << helix[i]->ycenter() << endl;
      sum += o;
   }
   pair<Double_t, Double_t> paths = helix[0]->pathLengths(*helix[1]);

   // cout << " Intersection: " << paths.first << " " << paths.second << endl;
   TVector3 p1 = helix[0]->at(paths.first);
   TVector3 p2 = helix[1]->at(paths.second);
   p1 += p2;
   p1 *= 0.5;
   if (TMath::Abs(paths.first) > 50) p1.SetXYZ(0, 0, 0);
   p1.GetXYZ(fVtx);
   for (Int_t i = 0; i < 2; ++i) delete helix[i];
   sum *= 0.5;
   // sum.GetXYZ(fVtx); // half-sum
}

//__________________________________________________________________________
Double_t MpdKfV0Fitter::FindVertex(MpdKalmanTrack *tr1, MpdKalmanTrack *tr2, TVector3 &vtx)
{
   /// Kalman filter based secondary vertex finder (fitter)

   const Int_t    nPass   = 3;     // number of iterations
   const Double_t cutChi2 = 1000.; // chi2-cut

   MpdKalmanTrack *tracks[2] = {tr1, tr2};

   EvalVertex(tracks);

   Int_t        nPrim = 0;
   TMatrixD     c(3, 3), cov(3, 3), xk0(3, 1), xk(3, 1), ck0(5, 1);
   TMatrixD     a(5, 3), b(5, 3);
   MpdKalmanHit hit;
   Double_t     chi2 = 0; //, chi2o = 0;

   xk.SetMatrixArray(fVtx);
   c(0, 0) = c(1, 1) = c(2, 2) = 1;

   for (Int_t ipass = 0; ipass < nPass; ++ipass) {

      // chi2o = chi2;
      chi2 = 0.;
      if (ipass == 0) cov = TMatrixD(TMatrixD::kInverted, c);

      for (Int_t itr = 0; itr < 2; ++itr) {
         xk0 = xk; // xk0 stands for x(k-1)
         cov = TMatrixD(TMatrixD::kInverted, c);

         MpdKalmanTrack *track = tracks[itr];

         // Select primaries
         // MpdMCTrack *mcTr = (MpdMCTrack*) fMCTracks->UncheckedAt(track->GetTrackID());
         // if (mcTr->GetMotherId() >= 0) continue; // secondary
         // Double_t th = TMath::PiOver2() - track->GetParam(3);
         // if (TMath::Abs(TMath::Log(TMath::Tan(th/2))) > 1.) continue; // eta-cut
         // if (1./TMath::Abs(track->GetParam(4)) < 0.2) continue; // pt-cut
         ++nPrim;

         MpdKalmanTrack track1 = *track;
         track1.SetParamNew(*track1.GetParam());
         track1.SetPos(track1.GetPosNew());
         track1.ReSetWeight();
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
         TMatrixD tmp42(cov, TMatrixD::kMult, dx);
         TMatrixD chi22(dx, TMatrixD::kTransposeMult, tmp42);
         // chi22.Print();
         if (chi21(0, 0) < 0 || chi22(0, 0) < 0) {
            cout << " !!! Chi2 < 0 " << chi21(0, 0) << " " << chi22(0, 0) << " " << ipass << " " << itr << " "
                 << track1.GetParamNew(4) << endl;
            // exit(0);
         }
         chi21 += chi22;
         chi2 += chi21(0, 0);
         // if (chi2 > cutChi2) {
         if (chi2 > cutChi2 || chi21(0, 0) < 0 || chi22(0, 0) < 0) {
            for (Int_t i = 0; i < 3; ++i) vtx[i] = fVtx[i];
            for (Int_t i = 0; i < 2; ++i) {
               // MpdTpcKalmanTrack trTmp = *tracks[i];
               MpdKalmanTrack trTmp = *tracks[i];
               trTmp.SetPos(trTmp.GetPosNew());
               trTmp.SetParamNew(*trTmp.GetParam());
               if (trTmp.GetNode() != "") trTmp.SetNode(""); // 25-jul-2012
               MpdKalmanFilter::Instance()->FindPca(&trTmp, fVtx);
               tracks[i]->SetParam(*trTmp.GetParamNew());
               tracks[i]->SetPosNew(trTmp.GetPosNew());
            }
            // AZ-Debug
            // cout << " !!! Too high chi2: " << ipass << " " << itr << " " << chi2 << " " << chi22(0,0) << " " <<
            // chi21(0,0) << " " << vtx[0] << " " << vtx[1] << " " << vtx[2] << " " << fVtx[0] << " " << fVtx[1] << " "
            // << fVtx[2] << " " << tracks[0]->GetNode() << " " << tracks[1]->GetNode() << " ids: " <<
            // tracks[0]->GetTrackID() << " " << tracks[1]->GetTrackID() << " " << xk0(0,0) << " " << xk0(1,0) << " " <<
            // xk0(2,0) << " " << xk(0,0) << " " << xk(1,0) << " " << xk(2,0) << endl;
            return cutChi2;
         }
      } // for (Int_t itr = 0; itr < 2;
   }    // for (Int_t ipass = 0; ipass < nPass;

   for (Int_t i = 0; i < 3; ++i) vtx[i] = fVtx[i] = xk(i, 0);
   // Covariance matrix
   fCovar.ResizeTo(3, 3);
   fCovar = c;

   Smooth(tracks);
   return chi2;
}

//__________________________________________________________________________
void MpdKfV0Fitter::Smooth(MpdKalmanTrack *tr[2])
{
   /// Smooth tracks from primary vertex (update momentum and track length -
   /// covariance matrix is not updated !!!)

   MpdKalmanHit hit;
   TMatrixD     c(3, 3), xk(3, 1), ck0(5, 1);
   TMatrixD     a(5, 3), b(5, 3);
   Double_t     rad = 0.;
   for (Int_t i = 0; i < 3; ++i) {
      xk(i, 0) = fVtx[i];
      if (i < 2) rad += fVtx[i] * fVtx[i];
   }
   rad          = TMath::Sqrt(rad);
   Double_t phi = TMath::ATan2(fVtx[1], fVtx[0]);

   for (Int_t itr = 0; itr < 2; ++itr) {
      MpdKalmanTrack *track = tr[itr];

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
      TMatrixD par = *track->GetParam();
      for (Int_t i = 0; i < 3; ++i) par(i + 2, 0) = qk(i, 0);
      par(0, 0) = rad * phi;
      par(1, 0) = fVtx[2];
      track->SetParam(par);
      track->SetPosNew(rad);

      // Update track length
      Double_t dLeng = track1.GetLength(); // track length from DCA to last saved position
      track1.SetParam(par);
      track1.SetPos(rad);
      track1.SetLength(0.);
      if (track->GetNode() == "")
         MpdKalmanFilter::Instance()->PropagateParamR(&track1, &hit, kTRUE);
      else
         MpdKalmanFilter::Instance()->PropagateParamP(&track1, &hit, kTRUE, kTRUE);

      track->SetLength(track->GetLength() - dLeng + track1.GetLength());
   } // for (Int_t itr = 0; itr < 2;
}

//__________________________________________________________________________
void MpdKfV0Fitter::ComputeAandB(TMatrixD &xk0, const MpdKalmanTrack *track, const MpdKalmanTrack &trackM, TMatrixD &a,
                                 TMatrixD &b, TMatrixD &ck0)
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
      // if (j == 4) shift *= TMath::Sign(1.,track0.GetParamNew(j)); // 1/p
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
void MpdKfV0Fitter::Proxim(const MpdKalmanTrack &track0, MpdKalmanTrack &track)
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

ClassImp(MpdKfV0Fitter);
