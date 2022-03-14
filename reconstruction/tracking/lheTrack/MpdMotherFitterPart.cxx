/// \class MpdMotherFitterPart
///
/// Kalman filter mother particle fitter for the MPD detector (using MpdParticle)
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdMotherFitterPart.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdHelix.h"
#include "MpdVertex.h"
#include "MpdMCTrack.h"

#include "FairField.h"
#include "FairRunAna.h"
#include "FairTask.h"

#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>

#include <iostream>
using std::cout;
using std::endl;

MpdMotherFitterPart *MpdMotherFitterPart::fgMF = 0x0;

//__________________________________________________________________________
MpdMotherFitterPart::MpdMotherFitterPart() : FairTask()
{
   /// Default constructor

   fCovar.ResizeTo(3, 3);
}

//__________________________________________________________________________
MpdMotherFitterPart::MpdMotherFitterPart(const char *name, const char *title) : FairTask(name)
{
   /// Constructor

   fgMF = this;
   fCovar.ResizeTo(3, 3);
}

//__________________________________________________________________________
MpdMotherFitterPart *MpdMotherFitterPart::Instance()
{
   /// Get pointer to the mother fitter singleton object
   if (!fgMF) {
      fgMF = new MpdMotherFitterPart;
      fgMF->Init();
      // automatic destruction is forced
      std::atexit(DestroyInstance);
   }
   return fgMF;
}

//__________________________________________________________________________
MpdMotherFitterPart *MpdMotherFitterPart::Instance(const char *name, const char *title)
{
   /// Get pointer to the mother fitter singleton object
   if (!fgMF) {
      fgMF = new MpdMotherFitterPart(name, title);
      fgMF->Init();
      // automatic destruction is forced
      std::atexit(DestroyInstance);
   }
   return fgMF;
}

//__________________________________________________________________________
MpdMotherFitterPart::~MpdMotherFitterPart()
{
   /// Destructor
   // FairRootManager *manager = FairRootManager::Instance();
   // manager->Write();
}

//__________________________________________________________________________
InitStatus MpdMotherFitterPart::Init()
{

   cout << "InitStatus MpdMotherFitterPart::Init\n\n";
   Double_t   bZ       = 5.;
   FairField *magField = FairRunAna::Instance()->GetField();
   if (!magField || TMath::Abs(magField->GetBz(0, 0, 0)) < 0.01) {
      cout << " -I- MpdMotherFitterPart::Init - Using the default constant magnetic field Bz = 5 kG " << endl;
   } else
      bZ = TMath::Abs(magField->GetBz(0, 0, 0));
   fieldConst = 0.3 * 0.01 * bZ / 10;
   return kSUCCESS;
}

//__________________________________________________________________________
InitStatus MpdMotherFitterPart::ReInit()
{
   // fMagField = FairRunAna::Instance()->GetField(); // !!! interim solution !!!
   return kSUCCESS;
}

//__________________________________________________________________________
void MpdMotherFitterPart::Reset()
{
   ///
   // cout << " MpdMotherFitterPart::Reset  " << endl;
}

//__________________________________________________________________________
void MpdMotherFitterPart::Register()
{
   ///
   // FairRootManager::Instance()->
   // Register("TpcLheKalmanTrack","TpcLheKalmanTrack", fTracks, kFALSE);
}

//__________________________________________________________________________
void MpdMotherFitterPart::Finish()
{
   ///
}

//__________________________________________________________________________
void MpdMotherFitterPart::Exec(Option_t *option)
{
   ///
}

//__________________________________________________________________________

Double_t MpdMotherFitterPart::BuildMother(MpdParticle *mother, vector<MpdParticle *> &vDaught)
{
   /// Build mother particle from daughters which were smoothed
   /// according to the decay vertex constraint (after FindVertex)

   TVector3 vtx;
   Double_t chi2 = MpdMotherFitterPart::Instance()->FindVertex(vDaught, vtx);
   if (chi2 < -1.0) return chi2; // failed to find decay vertex (too high chi2)
   Int_t nDaught = vDaught.size();
   // Check for sanity
   for (Int_t id = 0; id < nDaught; ++id) {
      Double_t theta = vDaught[id]->Theta();
      if (theta < 0 || theta > TMath::Pi()) return -chi2; // weird track
   }

   mother->SetChi2(chi2);

   TVector3 mom3;
   Int_t    charge = 0;
   Double_t energy = 0.0;

   for (Int_t i = 0; i < nDaught; ++i) {
      MpdParticle *part = vDaught[i];
      charge += part->GetCharge();
      mom3 += part->Momentum3();
      part->FillJ();
      Double_t ptot = part->Momentum3().Mag();
      energy += TMath::Sqrt(part->GetMass() * part->GetMass() + ptot * ptot);
      mother->AddDaughter(part->GetIndx());
   }
   mother->SetCharge(charge);
   mother->SetMass(TMath::Sqrt(energy * energy - mom3.X() * mom3.X() - mom3.Y() * mom3.Y() - mom3.Z() * mom3.Z()));
   mother->Setx(vDaught[0]->Getx());
   TMatrixD qm(3, 1);
   qm(0, 0) = mom3.Phi();
   qm(1, 0) = mom3.Theta();
   if (charge == 0)
      qm(2, 0) = mom3.Mag();
   else
      qm(2, 0) = -fieldConst / mom3.Pt() * TMath::Abs(charge);
   mother->Setq(qm);

   ParamsAtDca(mother); // compute params at DCA
   // return chi2; //

   mother->FillJinv(mom3);

   // Compute covariance matrix
   TMatrixD en(3, 3);
   TMatrixD ec(3, 3);
   for (Int_t i = 0; i < nDaught; ++i) {
      MpdParticle *part = vDaught[i];
      // E += E*Jt
      TMatrixD jt  = TMatrixD(TMatrixD::kTransposed, part->GetJ());
      TMatrixD tmp = TMatrixD(part->GetE(), TMatrixD::kMult, jt);
      if (part->GetCharge()) {
         // Charged track
         ec += tmp;
      } else {
         // Neutral
         en += tmp;
      }
   }

   TMatrixD etot = ec;
   etot += en;
   TMatrixD c    = GetCovariance();
   TMatrixD qtot = ComputeQmatr(vDaught);

   TMatrixD ck0(5, 1);
   // ComputeAandB(mother->Getx(), *mother, fA, fB, ck0);
   ComputeAandB(mother->Getx(), *mother, mother->GetA(), mother->GetB(), ck0);

   // Covar. matrix
   TMatrixD at(TMatrixD::kTransposed, mother->GetA());
   TMatrixD tmp11(c, TMatrixD::kMult, at);
   TMatrixD tmp12(mother->GetA(), TMatrixD::kMult, tmp11);

   TMatrixD bt(TMatrixD::kTransposed, mother->GetB());
   TMatrixD tmp21(mother->GetJinv(), TMatrixD::kTransposeMult, bt);
   TMatrixD tmp22(etot, TMatrixD::kMult, tmp21);
   TMatrixD tmp23(mother->GetA(), TMatrixD::kMult, tmp22);

   TMatrixD tmp31(etot, TMatrixD::kTransposeMult, at);
   TMatrixD tmp32(mother->GetJinv(), TMatrixD::kMult, tmp31);
   TMatrixD tmp33(mother->GetB(), TMatrixD::kMult, tmp32);

   TMatrixD tmp41(mother->GetJinv(), TMatrixD::kTransposeMult, bt);
   TMatrixD tmp42(qtot, TMatrixD::kMult, tmp41);
   TMatrixD tmp43(mother->GetJinv(), TMatrixD::kMult, tmp42);
   TMatrixD tmp44(mother->GetB(), TMatrixD::kMult, tmp43);

   TMatrixD gm(5, 5);
   gm = tmp12;
   gm += tmp23;
   gm += tmp33;
   gm += tmp44;
   // cout << " Mother covariance " << endl;
   // fG.Print();
   gm.Invert(); // mother weight
   mother->SetG(gm);

   return chi2;
}

//__________________________________________________________________________

Double_t MpdMotherFitterPart::BuildMother(MpdParticle *mother, vector<MpdKalmanTrack *> &vTracks,
                                          vector<MpdParticle *> &vDaught)
{
   /// Build mother particle from daughters which were smoothed
   /// according to the decay vertex constraint (after FindVertex).
   /// Daughters are built from tracks and parametrized at their
   /// intersection point.

   vDaught.clear();
   vector<MpdHelix> vhel;

   // Create 2 helices and cross them
   for (Int_t itr = 0; itr < 2; ++itr) {
      Double_t rad = vTracks[itr]->GetPosNew();
      Double_t phi = 0.0;
      if (rad > 1.e-6) phi = (*vTracks[itr]->GetParam())(0, 0) / rad;
      vhel.push_back(MpdHelix(vTracks[itr]->Momentum3(),
                              TVector3(rad * TMath::Cos(phi), rad * TMath::Sin(phi), (*vTracks[itr]->GetParam())(1, 0)),
                              vTracks[itr]->Charge()));
   }
   pair<Double_t, Double_t> paths = vhel[0].pathLengths(vhel[1]);
   TVector3                 cross;
   Double_t                 xyz[3] = {0};
   Int_t                    ntr    = vTracks.size();

   if (paths.first < 100.0 || paths.second < 100.0) {
      // MpdKalmanHit hit;
      // hit.SetType(MpdKalmanHit::kFixedR);
      cross = vhel[0].at(paths.first);
      cross += vhel[1].at(paths.second);
      cross *= 0.5;
      cross.GetXYZ(xyz);
      /*
      for (Int_t itr = 0; itr < ntr; ++itr) {
        //Double_t s = vhel[itr].pathLength(cross);
        //TVector3 pca = vhel[itr].at(s);
        //hit.SetPos(pca.Pt());
        //MpdKalmanFilter::Instance()->PropagateToHit(vTacks[itr],&hit,kFALSE);
      }
      */
   }

   // xyz[0] = xyz[1] = xyz[2] = 0;
   for (Int_t itr = 0; itr < ntr; ++itr) {
      MpdKalmanFilter::Instance()->FindPca(vTracks[itr], xyz);
      vTracks[itr]->SetParam(*vTracks[itr]->GetParamNew());
      vTracks[itr]->Weight2Cov(); // AZ
      vDaught.push_back(new MpdParticle(*vTracks[itr], -1, 0.1396, xyz));
   }
   Double_t chi2 = BuildMother(mother, vDaught);
   // Bring back to master frame in transverse plane
   for (Int_t j = 0; j < 2; ++j) mother->Getx()(j, 0) += xyz[j];
   return chi2;
}

//__________________________________________________________________________

void MpdMotherFitterPart::EvalVertex(vector<MpdParticle *> vDaught)
{
   /// Evaluate decay vertex position as PCA of helices

   // Create MpdHelix's
   Int_t      nDaught = vDaught.size(), i0 = 0;
   MpdHelix **helix = new MpdHelix *[nDaught];
   TVector3   sum, p1, p2, p0;
   for (Int_t i = 0; i < nDaught; ++i) {
      Double_t dip = TMath::PiOver2() - vDaught[i]->Theta();
      Double_t cur = TMath::Abs(vDaught[i]->GetMeas(4));
      if (vDaught[i]->GetCharge() == 0) cur = numeric_limits<double>::epsilon();
      Int_t h = (Int_t)TMath::Sign(1.1, vDaught[i]->GetMeas(4));
      // Double_t dca = vDaught[i]->Dca();
      Double_t phase = vDaught[i]->GetMeas(2) - TMath::PiOver2() * h;
      // Double_t x = dca * TMath::Sin(vDaught[i]->Phi());
      // Double_t y = -dca * TMath::Cos(vDaught[i]->Phi());
      Double_t x = vDaught[i]->GetXY(0);
      Double_t y = vDaught[i]->GetXY(1);
      TVector3 o(x, y, vDaught[i]->GetMeas(1));
      helix[i] = new MpdHelix(cur, dip, phase, o, h);
      // o.Print();
      //*
      if (vDaught[i]->GetCharge() == 0) {
         p0.SetXYZ(vDaught[i]->Getx()(0, 0), vDaught[i]->Getx()(1, 0), vDaught[i]->Getx()(2, 0));
         i0 = 1;
         // p0.SetXYZ(0,0,0);
         //*
         TVector3 dcaV(x, y, vDaught[i]->GetMeas(1));
         // p0.Print();
         // dcaV.Print();
         p0 += dcaV;
         p0 *= 0.5;
         //*/
      }
      //*/
   }

   Int_t    ncombs = 0, nD1 = nDaught - 1;
   Double_t pathMax = 0.0;
   for (Int_t i = 0; i < nD1; ++i) {
      for (Int_t j = i + 1; j < nDaught; ++j) {
         // if (TMath::Abs(helix[j]->curvature()) <= numeric_limits<double>::epsilon()) {
         // straight line
         // Double_t dx = helix[j]
         pair<Double_t, Double_t> paths = helix[i]->pathLengths(*helix[j]);
         p1                             = helix[i]->at(paths.first);
         p2                             = helix[j]->at(paths.second);
         // cout << " Intersection: " << helix[i]->period() << " " << paths.first << " " << helix[j]->period() << " " <<
         // paths.second << " " << (p1-p2).Mag() << endl;
         sum += (p1 + p2);
         ncombs += 2;
         pathMax = TMath::Max(pathMax, TMath::Abs(paths.first));
         pathMax = TMath::Max(pathMax, TMath::Abs(paths.second));
      }
   }

   if (pathMax > 80.)
      p1.SetXYZ(0, 0, 0);
   else
      p1 = sum * (1. / ncombs);
   if (i0)
      p0.GetXYZ(fVtx);
   else
      p1.GetXYZ(fVtx);
   for (Int_t i = 0; i < nDaught; ++i) delete helix[i];
   // sum *= 0.5;
   // sum.GetXYZ(fVtx); // half-sum
   delete[] helix;
}

//__________________________________________________________________________
Double_t MpdMotherFitterPart::FindVertex(vector<MpdParticle *> vDaught, TVector3 &vtx)
{
   /// Kalman filter based secondary vertex fitter

   const Int_t    nPass   = 3;     // number of iterations
   const Double_t cutChi2 = 1000.; // chi2-cut

   EvalVertex(vDaught);

   Int_t        nDaught = vDaught.size();
   TMatrixD     c(3, 3), cov(3, 3), xk0(3, 1), xk(3, 1), ck0(5, 1);
   TMatrixD     a(5, 3), b(5, 3);
   MpdKalmanHit hit;
   Double_t     chi2 = 0; //, chi2o = 0;

   xk.SetMatrixArray(fVtx);
   // xk += 1.0; // test
   c(0, 0) = c(1, 1) = c(2, 2) = 1;

   for (Int_t ipass = 0; ipass < nPass; ++ipass) {

      // chi2o = chi2;
      chi2 = 0.;
      // c.Zero(); // new 05.06.17
      c(0, 0) = c(1, 1) = c(2, 2) = 100; // new 05.06.17
      if (ipass == 0) cov = TMatrixD(TMatrixD::kInverted, c);

      Int_t ibeg = 0, iend = nDaught, istep = 1;
      if (ipass % 2 > 0) {
         ibeg  = nDaught - 1;
         iend  = -1;
         istep = -1;
      }

      for (Int_t itr = ibeg; itr != iend; itr += istep) {
         xk0 = xk; // xk0 stands for x(k-1)
         cov = TMatrixD(TMatrixD::kInverted, c);

         MpdParticle *part = vDaught[itr];

         ComputeAandB(xk0, *part, a, b, ck0); // compute matrices of derivatives
         TMatrixD g = part->GetG();

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
         TMatrixD m = part->GetMeas();
         // m.Print();
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
         // Double_t dxMax = dx.Max();
         // if (dxMax > 7*TMath::Sqrt(cov.Max())) return -cutChi2; //AZ-050121
         TMatrixD tmp42(cov, TMatrixD::kMult, dx);
         TMatrixD chi22(dx, TMatrixD::kTransposeMult, tmp42);
         // chi22.Print();
         if (chi21(0, 0) < 0 || chi22(0, 0) < 0) {
            cout << " !!! Chi2 < 0 " << chi21(0, 0) << " " << chi22(0, 0) << " " << ipass << " " << itr << " "
                 << part->GetMeas(4) << endl;
            // exit(0);
         }
         chi21 += chi22;
         chi2 += chi21(0, 0);
         // if (chi2 > cutChi2) {
         if (chi2 > cutChi2 || chi21(0, 0) < 0 || chi22(0, 0) < 0) {
            for (Int_t i = 0; i < 3; ++i) vtx[i] = fVtx[i];
            /*
            for (Int_t i = 0; i < nDaught; ++i) {
              MpdKalmanTrack trTmp = *vDaught[i];
              trTmp.SetPos(trTmp.GetPosNew());
              trTmp.SetParamNew(*trTmp.GetParam());
              if (trTmp.GetNode() != "") trTmp.SetNode(""); // 25-jul-2012
              MpdKalmanFilter::Instance()->FindPca(&trTmp,fVtx);
              vDaught[i]->SetParam(*trTmp.GetParamNew());
              vDaught[i]->SetPosNew(trTmp.GetPosNew());
            }
            */
            // AZ-Debug
            // cout << " !!! Too high chi2: " << ipass << " " << itr << " " << chi2 << " " << chi22(0,0) << " " <<
            // chi21(0,0) << " " << vtx[0] << " " << vtx[1] << " " << vtx[2] << " " << fVtx[0] << " " << fVtx[1] << " "
            // << fVtx[2] << " " << tracks[0]->GetNode() << " " << tracks[1]->GetNode() << " ids: " <<
            // tracks[0]->GetTrackID() << " " << tracks[1]->GetTrackID() << " " << xk0(0,0) << " " << xk0(1,0) << " " <<
            // xk0(2,0) << " " << xk(0,0) << " " << xk(1,0) << " " << xk(2,0) << endl;
            return -cutChi2;
         }
         // cout << ipass << " " << itr << " " << chi2 << endl;
      } // for (Int_t itr = 0; itr < nDaught;
   }    // for (Int_t ipass = 0; ipass < nPass;

   for (Int_t i = 0; i < 3; ++i) vtx[i] = fVtx[i] = xk(i, 0);
   fCovar = c;

   Smooth(vDaught);
   return chi2;
}

//__________________________________________________________________________
void MpdMotherFitterPart::Smooth(vector<MpdParticle *> vDaught)
{
   /// Smooth particle momentum and corresponding covariance matrices

   TMatrixD c(3, 3), xk(3, 1), ck0(5, 1);
   TMatrixD a(5, 3), b(5, 3);

   xk.SetMatrixArray(fVtx);
   Int_t nDaught = vDaught.size();

   for (Int_t itr = 0; itr < nDaught; ++itr) {

      MpdParticle *part = vDaught[itr];

      TMatrixD g = part->GetG();          // track weight matrix
      ComputeAandB(xk, *part, a, b, ck0); // compute matrices of derivatives

      // W = (Bt*G*B)'
      TMatrixD tmp(g, TMatrixD::kMult, b);
      TMatrixD w(b, TMatrixD::kTransposeMult, tmp);
      w.Invert();

      TMatrixD m = part->GetMeas();
      m -= ck0; // m-ck0

      // qk = W*Bt*G*(m-ck0-A*xk)
      TMatrixD tmp21(a, TMatrixD::kMult, xk);
      tmp21 *= -1;
      tmp21 += m; // m-ck0-A*xk
      TMatrixD tmp22(g, TMatrixD::kMult, tmp21);
      TMatrixD tmp23(b, TMatrixD::kTransposeMult, tmp22);
      TMatrixD qk(w, TMatrixD::kMult, tmp23);

      // Update momentum and last coordinates
      /*
      TMatrixD par = *track->GetParam();
      for (Int_t i = 0; i < 3; ++i) par(i+2,0) = qk(i,0);
      par(0,0) = rad * phi;
      par(1,0) = fVtx[2];
      track->SetParam(par);
      track->SetPosNew(rad);

      // Update track length
      Double_t dLeng = track1.GetLength(); // track length from DCA to last saved position
      track1.SetParam(par);
      track1.SetPos(rad);
      track1.SetLength(0.);
      if (track->GetNode() == "") MpdKalmanFilter::Instance()->PropagateParamR(&track1,&hit,kTRUE);
      else MpdKalmanFilter::Instance()->PropagateParamP(&track1,&hit,kTRUE,kTRUE);

      track->SetLength (track->GetLength() - dLeng + track1.GetLength());
      */

      // Save all matrices for further usage
      part->Setq(qk);
      part->Setx(xk);
      part->SetA(a);
      part->SetB(b);
      part->SetC(fCovar);
      part->SetW(w);

      // E = -C*At*G*B*W
      TMatrixD tmp31(b, TMatrixD::kMult, w);
      TMatrixD tmp32(g, TMatrixD::kMult, tmp31);
      TMatrixD tmp33(a, TMatrixD::kTransposeMult, tmp32);
      TMatrixD e(fCovar, TMatrixD::kMult, tmp33);

      // D = W + W*Bt*G*A*C*At*G*B*W = W + W*Bt*G*A*(-E)
      TMatrixD tmp41(a, TMatrixD::kMult, e);
      e *= -1.0;
      part->SetCovE(e);
      TMatrixD tmp42(g, TMatrixD::kMult, tmp41);
      TMatrixD tmp43(b, TMatrixD::kTransposeMult, tmp42);
      TMatrixD d(w, TMatrixD::kMult, tmp43);
      d += w;
      part->SetCovD(d);

   } // for (Int_t itr = 0; itr < nDaught;
}

//__________________________________________________________________________
void MpdMotherFitterPart::ComputeAandB(const TMatrixD &xk0, const MpdParticle &part, TMatrixD &a, TMatrixD &b,
                                       TMatrixD &ck0, Bool_t flag)
{
   /// Compute matrices of derivatives w.r.t. vertex coordinates and particle momentum

   if (part.GetCharge() == 0) {
      // Neutral particle
      ComputeAB0(xk0, part, a, b, ck0, flag);
      return;
   }

   // Put track at xk0 - find rotation angle from xk0 to DCA (i.e. gamma)
   // Center of original circle:
   Double_t mx = part.GetXY(0) - TMath::Sin(part.GetMeas(2)) / part.GetMeas(4);
   Double_t my = part.GetXY(1) + TMath::Cos(part.GetMeas(2)) / part.GetMeas(4);

   Double_t x0    = xk0(0, 0);
   Double_t y0    = xk0(1, 0);
   Double_t phDca = TMath::ATan2(part.GetXY(1) - my, part.GetXY(0) - mx);
   Double_t phXk  = TMath::ATan2(y0 - my, x0 - mx);
   Double_t gam   = phDca - phXk;
   gam            = TMath::ASin(TMath::Sin(gam)); // AZ-230821
   Double_t ph0   = part.GetMeas(2) - gam;        // phi at xk0

   // Center of shifted circle
   Double_t cosph = TMath::Cos(ph0);
   Double_t sinph = TMath::Sin(ph0);
   mx             = x0 - sinph / part.GetMeas(4);
   my             = y0 + cosph / part.GetMeas(4);
   // Recompute gamma angle
   Double_t kinv = 1. / part.GetMeas(4);
   // Double_t gam1 = TMath::ATan ((x0 * cosph + y0 * sinph) / (x0 * sinph - y0 * cosph - kinv));

   // Gamma derivatives
   Double_t mxmy2 = mx * mx + my * my;
   Double_t dgdx  = -my / mxmy2;
   Double_t dgdy  = mx / mxmy2;
   Double_t dgdk  = -(x0 * cosph + y0 * sinph) / mxmy2 * kinv * kinv;
   Double_t dgdph = y0 * dgdx - x0 * dgdy;

   // A = dHc / dx
   a               = 0.0;
   Double_t cosgam = TMath::Cos(gam);
   Double_t singam = TMath::Sin(gam);
   Double_t val    = kinv - x0 * sinph + y0 * cosph;
   Double_t costh  = TMath::Cos(part.GetMeas(3));
   Double_t sinth  = TMath::Sin(part.GetMeas(3));

   a(0, 0) = (cosgam * sinph - singam * dgdx * val) / cosgam / cosgam;  // d(dca) / dx0
   a(0, 1) = -(cosgam * cosph + singam * dgdy * val) / cosgam / cosgam; // d(dca) / dy0

   a(1, 0) = dgdx * costh / sinth * kinv; // dz / dx
   a(1, 1) = dgdy * costh / sinth * kinv; // dz / dy
   a(1, 2) = 1.0;                         // dz / dz

   a(2, 0) = dgdx; // d(phi) / dx
   a(2, 1) = dgdy; // d(phi) / dy

   // B = dHc / dq
   b       = 0.0;
   b(0, 0) = ((x0 * cosph + y0 * sinph) * cosgam - singam * dgdph * val) / cosgam / cosgam; // d(dca) / dphi
   b(0, 2) = -kinv * kinv + (kinv * kinv * cosgam - singam * dgdk * val) / cosgam / cosgam; // d(dca) / dk

   b(1, 0) = dgdph * costh / sinth * kinv;               // dz / dphi
   b(1, 1) = -gam / sinth / sinth * kinv;                // dz / dtheta
   b(1, 2) = kinv * cosph / sinph * (dgdk - kinv * gam); // dz / dk

   b(2, 0) = 1 + dgdph; // dphi / dphi
   b(2, 2) = dgdk;      // dphi / dk

   b(3, 1) = 1.0; // d(theta) / d(theta)

   b(4, 2) = 1.0; // dk / dk
   if (!flag) return;

   // ck0
   ck0(0, 0) = kinv - val / cosgam;
   ck0(1, 0) = xk0(2, 0) + gam * costh / sinth * kinv;
   // ck0(2,0) = ph0 + gam1;
   ck0(2, 0) = part.GetMeas(2);
   ck0(3, 0) = part.GetMeas(3);
   ck0(4, 0) = part.GetMeas(4);

   TMatrixD qk0(3, 1);
   qk0(0, 0) = ph0;
   qk0(1, 0) = part.GetMeas(3);
   qk0(2, 0) = part.GetMeas(4);
   ck0 -= TMatrixD(a, TMatrixD::kMult, xk0);
   ck0 -= TMatrixD(b, TMatrixD::kMult, qk0);
}

//__________________________________________________________________________
void MpdMotherFitterPart::ComputeAB0(const TMatrixD &xk0, const MpdParticle &part, TMatrixD &a, TMatrixD &b,
                                     TMatrixD &ck0, Bool_t flag)
{
   /// Compute matrices of derivatives w.r.t. vertex coordinates and particle momentum
   /// (for neutrals)

   // A = dHn / dx
   a              = 0.0;
   Double_t cosph = TMath::Cos(part.Phi());
   Double_t sinph = TMath::Sin(part.Phi());
   Double_t sinth = TMath::Sin(part.Theta());
   // Double_t tanth = TMath::Tan(part.Theta());
   Double_t costh = TMath::Cos(part.Theta());
   Double_t coTan = costh / sinth;

   a(0, 0) = sinph;  // d(dca) / dx0
   a(0, 1) = -cosph; // d(dca) / dy0

   // a(1,0) = -cosph / tanth; // dz / dx0
   // a(1,1) = -sinph / tanth; // dz / dy0
   a(1, 0) = -cosph * coTan; // dz / dx0
   a(1, 1) = -sinph * coTan; // dz / dy0
   a(1, 2) = 1.0;            // dz / dz0

   // B = dHn / dq
   Double_t x0 = xk0(0, 0);
   Double_t y0 = xk0(1, 0);
   Double_t r  = TMath::Sqrt(x0 * x0 + y0 * y0);
   // Double_t phi1 = TMath::ACos (x0 / r);
   Double_t phi1   = TMath::ATan2(y0, x0);
   Double_t ksi    = part.Phi() - phi1;
   Double_t cosksi = TMath::Cos(ksi);
   Double_t sinksi = TMath::Sin(ksi);

   b       = 0.0;
   b(0, 0) = r * cosksi; // d(dca) / dphi

   // b(1,0) = r * sinksi / tanth; // dz / dphi
   b(1, 0) = r * sinksi * coTan; // dz / dphi
   // b(1,1) = r * cosph / sinth / sinth; // dz / d(theta)
   b(1, 1) = r * cosksi / sinth / sinth; // dz / d(theta)

   b(2, 0) = 1.0; // dphi / dphi

   b(3, 1) = 1.0; // d(theta) / d(theta)

   b(4, 2) = 1.0; // dp / dp
   if (!flag) return;

   // ck0
   ck0(0, 0) = r * sinksi;
   // ck0(1,0) = xk0(2,0) - r * cosksi / tanth;
   ck0(1, 0) = xk0(2, 0) - r * cosksi * coTan;
   ck0(2, 0) = part.GetMeas(2);
   ck0(3, 0) = part.GetMeas(3);
   ck0(4, 0) = part.GetMeas(4);

   TMatrixD qk0(3, 1);
   qk0(0, 0) = part.GetMeas(2);
   qk0(1, 0) = part.GetMeas(3);
   qk0(2, 0) = part.GetMeas(4);
   ck0 -= TMatrixD(a, TMatrixD::kMult, xk0);
   ck0 -= TMatrixD(b, TMatrixD::kMult, qk0);
}

//__________________________________________________________________________
void MpdMotherFitterPart::ParamsAtDca(MpdParticle *part)
{
   /// Compute particle parameters at DCA

   TMatrixD meas(5, 1);

   if (part->GetCharge() == 0) {
      // Neutral
      meas(2, 0)  = part->Phi();
      meas(3, 0)  = part->Theta();
      meas(4, 0)  = part->Momentum();
      Double_t x0 = part->Getx()(0, 0);
      Double_t y0 = part->Getx()(1, 0);
      Double_t r  = TMath::Sqrt(x0 * x0 + y0 * y0);
      // Double_t phi1 = TMath::ACos (x0 / r);
      Double_t phi1 = TMath::ATan2(y0, x0);
      Double_t ksi  = meas(2, 0) - phi1;
      meas(0, 0)    = r * TMath::Sin(ksi);
      // meas(1,0) = part->Getx()(2,0) - r * TMath::Cos(ksi) / TMath::Tan(meas(3,0));
      meas(1, 0) = part->Getx()(2, 0) - r * TMath::Cos(ksi) * TMath::Cos(meas(3, 0)) / TMath::Sin(meas(3, 0));
      part->SetMeas(meas);
      part->SetXY(meas(0, 0) * TMath::Sin(meas(2, 0)), -meas(0, 0) * TMath::Cos(meas(2, 0)));
      return;
   }

   // Charged
   meas(3, 0)     = part->Theta();
   meas(4, 0)     = part->Getq()(2, 0);
   Double_t cosph = TMath::Cos(part->Phi());
   Double_t sinph = TMath::Sin(part->Phi());
   Double_t kinv  = 1. / meas(4, 0);
   Double_t x0 = part->Getx()(0, 0), y0 = part->Getx()(1, 0);
   Double_t expr = x0 * sinph - y0 * cosph - kinv;
   Double_t gam  = TMath::ATan((x0 * cosph + y0 * sinph) / expr);
   meas(2, 0)    = part->Phi() + gam;
   meas(0, 0)    = kinv + expr / TMath::Cos(gam);
   meas(1, 0)    = part->Getx()(2, 0) + gam * kinv * TMath::Cos(meas(3, 0)) / TMath::Sin(meas(3, 0));
   part->SetMeas(meas);
   part->SetXY(-meas(0, 0) * TMath::Cos(meas(2, 0) + TMath::PiOver2()),
               -meas(0, 0) * TMath::Sin(meas(2, 0) + TMath::PiOver2()));
}

//__________________________________________________________________________
void MpdMotherFitterPart::Proxim(const MpdKalmanTrack &track0, MpdKalmanTrack &track)
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
TMatrixD MpdMotherFitterPart::ComputeQmatr(vector<MpdParticle *> vDaught)
{
   /// Compute matrix Q = covariance cov(qk,qj)
   /// Qkj = Wk*Bkt*Gk*Ak*C*Ajt*Gj*Bj*Wj = Wk*Bkt*Gk*Ak*(-E), k!=j
   /// Qii = D

   TMatrixD qc(3, 3), qn(3, 3), qnc(3, 3);
   Int_t    nDaught = vDaught.size();

   for (Int_t i = 0; i < nDaught; ++i) {
      MpdParticle *part = vDaught[i];

      for (Int_t j = 0; j < nDaught; ++j) {
         MpdParticle *part1 = vDaught[j];

         TMatrixD q(3, 3);
         if (j == i) {
            q = part->GetD();
         } else {
            TMatrixD tmp = part1->GetE();
            tmp *= -1.0;
            TMatrixD tmp1 = TMatrixD(part->GetA(), TMatrixD::kMult, tmp);
            TMatrixD tmp2 = TMatrixD(part->GetG(), TMatrixD::kMult, tmp1);
            TMatrixD tmp3 = TMatrixD(part->GetB(), TMatrixD::kTransposeMult, tmp2);
            q             = TMatrixD(part->GetW(), TMatrixD::kMult, tmp3);
         }

         TMatrixD jt    = TMatrixD(TMatrixD::kTransposed, part1->GetJ());
         TMatrixD tmp11 = TMatrixD(q, TMatrixD::kMult, jt);
         TMatrixD tmp12 = TMatrixD(part->GetJ(), TMatrixD::kMult, tmp11);
         if (part->GetCharge() && part1->GetCharge())
            qc += tmp12;
         else if (part->GetCharge() == 0 && part1->GetCharge() == 0)
            qn += tmp12;
         else
            qnc += tmp12;
      }
   }

   qc += qn;
   qc += qnc;
   return qc;
}

//__________________________________________________________________________
Double_t MpdMotherFitterPart::Chi2Vertex(MpdParticle *part, const MpdVertex *vtx)
{
   /// Compute Chi2 w.r.t. vertex

   Double_t    vpos[3] = {vtx->GetX(), vtx->GetY(), vtx->GetZ()};
   TMatrixD    xk(3, 1), xk0(3, 1), ck0(5, 1), a(5, 3), b(5, 3), c(3, 3);
   TMatrixFSym cov(3);
   cov.SetTol(1.e-10); // AZ-311220

   xk0.SetMatrixArray(vpos);
   vtx->CovMatrix(cov);
   cov.Invert();

   ComputeAandB(xk0, *part, a, b, ck0, kTRUE); // compute matrices of derivatives
   TMatrixD g = part->GetG();

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
   TMatrixD m = part->GetMeas();
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
   chi21 += chi22;

   return chi21(0, 0);
}

ClassImp(MpdMotherFitterPart)
