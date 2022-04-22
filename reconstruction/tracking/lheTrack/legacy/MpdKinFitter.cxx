/// \class MpdKinFitter
///
/// MPD kinematic fitter.
///
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)
/// (original version by N.Geraksiev, formalism by A.Belyaev)

#include "MpdKinFitter.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanTrack.h"
//#include "MpdMotherFitterPart.h"
#include "MpdParticle.h"
#include "FairMCTrack.h"

#include <TRandom3.h>
#include <TVector3.h>
#include <TDatabasePDG.h>
//__________________________________________________________________________

MpdKinFitter::MpdKinFitter(vector<TObject *> &daughts) : TObject(), fChi2(0.0), fieldConst(0.0)
{
   /// Constructor

   fDaughts = daughts;
   Init();
}
//__________________________________________________________________________

MpdKinFitter::MpdKinFitter(vector<MpdParticle *> &daughts) : TObject(), fChi2(0.0), fieldConst(0.0)
{
   /// Constructor

   Int_t nDau = daughts.size();
   fDaughts.clear();
   for (Int_t j = 0; j < nDau; ++j) fDaughts.push_back(daughts[j]);
   Init();
}
//__________________________________________________________________________

void MpdKinFitter::Init()
{
   /// Initialization

   Double_t pos[3] = {0.0}, magf[3] = {0.0};
   MpdKalmanFilter::Instance()->GetField(pos, magf);
   fieldConst = 0.3 * 0.01 * magf[2] / 10;

   Int_t nDau = fDaughts.size();
   Int_t npar = nDau * 3;
   fx.ResizeTo(npar, 1);
   fF.ResizeTo(1, 1);
   fdFdx.ResizeTo(npar, 1);
   fdydx.ResizeTo(npar, 3);
   fVm1.ResizeTo(npar, npar);
   fUm1.ResizeTo(3, 3);
}
//__________________________________________________________________________

MpdKinFitter::~MpdKinFitter()
{
   /// Destructor
}
//__________________________________________________________________________

Double_t MpdKinFitter::DoKinFit(Int_t pdg0)
{
   // Kinematic fit with energy conservation constrain (1C fit), minimization with lagrangian multipliers
   // pdg0 - mother PDG code

   Double_t m0 = TDatabasePDG::Instance()->GetParticle(pdg0)->Mass();
   return DoKinFit(m0);
}
//__________________________________________________________________________

Double_t MpdKinFitter::DoKinFit(Double_t m0)
{
   // Kinematic fit with energy conservation constrain (1C fit), minimization with lagrangian multipliers
   // m0 - mother mass
   //
   // parameters: pt, theta, phi
   // TODO only diagonal covariance elements used
   // TODO parameters:phi, theta, fieldConst/pt
   // TODO kalman track doesnt have a pdg code can be set as input parameters
   // TODO setters getters for parameters

   /*
   if (vDaught.size() != 2) {
     cout << "DoKinFit requires two daughters" << endl;
     return -1;
   }
   */
   // constants
   const Int_t    iter_limiter = 10;    // limit iterations, for pt about 10 should be enough
   const Double_t epsLa        = 1e-13; // 1e-13; convergence limiters, different values can be used
   const Double_t epsX         = 1e-14; // 1e-14; TODO use them as input parameters of function or setters?

   Int_t    pdg[9];
   Int_t    defpdg[2] = {2212, 211}; // proton, pion//default values for kalman tracks
   Double_t mD[9];
   // TMatrixD xm(6,1);
   // TMatrixD Gm1(6,6);
   TMatrixD xm(fx);
   TMatrixD Gm1(fVm1);
   // for mc tracks and smearing
   TVector3 mom;
   Double_t gamma = 0.;
   // TMatrixD sigma(6,1);
   TMatrixD sigma(fx);
   TRandom3 r(0);
   Double_t c1, c2, c3, cc1, cc2, cc3; // TODO  different coefficients for parameters?
   c1 = c2 = c3 = 0.005;
   cc1 = cc2 = cc3 = 0.1;

   // Setting parameters based on input classes: MpdParticle, MpdTpcKalmanTrack or FairMCTrack
   for (Int_t i = 0; i < fDaughts.size(); i++) {
      Int_t    j  = 3 * i;       // 3 params per daughter use j for first daughter j = 0,1,2, for second j = 3,4,5
      TObject *od = fDaughts[i]; // object daughter

      if (od->InheritsFrom("MpdParticle")) {
         MpdParticle *d1 = (MpdParticle *)od;
         pdg[i]          = d1->GetPdg();
         xm(j, 0)        = fieldConst / TMath::Abs(d1->Getq()(2, 0)); // pt
         xm(j + 1, 0)    = d1->Getq()(1, 0);                          // theta
         xm(j + 2, 0)    = d1->Getq()(0, 0);                          // phi
         TMatrixD D      = d1->GetD();
         // AZ Gm1(j,j) = D(2,2)*xm(0,0)*xm(0,0)/fieldConst/fieldConst;
         // Gm1(j,j) = D(2,2)*xm(j,0)*xm(j,0)/fieldConst/fieldConst; //AZ ????
         Gm1(j, j)         = D(2, 2) * xm(j, 0) * xm(j, 0) / d1->Getq()(2, 0) / d1->Getq()(2, 0); // AZ
         Gm1(j + 1, j + 1) = D(1, 1);
         Gm1(j + 2, j + 2) = D(0, 0);
      } else if (od->InheritsFrom("MpdKalmanTrack")) {
         MpdKalmanTrack *d1 = (MpdKalmanTrack *)od;
         pdg[i]             = defpdg[i];            // FIXME default pdg code
         xm(j, 0)           = TMath::Abs(d1->Pt()); // pt
         xm(j + 1, 0)       = d1->Theta();
         xm(j + 2, 0)       = d1->Phi();
         TMatrixD D         = *d1->GetCovariance();
         // AZ Gm1(j,j) = D(4,4)*xm(0,0)*xm(0,0);
         Gm1(j, j)         = D(4, 4) * xm(j, 0) * xm(j, 0);
         Gm1(j + 1, j + 1) = D(3, 3);
         Gm1(j + 2, j + 2) = D(2, 2);
      } else if (od->InheritsFrom("FairMCTrack")) { // mc tracks with smearing
         FairMCTrack *d1 = (FairMCTrack *)od;
         pdg[i]          = d1->GetPdgCode();
         d1->GetMomentum(mom);
         // smearing			// get a new random value for each param
         gamma             = r.Gaus(0, 1);
         sigma(j, 0)       = c1 * mom.Perp() * (1 + cc1 * gamma);
         gamma             = r.Gaus(0, 1);
         sigma(j + 1, 0)   = c2 * (1 + cc2 * gamma);
         gamma             = r.Gaus(0, 1);
         sigma(j + 2, 0)   = c3 * (1 + cc3 * gamma);
         gamma             = r.Gaus(0, 1);
         xm(j, 0)          = mom.Perp() + sigma(j, 0) * gamma;
         gamma             = r.Gaus(0, 1);
         xm(j + 1, 0)      = mom.Theta() + sigma(j + 1, 0) * gamma;
         gamma             = r.Gaus(0, 1);
         xm(j + 2, 0)      = mom.Phi() + sigma(j + 2, 0) * gamma;
         Gm1(j, j)         = sigma(j, 0) * sigma(j, 0);
         Gm1(j + 1, j + 1) = sigma(j + 1, 0) * sigma(j + 1, 0);
         ;
         Gm1(j + 2, j + 2) = sigma(j + 2, 0) * sigma(j + 2, 0);
         ;
      }

      mD[i] = TDatabasePDG::Instance()->GetParticle(pdg[i])->Mass();
   } // for n daughters

   if (Gm1(0, 0) == 0) return -2; // vertex fit has failed

   fx = xm;                    // initial value x0
   CalculateFdFdx(fx, m0, mD); // Calc F dFdx with new parameters

   Int_t    niter        = 0;
   Bool_t   twice        = kFALSE;
   Bool_t   convergence  = kFALSE;
   Bool_t   loopPossible = kTRUE;
   TMatrixD H(1, 1);
   TMatrixD Hm1(1, 1); // inverted H
   TMatrixD la(1, 1);

   while (loopPossible) {
      niter++;
      if (niter >= iter_limiter) {
         loopPossible = kFALSE;
         return -1;
      }
      TMatrixD dFdx0 = fdFdx; // save old value of dFdx for psi calculation
      // E = Gm1 * dFdx;
      TMatrixD E(Gm1, TMatrixD::kMult, fdFdx);
      // H = Et*dFdx;
      H             = TMatrixD(E, TMatrixD::kTransposeMult, fdFdx);
      Double_t detH = H.Determinant();
      if (detH == 0) { // has inverse if DetA != 0
         loopPossible = kFALSE;
         return -1;
      }
      Hm1 = TMatrixD(TMatrixD::kInverted, H);
      // b = F + (xm - x).T()*dFdx;
      TMatrixD difx = xm;
      difx -= fx;
      TMatrixD tmp0(difx, TMatrixD::kTransposeMult, fdFdx);
      TMatrixD b = fF;
      b += tmp0;
      TMatrixD bt(TMatrixD::kTransposed, b);
      // cout<<"bt2";bt.Print();
      // la = Hm1*bt;
      la = TMatrixD(Hm1, TMatrixD::kMult, bt);
      // x = xm - E*la;
      TMatrixD tmp1(E, TMatrixD::kMult, la);
      fx = xm;
      fx -= tmp1;
      // Calc F dFdx with new parameters
      CalculateFdFdx(fx, m0, mD);
      // dLla = F*Hm1*Ft;
      TMatrixD Ft(TMatrixD::kTransposed, fF);
      TMatrixD tmp2(Hm1, TMatrixD::kMult, Ft);
      TMatrixD dLla(fF, TMatrixD::kMult, tmp2);
      if (dLla(0, 0) < epsLa) {
         TMatrixD Psi = fdFdx;
         Psi -= dFdx0;
         Psi *= la;
         TMatrixD tmp3(Gm1, TMatrixD::kMult, Psi);
         TMatrixD dLx(Psi, TMatrixD::kTransposeMult, tmp3);
         if (dLx(0, 0) >= epsX) {
            twice = kFALSE;
         } else {
            if (twice) {
               loopPossible = kFALSE;
               convergence  = kTRUE;
               break;
            } else {
               twice = kTRUE;
            }
         }
      } else {
         twice = kFALSE;
      } // epsLa
   }    // while loop

   // chi2 = lat*H*la;
   TMatrixD tmp4(H, TMatrixD::kMult, la);
   TMatrixD kinchi2(la, TMatrixD::kTransposeMult, tmp4);
   if (kinchi2(0, 0) == 0.) return -1;
   // E = Gm1*dFdx;
   TMatrixD E(Gm1, TMatrixD::kMult, fdFdx); // E with last dFdx, recalculate H, Hm1 or use previous?
   // Vm1 = Gm1 - E*Hm1*Et;  //updated covariance matrix of daughter parameters
   TMatrixD tmp5(TMatrixD::kTransposed, E);
   TMatrixD tmp6(Hm1, TMatrixD::kMult, tmp5);
   TMatrixD tmp7(E, TMatrixD::kMult, tmp6);
   fVm1 = Gm1;
   fVm1 -= tmp7;
   // Um1 = dydxt*Vm1*dydx;	 //updated covariance matrix of mother parameters
   TMatrixD tmp8(fVm1, TMatrixD::kMult, fdydx);
   fUm1 = TMatrixD(fdydx, TMatrixD::kTransposeMult, tmp8);

   return kinchi2(0, 0);

   // TODO non-diagonal elements ? convert C/pt to pt?
   /*
      Gm1(0,1) = part1D(2,1); //Cov(C/pt and theta) multiply by xm/fieldConst ?
      Gm1(0,2) = part1D(2,0); //Cov(C/pt phi)
      Gm1(1,0) = part1D(1,2); //Cov(theta C/pt)
      Gm1(1,2) = part1D(1,0); //Cov(theta phi)
      Gm1(2,0) = part1D(0,2); //Cov(phi C/pt)
      Gm1(2,1) = part1D(0,1); //Cov(phi, theta)
      Gm1(3,4) = part2D(2,1);
      Gm1(3,5) = part2D(2,0);
      Gm1(4,3) = part2D(1,2);
      Gm1(4,5) = part2D(1,0);
      Gm1(5,4) = part2D(0,2);
      Gm1(5,4) = part2D(0,1);
   */
}

//__________________________________________________________________________
void MpdKinFitter::CalculateFdFdx(TMatrixD mtr, Double_t m0, Double_t *mD)
{
   // Calculate F and derivatives, 1C - fit

   Int_t nDau = fDaughts.size(), ndim = nDau * 3;
   ;
   Double_t px0 = 0, py0 = 0, pz0 = 0, esum = 0;
   // Matrix containing partial derivatives of momentum x,y,z x dxi = 6x3
   TMatrixD pDp(ndim, 3), deidx(ndim, 1);

   for (Int_t i = 0; i < nDau; ++i) {
      Int_t    j   = i * 3;
      Double_t pt  = mtr(j, 0);
      Double_t th  = mtr(j + 1, 0);
      Double_t phi = mtr(j + 2, 0);

      Double_t costh  = TMath::Cos(th);
      Double_t cosphi = TMath::Cos(phi);
      Double_t sinth  = TMath::Sin(th);
      Double_t sinphi = TMath::Sin(phi);
      Double_t cotth  = costh / sinth;

      // Double_t cosphi1Mphi2 = TMath::Cos(phi1 - phi2);

      Double_t p  = pt / sinth;
      Double_t p2 = p * p;
      // Double_t pt2sin2th = pt * pt / (sinth * sinth);
      Double_t e = TMath::Sqrt(mD[i] * mD[i] + p2);
      esum += e;

      Double_t px = pt * cosphi;
      Double_t py = pt * sinphi;
      Double_t pz = pt * cotth;
      px0 += px;
      py0 += py;
      pz0 += pz;

      pDp(j + 0, 0) = cosphi;
      pDp(j + 0, 1) = sinphi;
      pDp(j + 0, 2) = cotth;
      // pDp(j+1,0) = 0;		    pDp(j+1,1) = 0;		pDp(j+1,2) = - pt12sin2th1/pt1;
      pDp(j + 1, 0) = 0;
      pDp(j + 1, 1) = 0;
      pDp(j + 1, 2) = -p / sinth;
      pDp(j + 2, 0) = -py;
      pDp(j + 2, 1) = px;
      pDp(j + 2, 2) = 0;

      // Matrix containing partial derivatives of e
      // deidx(j+0,0) = pt12sin2th1/(pt1*e1);
      // deidx(j+1,0) = - pt12sin2th1*cotth1/e1;
      deidx(j + 0, 0) = p / (sinth * e);
      deidx(j + 1, 0) = -p2 * cotth / e;
      deidx(j + 2, 0) = 0;
   }

   Double_t e0 = TMath::Sqrt(m0 * m0 + px0 * px0 + py0 * py0 + pz0 * pz0);
   fF(0, 0)    = e0 - esum;
   fdydx       = pDp;

   for (Int_t i = 0; i < ndim; i++)
      fdFdx(i, 0) = (px0 * pDp(i, 0) + py0 * pDp(i, 1) + pz0 * pDp(i, 2)) / e0 - deidx(i, 0);
}
//__________________________________________________________________________

ClassImp(MpdKinFitter);
