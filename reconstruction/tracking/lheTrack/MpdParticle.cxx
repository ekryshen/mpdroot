/// \class MpdParticle
///
/// Particle object in MPD (to work with decays).
///
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdParticle.h"
#include "MpdHelix.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
#include "MpdMotherFitterPart.h"

#include <TDatabasePDG.h>

//__________________________________________________________________________
MpdParticle::MpdParticle()
   : TObject(), fIndx(-1), fCharge(0), fMass(TDatabasePDG::Instance()->GetParticle("pi+")->Mass()), fChi2(0.0),
     fChi2ver(-1.0), fFlag(0), fPoint00(kTRUE)
{
   /// Default constructor

   fMeas.ResizeTo(5, 1);
   fq.ResizeTo(3, 1);
   fx.ResizeTo(3, 1);
   fJ.ResizeTo(3, 3);
   fJinv.ResizeTo(3, 3);
   fD.ResizeTo(3, 3);
   fE.ResizeTo(3, 3);
   // fQ.ResizeTo(3,3);
   fA.ResizeTo(5, 3);
   fB.ResizeTo(5, 3);
   fC.ResizeTo(3, 3);
   fG.ResizeTo(5, 5);
   fW.ResizeTo(3, 3);

   fDaughtersInds.clear();
   Double_t pos[3] = {0.0}, magf[3] = {0.0};
   MpdKalmanFilter::Instance()->GetField(pos, magf);
   fieldConst = 0.3 * 0.01 * magf[2] / 10;
   fXY0[0] = fXY0[1] = 0.0;
}

//__________________________________________________________________________
MpdParticle::MpdParticle(const MpdKalmanTrack &track, Int_t indx, Double_t mass, Double_t *orig)
   : TObject(track), fIndx(indx), fCharge(track.Charge()), fMass(mass), fChi2(0.0), fChi2ver(-1.0), fFlag(0),
     fPoint00(kTRUE)
{
   /// Constructor from Kalman track

   fMeas.ResizeTo(5, 1);
   fq.ResizeTo(3, 1);
   fx.ResizeTo(3, 1);
   fJ.ResizeTo(3, 3);
   fJinv.ResizeTo(3, 3);
   fD.ResizeTo(3, 3);
   fE.ResizeTo(3, 3);
   // fQ.ResizeTo(3,3);
   fA.ResizeTo(5, 3);
   fB.ResizeTo(5, 3);
   fC.ResizeTo(3, 3);
   fG.ResizeTo(5, 5);
   fW.ResizeTo(3, 3);

   fDaughtersInds.clear();
   AddDaughter(indx);
   Double_t pos[3] = {0.0}, magf[3] = {0.0};
   MpdKalmanFilter::Instance()->GetField(pos, magf);
   fieldConst = 0.3 * 0.01 * magf[2] / 10;
   fXY0[0] = fXY0[1] = 0.0;
   // AZ if (orig) fPoint00 = kFALSE;

   Track2Part(track, kTRUE, orig);
}

//__________________________________________________________________________
MpdParticle::MpdParticle(const MpdParticle &part)
   : TObject(part), fIndx(part.fIndx), fPdg(part.fPdg), fCharge(part.fCharge), fMass(part.fMass),
     fieldConst(part.fieldConst), fMeas(part.fMeas), fq(part.fq), fx(part.fx), fJ(part.fJ), fJinv(part.fJinv),
     fD(part.fD), fE(part.fE), fA(part.fA), fB(part.fB), fC(part.fC), fG(part.fG), fW(part.fW), fChi2(part.fChi2),
     fChi2ver(part.fChi2ver), fFlag(part.fFlag), fPoint00(part.fPoint00)
{
   /// copy constructor

   fXY0[0]        = part.fXY0[0];
   fXY0[1]        = part.fXY0[1];
   fDaughtersInds = part.fDaughtersInds;
}

//__________________________________________________________________________
MpdParticle &MpdParticle::operator=(const MpdParticle &part)
{
   /// Asignment operator

   // check assignement to self
   if (this == &part) return *this;

   // base class assignement
   TObject::operator=(part);

   fIndx      = part.fIndx;
   fPdg       = part.fPdg;
   fCharge    = part.fCharge;
   fMass      = part.fMass;
   fieldConst = part.fieldConst;
   fMeas      = part.fMeas;
   fq         = part.fq;
   fx         = part.fx;
   fJ         = part.fJ;
   fJinv      = part.fJinv;
   fD         = part.fD;
   fE         = part.fE;
   fA         = part.fA;
   fB         = part.fB;
   fC         = part.fC;
   fG         = part.fG;
   fW         = part.fW;
   fChi2      = part.fChi2;
   fChi2ver   = part.fChi2ver;
   fFlag      = part.fFlag;

   fXY0[0]        = part.fXY0[0];
   fXY0[1]        = part.fXY0[1];
   fDaughtersInds = part.fDaughtersInds;

   return *this;
}

//__________________________________________________________________________
MpdParticle::~MpdParticle()
{
   /// Destructor
}

/*
//__________________________________________________________________________
Int_t MpdParticle::Compare(const TObject* part) const
{
/// "Compare" function to sort in descending order in pt

  MpdKalmanTrack *trackKF = (MpdKalmanTrack*) track;
  Double_t pt = 1. / TMath::Abs(trackKF->GetParam(4));
  Double_t ptthis = 1. / TMath::Abs((*fParam)(4,0));
  if (ptthis < pt) return 1;
  else if (ptthis > pt) return -1;
  return 0;
}
*/
//__________________________________________________________________________
void MpdParticle::SetMass(Double_t mass)
{
   /// Set particle mass (if negative, use PDG table value)

   fMass = (mass > -1.0) ? mass : TDatabasePDG::Instance()->GetParticle(fPdg)->Mass();
}

//__________________________________________________________________________
void MpdParticle::Track2Part(const MpdKalmanTrack &track, Bool_t setWeight, Double_t *orig)
{
   /// Different parameterization for MpdParticle as compared to MpdKalmanTrack
   /// (see R.Luchsinger, Ch.Grab CPC 76 (1993) 263-280), but also different
   /// from the described in the paper - different order of parameters in order to be
   /// more consistent with MpdKalmanTrack

   // !!! Check for modular tracking (cases with local coordinates (nodes)) !!!

   Double_t vert[3] = {0};
   if (orig == NULL) orig = vert;

   fMeas(4, 0)  = track.GetParam(4) * fieldConst; // k
   Double_t phi = track.GetParam(2);
   fMeas(2, 0)  = TMath::ATan2(TMath::Sin(phi), TMath::Cos(phi)); // phi: -pi:pi
   fMeas(3, 0)  = track.Theta();                                  // theta
   fMeas(0, 0)  = track.GetPosNew();                              // DCA w.r.t. (0,0)
   fMeas(1, 0)  = track.GetParam(1);                              // z at DCA
   if (fMeas(0, 0) > 1.e-7) {
      // Compute sign of DCA (DCA x Pt)
      phi         = track.GetParam(0) / fMeas(0, 0);
      fXY0[0]     = fMeas(0, 0) * TMath::Cos(phi);
      fXY0[1]     = fMeas(0, 0) * TMath::Sin(phi);
      Double_t dx = fXY0[0] - orig[0];
      Double_t dy = fXY0[1] - orig[1];
      TVector3 dcaV(dx, dy, 0.0);
      TVector3 momV(TMath::Cos(fMeas(2, 0)), TMath::Sin(fMeas(2, 0)), 0.0);
      fMeas(0, 0) = TMath::Sqrt(dx * dx + dy * dy) * TMath::Sign(1.0, dcaV.Cross(momV).Z()); // signed DCA w.r.t. "orig"
      // AZ
      fXY0[0] = dx;
      fXY0[1] = dy;
   }
   // cout << fXY0[0] << " " << -fMeas(0,0)*TMath::Cos(fMeas(2,0)+TMath::PiOver2()) << " "
   //    << fXY0[1] << " " << -fMeas(0,0)*TMath::Sin(fMeas(2,0)+TMath::PiOver2()) << endl;

   if (!setWeight) return; // skip all the rest

   fq(0, 0) = fMeas(2, 0);
   fq(1, 0) = fMeas(3, 0);
   fq(2, 0) = fMeas(4, 0);

   // Obtain weight matrix at DCA
   MpdKalmanTrack track1 = track;
   if (fPoint00) {
      // For tracks extrapolated to the point (0,0)
      track1.SetParamNew(*track1.GetParam());
      track1.SetPos(track1.GetPosNew());
      track1.ReSetWeight();
      TMatrixD g = *track1.GetWeight(); // track weight matrix
      // cout << " track Weight " << endl;
      // g.Print();

      // Propagate to the point where track covar. matrix was stored
      MpdKalmanHit hit;
      if (track.GetNode() == "") {
         hit.SetType(MpdKalmanHit::kFixedR);
         hit.SetPos(track.GetPos());
      } else {
         // cout << " !!! Not implemented for local coordinates yet !!! " << endl;
         // exit(0);
         hit.SetType(MpdKalmanHit::kFixedP);
         TString detName = track.GetNode();
         if (track.GetUniqueID()) {
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
      track1.SetWeight(g);

      /*
        TMatrixD jacob = MpdKalmanFilter::Instance()->GetJacob(); // propagation Jacobian
        jacob.Invert();

        // Weight propagation to DCA
        TMatrixD tmp(g,TMatrixD::kMult,jacob); // WD
        TMatrixD weight1(jacob,TMatrixD::kTransposeMult,tmp); // DtWD
        g.Print();
        weight1.Print();
        // Apply Jacobian
        */
   } // if (fPoint00)

   track1.SetPos(track1.GetPosNew());
   // MpdMotherFitterPart::Instance()->WeightAtDca(this, track1, orig);
   WeightAtDca(track1, orig);
   // cout << " Weight at DCA " << endl;
   // fG.Print();
}

//__________________________________________________________________________

void MpdParticle::WeightAtDca(MpdKalmanTrack &tr, Double_t *vert)
{
   // Obtain MpdParticle weight at DCA from MpdKalmanTrack weight

   TMatrixDSym *covar  = tr.Weight2Cov(); // get covariance matrix
   TMatrixD     g      = *tr.GetWeight(); // track weight matrix
   TMatrixD     meas0  = GetMeas();
   TMatrixD     param0 = *tr.GetParamNew();
   if (tr.GetNodeNew() != "") {
      // Local coordinates
      tr.SetNode(tr.GetNodeNew());
      tr.SetNodeNew("");
   }
   TMatrixD jacob(5, 5);
   // Double_t vert[3] = {0.0};

   Double_t dPar;
   // Loop over parameters to find change of the propagated vs initial ones
   Bool_t ok = kTRUE;
   for (Int_t i = 0; i < 5; ++i) {
      dPar = TMath::Sqrt((*covar)(i, i));
      if (i < 4)
         dPar = TMath::Min(dPar, 0.1);
      else
         dPar = TMath::Min(dPar, 0.1 * TMath::Abs(param0(4, 0)));
      tr.SetParam(param0);
      if (i == 4) dPar *= TMath::Sign(1., -param0(4, 0)); // 1/p
      // else if (i == 2) dPar *= sign;
      else if (i == 3)
         dPar *= TMath::Sign(1., -param0(3, 0)); // dip-angle
      tr.SetParam(i, param0(i, 0) + dPar);
      ok = MpdKalmanFilter::Instance()->FindPca(&tr, vert);
      tr.SetParam(*tr.GetParamNew());
      MpdParticle tmpPart;
      tmpPart.Track2Part(tr, kFALSE, vert);

      Double_t meas = 0.0;
      for (Int_t j = 0; j < 5; ++j) {
         if (j == 2)
            meas = MpdKalmanFilter::Instance()->Proxim(meas0(j, 0), tmpPart.GetMeas(j));
         else
            meas = tmpPart.GetMeas(j);
         jacob(j, i) = (meas - meas0(j, 0)) / dPar;
         // cout << i << " " << j << " " << dPar << " " << track->GetParamNew(j) << " " << paramNew0(j,0)
         //    << " " << track->GetPos() << " " << track->GetPosNew() << " " << jacob(j,i) << endl;
      }
   }

   // jacob.Print();
   jacob.Invert();

   TMatrixD tmp(g, TMatrixD::kMult, jacob);                // WD
   TMatrixD weight1(jacob, TMatrixD::kTransposeMult, tmp); // DtWD
   // part->SetG (weight1);
   /*AZ
   if (part->Point00()) {
     // Add multiple scattering
     weight1.Invert(); // covar
     Double_t step = tr.GetPos(); // radius of TPC inner layer
     Double_t x0 = 13363.6; // rad. length - TPCMixture
     TString mass2;
     mass2 += (part->GetMass() * part->GetMass());
     Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(&tr, x0, step, mass2);
     Double_t th = tr.GetParamNew(3);
     Double_t cosTh = TMath::Cos(th);
     weight1(2,2) += (angle2 / cosTh / cosTh);
     weight1(3,3) += angle2;
     weight1.Invert(); // weight
   }
   */
   SetG(weight1);
}

//__________________________________________________________________________
Double_t MpdParticle::BuildMother(vector<MpdParticle *> &vDaught)
{
   /// Build mother particle from daughters which were smoothed
   /// according to the decay vertex constraint (after FindVertex)

   return MpdMotherFitterPart::Instance()->BuildMother(this, vDaught);
}

//__________________________________________________________________________

Double_t MpdParticle::BuildMother(vector<MpdKalmanTrack *> &vTracks, vector<MpdParticle *> &vDaught)
{
   /// Build mother particle from daughters which were smoothed
   /// according to the decay vertex constraint (after FindVertex).
   /// Daughters are built from tracks and parametrized at their
   /// intersection point.

   return MpdMotherFitterPart::Instance()->BuildMother(this, vTracks, vDaught);
}

//__________________________________________________________________________

void MpdParticle::FillJ()
{
   /// Fill Jacobian matrix fJ

   Double_t p0    = Momentum();
   Double_t pt    = Pt();
   Double_t ph0   = Phi();
   Double_t th0   = Theta();
   Double_t k0    = fq(2, 0);
   Double_t cosPh = TMath::Cos(ph0);
   Double_t sinPh = TMath::Sin(ph0);
   Double_t sinTh = TMath::Sin(th0);
   Double_t cosTh = TMath::Cos(th0);
   fJ             = 0.0;

   if (fCharge) {
      // Charged particle
      // d/d(phi)
      fJ(0, 0) = -pt * sinPh;
      fJ(1, 0) = pt * cosPh;
      // d/d(theta)
      fJ(2, 1) = -pt / sinTh / sinTh;
      // d/dk
      fJ(0, 2) = -pt / k0 * cosPh;
      fJ(1, 2) = -pt / k0 * sinPh;
      // fJ(2,2) = -pt / k0 / TMath::Tan(th0);
      Double_t coTan = cosTh / sinTh;
      fJ(2, 2)       = -pt / k0 * coTan;
   } else {
      // neutral - check !!!
      // d/d(phi)
      fJ(0, 0) = -p0 * sinTh * sinPh;
      fJ(1, 0) = p0 * sinTh * cosPh;
      // d/d(theta)
      fJ(0, 1) = p0 * cosTh * cosPh;
      fJ(1, 1) = p0 * cosTh * sinPh;
      fJ(2, 1) = -p0 * sinTh;
      // d/dp
      fJ(0, 2) = sinTh * cosPh;
      fJ(1, 2) = sinTh * sinPh;
      fJ(2, 2) = cosTh;
   }

   TMatrixD tmp(TMatrixD::kTransposed, fJ);
   // fJ = tmp;
}

//__________________________________________________________________________
void MpdParticle::FillJinv(TVector3 &mom3)
{
   /// Fill Jacobian matrix fJinv

   Double_t px0 = mom3.X();
   Double_t py0 = mom3.Y();
   Double_t pz0 = mom3.Z();
   Double_t p0  = mom3.Mag();
   Double_t pt0 = mom3.Pt();
   // Double_t ph0 = mom3.Phi();
   // Double_t th0 = mom3.Theta();
   // Double_t k0 = Theta();
   fJinv = 0.0;

   if (fCharge) {
      // Charged mother - check !!!
      // d/d(px)
      fJinv(0, 0) = -py0 / pt0 / pt0;
      fJinv(1, 0) = px0 * pz0 / p0 / p0 / pt0;
      fJinv(2, 0) = fCharge * fieldConst * px0 / pt0 / pt0 / pt0;
      // d/d(py)
      fJinv(0, 1) = px0 / pt0 / pt0;
      fJinv(1, 1) = py0 * pz0 / p0 / p0 / pt0;
      fJinv(2, 1) = fCharge * fieldConst * py0 / pt0 / pt0 / pt0;
      // d/d(pz)
      fJinv(1, 2) = -pt0 / p0 / p0;
   } else {
      // neutral mother
      // d/d(px)
      fJinv(0, 0) = -py0 / pt0 / pt0;
      fJinv(1, 0) = px0 * pz0 / p0 / p0 / pt0;
      fJinv(2, 0) = px0 / p0;
      // d/d(py)
      fJinv(0, 1) = px0 / pt0 / pt0;
      fJinv(1, 1) = py0 * pz0 / p0 / p0 / pt0;
      fJinv(2, 1) = py0 / p0;
      // d/d(pz)
      fJinv(1, 2) = -pt0 / p0 / p0;
      fJinv(2, 2) = pz0 / p0;
   }

   TMatrixD tmp(TMatrixD::kTransposed, fJinv);
   // fJinv = tmp;
}

//__________________________________________________________________________
Double_t MpdParticle::Chi2Vertex(MpdVertex *vtx)
{
   /// Compute Chi2 w.r.t. vertex

   fChi2ver = MpdMotherFitterPart::Instance()->Chi2Vertex(this, vtx);
   return fChi2ver;
}

//__________________________________________________________________________
Double_t MpdParticle::Energy() const
{
   /// Return particle energy

   Double_t mom = Momentum();
   return TMath::Sqrt(mom * mom + fMass * fMass);
}

//__________________________________________________________________________
Double_t MpdParticle::Rapidity() const
{
   /// Return particle rapidity

   TVector3 mom3 = Momentum3();
   Double_t mom  = Momentum();
   Double_t e    = TMath::Sqrt(mom * mom + fMass * fMass);
   Double_t pz   = mom3.Z();
   Double_t y    = 0.5 * TMath::Log((e + pz) / (e - pz));
   return y;
}

//__________________________________________________________________________
void MpdParticle::Print()
{
   /// Print particle info
}

ClassImp(MpdParticle)