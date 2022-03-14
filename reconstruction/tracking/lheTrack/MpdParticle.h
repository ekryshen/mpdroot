// -------------------------------------------------------------------------
// -----                      MpdParticle header file                  -----
// -----                  Created 21/01/13  by A. Zinchenko             -----
// -------------------------------------------------------------------------

/** MpdParticle.h
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Class for a particle in MPD (to work with decays).
 ** Data level: RECO
 **/

#ifndef MPDPARTICLE_H
#define MPDPARTICLE_H 1

//#include "TArrayI.h"
#include "TMatrixD.h"
//#include "TMatrixFSym.h"
//#include "TNamed.h"
#include "TMath.h"
#include "TVector3.h"
#include <vector>
class MpdKalmanTrack;
class MpdVertex;

using namespace std;

class MpdParticle : public TObject {

public:
   /** Default constructor  **/
   MpdParticle();                        ///< Default ctor
   MpdParticle(const MpdParticle &part); ///< copy constructor
   MpdParticle(const MpdKalmanTrack &track, Int_t indx = -1, Double_t mass = 0.1396,
               Double_t *orig = NULL); ///< ctor from Kalman track

   MpdParticle &operator=(const MpdParticle &part); ///< assignment operator

   /** Destructor **/
   virtual ~MpdParticle();

   /** Ouput to screen **/
   void Print();

   /** Accessors **/
   Int_t    GetIndx() const { return fIndx; }
   Int_t    GetPdg() const { return fPdg; }
   Int_t    GetCharge() const { return fCharge; }
   Double_t GetMass() const { return fMass; }
   Double_t GetMeas(Int_t i) const { return fMeas(i, 0); }
   Double_t GetXY(Int_t i) const { return fXY0[i]; }
   Int_t    GetFlag() const { return fFlag; }
   // Double_t Phi() const { return GetMeas(2); }
   // Double_t Pt() const { return TMath::Min (TMath::Abs(1./GetMeas(4)*fieldConst), 100.); }
   // Double_t Theta() const { return GetMeas(3); }
   // Double_t Momentum() const { return Pt() / TMath::Sin(Theta()); }
   // TVector3 Momentum3() const { return TVector3(Pt()*TMath::Cos(Phi()), Pt()*TMath::Sin(Phi()),
   //				       Momentum()*TMath::Cos(Theta())); }
   Double_t Phi() const { return fq(0, 0); } // smoothed value
   Double_t Pt() const
   {
      return fCharge == 0 ? TMath::Min(fq(2, 0) * TMath::Sin(Theta()), 100.)
                          : TMath::Min(TMath::Abs(fCharge / fq(2, 0) * fieldConst), 100.);
   }
   // Double_t Theta() const { return fq(1,0); }
   Double_t Theta() const;
   Double_t Momentum() const { return fCharge == 0 ? fq(2, 0) : Pt() / TMath::Sin(Theta()); }
   TVector3 Momentum3() const
   {
      return TVector3(Pt() * TMath::Cos(Phi()), Pt() * TMath::Sin(Phi()), Momentum() * TMath::Cos(Theta()));
   }
   Double_t             Energy() const;
   Double_t             Rapidity() const;
   Double_t             Dca() const { return fMeas(0, 0); } ///< signed DCA
   Int_t                Ndaughters() const { return fDaughtersInds.size(); }
   const vector<Int_t> &DaughterInds() const { return fDaughtersInds; }
   void Track2Part(const MpdKalmanTrack &track, Bool_t setWeight, Double_t *orig); // conversion from track to particle
   const Double_t Chi2Vertex() { return fChi2ver; }                                ///< return Chi2 w.r.t. vertex
   Double_t       Chi2Vertex(MpdVertex *vtx);                                      ///< compute Chi2 w.r.t. vertex
   Double_t       Chi2() const { return fChi2; }                                   ///< Chi2 of mother particle
   Bool_t         Point00() const { return fPoint00; } ///< flag for tracks extrapolated to (0,0)
   void           FillJ();                             // fill Jacobian matrix fJ
   void           FillJinv(TVector3 &mom3);            // fill Jacobian matrix fJinv

   TMatrixD &GetMeas() { return fMeas; }
   TMatrixD &GetJ() { return fJ; }
   TMatrixD &GetJinv() { return fJinv; }
   TMatrixD &GetD() { return fD; }
   TMatrixD &GetE() { return fE; }
   TMatrixD &GetA() { return fA; }
   TMatrixD &GetB() { return fB; }
   TMatrixD &GetC() { return fC; }
   TMatrixD &GetG() { return fG; }
   TMatrixD &GetW() { return fW; }
   TMatrixD &Getq() { return fq; }
   TMatrixD &Getx() { return fx; }
   // TMatrixD& GetW() const { return fW; }

   void SetIndx(Int_t indx) { fIndx = indx; }
   void SetPdg(Int_t pdg)
   {
      fPdg = pdg;
      SetMass();
   }
   void     SetCharge(Int_t charge) { fCharge = charge; }
   void     SetMass(Double_t mass = -2.0);
   void     AddDaughter(Int_t indx) { fDaughtersInds.push_back(indx); }
   Double_t BuildMother(vector<MpdParticle *> &vDaught);
   Double_t BuildMother(vector<MpdKalmanTrack *> &vTracks, vector<MpdParticle *> &vDaught);
   void     SetChi2(Double_t chi2) { fChi2 = chi2; }
   void     SetMeas(TMatrixD &matr) { fMeas = matr; }
   void     SetCovD(TMatrixD &matr) { fD = matr; }
   void     SetCovE(TMatrixD &matr) { fE = matr; }
   // void SetCovQ(TMatrixD &matr) { fQ = matr; }
   void SetA(TMatrixD &matr) { fA = matr; }
   void SetB(TMatrixD &matr) { fB = matr; }
   void SetC(TMatrixD &matr) { fC = matr; }
   void SetG(TMatrixD &matr) { fG = matr; }
   void SetW(TMatrixD &matr) { fW = matr; }
   void Setq(TMatrixD &matr) { fq = matr; }
   void Setx(TMatrixD &matr) { fx = matr; }
   void SetXY(Double_t x, Double_t y)
   {
      fXY0[0] = x;
      fXY0[1] = y;
   }
   void SetFlag(Int_t flag) { fFlag = flag; }

private:
   void WeightAtDca(MpdKalmanTrack &track, Double_t *vert); // obtain MpdParticle weight at DCA

   Int_t         fIndx;          // index of particle
   Int_t         fPdg;           // PDG hypothesis
   Int_t         fCharge;        // charge
   Double_t      fMass;          // particle mass (GeV)
   Double_t      fieldConst;     //! field constant
   vector<Int_t> fDaughtersInds; // indices of particles it is created from
   TMatrixD      fMeas;          // vector of measurements (params)
   TMatrixD      fq;             // geometrical momentum
   TMatrixD      fx;             // particle origin (production vertex)
   Double_t      fXY0[2];        // X and Y at DCA
   TMatrixD      fJ;             //! Jacobian matrix (from geometrical to kinematical momentum)
   TMatrixD      fJinv;          //! Jacobian matrix (from kinematical to geometrical momentum)
   TMatrixD      fD;             //! covariance cov(qk)
   TMatrixD      fE;             //! covariance cov(xk,qk)
   // TMatrixD fQ;                  //! covariance cov(qk,qj)
   TMatrixD fA;       //! derivatives
   TMatrixD fB;       //! derivatives
   TMatrixD fC;       //! covariance cov(xk)
   TMatrixD fG;       // covariance of params
   TMatrixD fW;       //! (Bt*G*B)
   Double_t fChi2;    // Chi2 of mother particle
   Double_t fChi2ver; // Chi2 of particle w.r.t. vertex
   Int_t    fFlag;    // status flag
   Bool_t   fPoint00; // flag for tracks extrapolated to (0,0)

   ClassDef(MpdParticle, 2);
};

//__________________________________________________________________________

inline Double_t MpdParticle::Theta() const
{
   // Theta angle

   if (TMath::Abs(fq(1, 0)) < TMath::PiOver2()) {
      if (TMath::Abs(fq(1, 0)) > 0.001) return fq(1, 0);
      return TMath::Sign(0.001, fq(1, 0));
   } else {
      if (TMath::Abs(fq(1, 0)) < TMath::Pi() - 0.001) return fq(1, 0);
      return TMath::Sign(TMath::Pi() - 0.001, fq(1, 0));
   }
}

//__________________________________________________________________________

#endif
