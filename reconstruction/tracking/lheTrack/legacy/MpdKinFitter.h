// -------------------------------------------------------------------------
// -----                     MpdKinFitter header file                  -----
// -----                  Created 11/09/17  by A.Zinchenko             -----
// -----   Original version by N.Geraksiev (formalism by A.Belyaev)    -----
// -------------------------------------------------------------------------

/** MpdKinFitter.h
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Class for kinematic fitting in MPD.
 ** Data level: RECO
 **/

#ifndef MPDKINFITTER_H
#define MPDKINFITTER_H 1

#include "TMatrixD.h"
//#include "TMath.h"
#include <vector>

class MpdParticle;

using namespace std;

class MpdKinFitter : public TObject {

public:
   /** Constructor  **/
   MpdKinFitter(std::vector<TObject *> &daughts); ///< Ctor

   /** Constructor  **/
   MpdKinFitter(std::vector<MpdParticle *> &daughts); ///< Ctor

   /** Destructor **/
   virtual ~MpdKinFitter();

   /** Accessors **/
   Double_t Chi2() const { return fChi2; } ///< Chi2 of mother particle

   // kinfit
   TMatrixD &GetKinx() { return fx; } // pt1,theta1,phi1,pt2,theta2,phi2
   TMatrixD &GetKinVm1() { return fVm1; }
   TMatrixD &GetKinUm1() { return fUm1; } //

   Double_t DoKinFit(Int_t pdg0);  // perform fitting
   Double_t DoKinFit(Double_t m0); // perform fitting

private:
   void Init();
   void CalculateFdFdx(TMatrixD mtr, Double_t m0, Double_t *mD);

   Double_t               fChi2;
   Double_t               fieldConst;
   std::vector<TObject *> fDaughts; // decay products
   TMatrixD               fx;       // 3*nD x 1 parameter matrix pt_i,theta_i,phi_i fitted parameters
   TMatrixD               fF;
   TMatrixD               fdFdx;
   TMatrixD               fdydx;
   TMatrixD               fVm1; // 3*nD x 3*nD covariance matrix of fitted daughter parameters
   TMatrixD               fUm1; // 3x3 covariance matrix of mother parameters: pt theta phi

   ClassDef(MpdKinFitter, 0);
};

#endif
