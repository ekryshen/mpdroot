/**
 * \class MpdMiniTrackCovMatrix
 * \brief Stores track covariance matrix
 *
 * The class holds parameters of the covariance matrix
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * email: nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniTrackCovMatrix_h
#define MpdMiniTrackCovMatrix_h

// ROOT headers
#include <TObject.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <vector>

//_________________

class MpdMiniTrackCovMatrix : public TObject {
public:
   /// Default constructor
   MpdMiniTrackCovMatrix();
   /// Copy constructor
   MpdMiniTrackCovMatrix(const MpdMiniTrackCovMatrix &matrix);
   /// Destructor
   virtual ~MpdMiniTrackCovMatrix();
   /// Print option
   virtual void Print(Char_t const *option = "") const;

   //
   // Getters
   //

   /// Return components of state vector

   Float_t rPhi0() { return rphi0; }

   Float_t Z() { return z; }

   Float_t Phi() { return phi; }

   Float_t Lambda() { return lambda; }

   Float_t MinusQoverPt() { return minusQoverPt; }

   Float_t R() { return r; }

   Float_t R0() { return r0; }

   /// Return pointer to the sigma array

   const Float_t *sigmas() const { return mSigma; }
   /// Return pointer to the correlation array

   const Float_t *correlations() const { return mCorr; }

   /// Return true, if all values 0. It corresponds to
   /// the case when track did not have a covariance
   /// matrix in MiniDst
   Bool_t isBadCovMatrix();

   //
   // Setters
   //

   /// Set state vector

   void SetStateVector(Float_t ele0, Float_t ele1, Float_t ele2, Float_t ele3, Float_t ele4)
   {
      rphi0        = ele0;
      z            = ele1;
      phi          = ele2;
      lambda       = ele3;
      minusQoverPt = ele4;
   }

   /// Set 5 sigma parameters

   void setSigmas(Float_t sigmas[5])
   {
      for (Int_t iIter = 0; iIter < 5; iIter++) {
         mSigma[iIter] = sigmas[iIter];
      }
   }
   /// Set 5 sigma parameters

   void setSigmas(std::vector<Double_t> sigmas)
   {
      Int_t iIter = 0;
      for (auto it : sigmas) {
         mSigma[iIter] = it;
         iIter++;
      }
   }
   /// Set 10 correlation parameters

   void setCorrelations(Float_t corr[10])
   {
      for (Int_t iIter = 0; iIter < 10; iIter++) {
         mCorr[iIter] = corr[iIter];
      }
   }
   /// Set 10 correlation parameters

   void setCorrelations(std::vector<Double_t> corr)
   {
      Int_t iIter = 0;
      for (auto it : corr) {
         mCorr[iIter] = it;
         iIter++;
      }
   }

   void SetR(Float_t R) { r = R; }

   void SetR0(Float_t R0) { r0 = R0; }

   /// Formers for state vector and covariance matrix (partially, to be used in KFParticle)

   TMatrixD stateVector()
   {
      const Double_t data[5] = {rphi0, z, phi, lambda, minusQoverPt};
      TMatrixD       matrix(5, 1, data);
      return matrix;
   }

   TMatrixDSym covarianceMatrix()
   {
      const Double_t data[25] = {mSigma[0], mCorr[0], mCorr[1], mCorr[3], mCorr[6],  mCorr[0],  mSigma[1],
                                 mCorr[2],  mCorr[4], mCorr[7], mCorr[1], mCorr[2],  mSigma[2], mCorr[5],
                                 mCorr[8],  mCorr[3], mCorr[4], mCorr[5], mSigma[3], mCorr[9],  mCorr[6],
                                 mCorr[7],  mCorr[8], mCorr[9], mSigma[4]};

      TMatrixDSym matrix(5, data);

      return matrix;
   }

private:
   // Track state vector and covariance matrix used in MPD look as follows:

   // State vector: (r\phi0, Z, \phi, \lambda, -q / Pt)^{T}

   // Covariance matrix: j    0     1     2     3     4 Cov(r\phi0, r\phi0) i 0  0(0) Cov(r\phi0, Z)   Cov(Z,Z) 1 1(0)
   // 2(1) Cov(r\phi0, \phi) Cov(Z, \phi) Cov(\phi, \phi) 2    3(1)  4(2)  5(2) Cov(r\phi0, \lambda) Cov(Z, \lambda)
   // Cov(\phi, \lambda) Cov(\lambda, \lambda)                                       3    6(3)  7(4)  8(5)  9(3)
   // Cov(r\phi0, (-q / Pt)) Cov(Z, (-q / Pt)) Cov(\phi, (-q / Pt)) Cov (\lambda, (-q / Pt)) Cov((-q / Pt), (-q / Pt))
   // 4    10(6) 11(7) 12(8) 13(9) 14(4)

   /* State vector: */
   Float_t rphi0;
   Float_t z;
   Float_t phi;
   Float_t lambda;
   Float_t minusQoverPt;

   /*Propagated radial position */
   Float_t r;
   Float_t r0;

   /* Covariance matrix: */
   /// Diagonal elements
   Float_t mSigma[5];
   /// Off-diagonal elements
   Float_t mCorr[10];

   ClassDef(MpdMiniTrackCovMatrix, 1);
};

#endif // #define MpdMiniTrackCovMatrix_h
