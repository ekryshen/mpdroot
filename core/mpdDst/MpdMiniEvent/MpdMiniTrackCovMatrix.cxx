//
// Holds information about track covariance matrix
//

// MiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniTrackCovMatrix.h"

ClassImp(MpdMiniTrackCovMatrix);

   //_________________
   MpdMiniTrackCovMatrix::MpdMiniTrackCovMatrix()
   : mSigma{}, mCorr{}
{
   rphi0        = 0;
   z            = 0.;
   phi          = 0;
   lambda       = 0.;
   minusQoverPt = 0.;
   r            = 0.;
   r0           = 0.;
   /* empty */
}

//_________________

MpdMiniTrackCovMatrix::MpdMiniTrackCovMatrix(const MpdMiniTrackCovMatrix &mtx) : TObject()
{
   for (Int_t iIter = 0; iIter < 5; iIter++) mSigma[iIter] = mtx.mSigma[iIter];

   for (Int_t iIter = 0; iIter < 10; iIter++) mCorr[iIter] = mtx.mCorr[iIter];
}

//_________________

MpdMiniTrackCovMatrix::~MpdMiniTrackCovMatrix()
{
   /* empty */
}

//_________________

void MpdMiniTrackCovMatrix::Print(Char_t const *option __attribute__((unused))) const
{
   const Float_t *lSigma = sigmas();
   const Float_t *lCorr  = correlations();
   LOG_INFO << "sigmas: \n"
            << lSigma[0] << "/" << lSigma[1] << "/" << lSigma[2] << "/" << lSigma[3] << "/" << lSigma[4] << "\n"
            << "correlations: \n"
            << lCorr[0] << "/" << lCorr[1] << "/" << lCorr[2] << "/" << lCorr[3] << "/" << lCorr[4] << "/" << lCorr[5]
            << "/" << lCorr[6] << "/" << lCorr[7] << "/" << lCorr[8] << "/" << lCorr[9] << endm;
}

//_________________

Bool_t MpdMiniTrackCovMatrix::isBadCovMatrix()
{
   return (mSigma[0] == 0 && mSigma[1] == 0 && mSigma[2] == 0 && mSigma[3] == 0 && mSigma[4] == 0 && mCorr[0] == 0 &&
           mCorr[1] == 0 && mCorr[2] == 0 && mCorr[3] == 0 && mCorr[4] == 0 && mCorr[5] == 0 && mCorr[6] == 0 &&
           mCorr[7] == 0 && mCorr[8] == 0 && mCorr[9] == 0);
}
