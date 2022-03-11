//
// MpdMiniBTofPidTraits keeps information about tracks that matched BTOF
//

// C++ headers
#include <limits>
#include <cmath>

// ROOT headers
#include "TMath.h"

// MiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniBTofPidTraits.h"

ClassImp(MpdMiniBTofPidTraits)

   //_________________
   MpdMiniBTofPidTraits::MpdMiniBTofPidTraits()
   : TObject(), fTrackIndex(-1), fHitIndex(-1), fBTofBeta(0), fPx(0), fPy(0), fPz(0)
{
   /* empty */
}

//_________________
MpdMiniBTofPidTraits::MpdMiniBTofPidTraits(const MpdMiniBTofPidTraits &traits) : TObject()
{
   fTrackIndex = traits.fTrackIndex;
   fHitIndex   = traits.fHitIndex;
   fBTofBeta   = traits.fBTofBeta;
   fPx         = traits.fPx;
   fPy         = traits.fPy;
   fPz         = traits.fPz;
}

//_________________
MpdMiniBTofPidTraits::~MpdMiniBTofPidTraits()
{
   /* empty */
}

//_________________
void MpdMiniBTofPidTraits::setBeta(Float_t beta)
{
   if (beta <= 0) {
      fBTofBeta = 0;
   } else {
      fBTofBeta =
         ((beta * 20000.) > std::numeric_limits<unsigned short>::max() ? std::numeric_limits<unsigned short>::max()
                                                                       : (UShort_t)(TMath::Nint(beta * 20000.)));
   }
}

//_________________
void MpdMiniBTofPidTraits::Print(const Char_t *option __attribute__((unused))) const
{
   LOG_INFO << "TrkId: " << fTrackIndex << " HitId: " << fHitIndex << " beta: " << beta() /*<< " p: " << p().Mag()
                                                                             << " m2: " << massSqr() */
            << endm;
}

//_________________
Float_t MpdMiniBTofPidTraits::massSqr() const
{
   return p().Mag2() * (1. / (beta() * beta()) - 1.);
}
