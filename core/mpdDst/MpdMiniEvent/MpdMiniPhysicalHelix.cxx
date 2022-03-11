//
// MpdMiniPhysicalHelix is a parametrization of a particle that moves along the helix
//

// C++ headers
#include <math.h>

// MiniDst headers
#include "MpdMiniHelix.h"
#include "MpdMiniPhysicalHelix.h"
#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

ClassImpT(MpdMiniPhysicalHelix, Double_t);

//_________________
MpdMiniPhysicalHelix::MpdMiniPhysicalHelix()
{
   /* no-op */
}

//_________________
MpdMiniPhysicalHelix::~MpdMiniPhysicalHelix()
{
   /* no-op */
}

//_________________
MpdMiniPhysicalHelix::MpdMiniPhysicalHelix(const TVector3 &p, const TVector3 &o, Double_t B, Double_t q)
{
   fH = (q * B <= 0) ? 1 : -1;
   if (p.y() == 0 && p.x() == 0) {
      setPhase((M_PI / 4) * (1 - 2. * fH));
   } else {
      setPhase(atan2(p.y(), p.x()) - fH * M_PI / 2);
   }
   setDipAngle(atan2(p.z(), p.Perp()));
   fOrigin = o;

#ifndef ST_NO_NAMESPACES
   {
      using namespace units;
#endif
      setCurvature(::fabs((c_light * nanosecond / meter * q * B / tesla) / (p.Mag() / GeV * fCosDipAngle) / meter));
#ifndef ST_NO_NAMESPACES
   }
#endif
}

//_________________
MpdMiniPhysicalHelix::MpdMiniPhysicalHelix(Double_t c, Double_t d, Double_t phase, const TVector3 &o, Int_t h)
   : MpdMiniHelix(c, d, phase, o, h)
{
   /* no-op */
}

//_________________
TVector3 MpdMiniPhysicalHelix::momentum(Double_t B) const
{

   if (fSingularity) {
      return (TVector3(0, 0, 0));
   } else {
#ifndef ST_NO_NAMESPACES
      {
         using namespace units;
#endif
         Double_t pt = GeV * fabs(c_light * nanosecond / meter * B / tesla) / (fabs(fCurvature) * meter);

         return (TVector3(pt * cos(fPhase + fH * M_PI / 2), // pos part pos field
                          pt * sin(fPhase + fH * M_PI / 2), pt * tan(fDipAngle)));
#ifndef ST_NO_NAMESPACES
      }
#endif
   } // else
}

//_________________
TVector3 MpdMiniPhysicalHelix::momentumAt(Double_t S, Double_t B) const
{
   // Obtain phase-shifted momentum from phase-shift of origin
   MpdMiniPhysicalHelix tmp(*this);
   tmp.moveOrigin(S);
   return tmp.momentum(B);
}

//_________________
Int_t MpdMiniPhysicalHelix::charge(Double_t B) const
{
   return (B > 0 ? -fH : fH);
}

//_________________
Double_t MpdMiniPhysicalHelix::geometricSignedDistance(Double_t x, Double_t y)
{
   // Geometric signed distance
   Double_t thePath       = this->pathLength(x, y);
   TVector3 DCA2dPosition = this->at(thePath);
   DCA2dPosition.SetZ(0);
   TVector3 position(x, y, 0);
   TVector3 DCAVec = (DCA2dPosition - position);
   TVector3 momVec;
   // Deal with straight tracks
   if (this->fSingularity) {
      momVec = this->at(1) - this->at(0);
      momVec.SetZ(0);
   } else {
      momVec = this->momentumAt(thePath, 1. / tesla); // Don't care about Bmag.  Helicity is what matters.
      momVec.SetZ(0);
   }

   Double_t cross   = DCAVec.x() * momVec.y() - DCAVec.y() * momVec.x();
   Double_t theSign = (cross >= 0) ? 1. : -1.;
   return theSign * DCAVec.Perp();
}

//_________________
Double_t MpdMiniPhysicalHelix::curvatureSignedDistance(Double_t x, Double_t y)
{
   // Protect against mH = 0 or zero field
   if (this->fSingularity || abs(this->fH) <= 0) {
      return (this->geometricSignedDistance(x, y));
   } else {
      return (this->geometricSignedDistance(x, y)) / (this->fH);
   }
}

//_________________
Double_t MpdMiniPhysicalHelix::geometricSignedDistance(const TVector3 &pos)
{
   Double_t sdca2d  = this->geometricSignedDistance(pos.x(), pos.y());
   Double_t theSign = (sdca2d >= 0) ? 1. : -1.;
   return (this->distance(pos)) * theSign;
}

//_________________
Double_t MpdMiniPhysicalHelix::curvatureSignedDistance(const TVector3 &pos)
{
   Double_t sdca2d  = this->curvatureSignedDistance(pos.x(), pos.y());
   Double_t theSign = (sdca2d >= 0) ? 1. : -1.;
   return (this->distance(pos)) * theSign;
}
