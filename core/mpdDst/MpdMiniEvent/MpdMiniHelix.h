/**
 * \class MpdMiniHelix
 * \brief Helix parametrization that uses ROOT TVector3
 *
 * Parametrization of a helix (modification of StHelix). Can also cope
 * with straight tracks, i.e. with zero curvature. This represents only
 * the mathematical model of a helix.
 *
 * \author Grigory Nigmatkulov
 * \email nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniHelix_h
#define MpdMiniHelix_h

// C++ headers
#include <math.h>
#include <utility>
#include <algorithm>

// ROOT headers
#include "TVector3.h"

// MiniDst headers
#include "SystemOfUnits.h"

// Declare C++ namespaces
#if !defined(ST_NO_NAMESPACES)
using std::max;
using std::pair;
using std::swap;
#endif

//_________________
class MpdMiniHelix {

public:
   /// Default constructor
   MpdMiniHelix();

   /// Constructor that takes next arguments:
   /// curvature, dip angle, phase, origin, h
   MpdMiniHelix(Double_t c, Double_t dip, Double_t phase, const TVector3 &o, Int_t h = -1);

   /// Copy constructor
   MpdMiniHelix(const MpdMiniHelix &);

   // Assignment operator (will use the one, provided by compiler)
   // MpdMiniHelix& operator=(const MpdMiniHelix&);

   /// Destructor
   virtual ~MpdMiniHelix();

   /// Return dip angle
   Double_t dipAngle() const;
   /// Return curvature: 1/R in xy-plane
   Double_t curvature() const;
   /// Return phase: aziumth in xy-plane measured from ring center
   Double_t phase() const;
   /// Return x-center of circle in xy-plane
   Double_t xcenter() const;
   /// Return y-center of circle in xy-plane
   Double_t ycenter() const;
   /// Return -sign(q*B);
   Int_t h() const;

   /// Return origin of the helix = starting point
   const TVector3 &origin() const;

   /// Set helix parameters
   void setParameters(Double_t c, Double_t dip, Double_t phase, const TVector3 &o, Int_t h);

   /// coordinates of helix at point s
   Double_t x(Double_t s) const;
   Double_t y(Double_t s) const;
   Double_t z(Double_t s) const;
   TVector3 at(Double_t s) const;

   /// pointing vector of helix at point s
   Double_t cx(Double_t s) const;
   Double_t cy(Double_t s) const;
   Double_t cz(Double_t s = 0) const;
   TVector3 cat(Double_t s) const;

   /// returns period length of helix
   Double_t period() const;

   /// path length at given r (cylindrical r)
   pair<Double_t, Double_t> pathLength(Double_t r) const;

   /// path length at given r (cylindrical r, cylinder axis at x,y)
   pair<Double_t, Double_t> pathLength(Double_t r, Double_t x, Double_t y);

   /// path length at distance of closest approach to a given point
   Double_t pathLength(const TVector3 &p, Bool_t scanPeriods = true) const;

   /// path length at intersection with plane
   Double_t pathLength(const TVector3 &r, const TVector3 &n) const;

   /// path length at distance of closest approach in the xy-plane to a given point
   Double_t pathLength(Double_t x, Double_t y) const;

   /// path lengths at dca between two helices
   pair<Double_t, Double_t> pathLengths(const MpdMiniHelix &, Double_t minStepSize = 10 * micrometer,
                                        Double_t minRange = 10 * centimeter) const;

   /// minimal distance between point and helix
   Double_t distance(const TVector3 &p, Bool_t scanPeriods = true) const;

   /// checks for valid parametrization
   Bool_t valid(Double_t world = 1.e+5) const { return !bad(world); }
   Int_t  bad(Double_t world = 1.e+5) const;

   /// Move the origin along the helix to s which becomes then s=0
   virtual void moveOrigin(Double_t s);

   static const Double_t NoSolution;

protected:
   /// Set curvature of the helix
   void setCurvature(Double_t); /// performs also various checks
   /// Set phase of the helix
   void setPhase(Double_t);
   /// Set dip angle of the helix
   void setDipAngle(Double_t);
   /// Value of S where distance in x-y plane is minimal
   Double_t fudgePathLength(const TVector3 &) const;

protected:
   /// true for straight line case (B=0)
   Bool_t fSingularity;
   /// starting point of a helix
   TVector3 fOrigin;
   /// Dip angle
   Double_t fDipAngle;
   /// Curvature = 1/R
   Double_t fCurvature;
   /// Phase
   Double_t fPhase;
   /// -sign(q*B);
   Int_t fH;

   /// Cos of dip angle
   Double_t fCosDipAngle;
   /// Sin of dip angle
   Double_t fSinDipAngle;
   /// Cos of phase
   Double_t fCosPhase;
   //// Sin of phase
   Double_t fSinPhase;

   ClassDef(MpdMiniHelix, 1)
};

//
//     Non-member functions
//
Int_t         operator==(const MpdMiniHelix &, const MpdMiniHelix &);
Int_t         operator!=(const MpdMiniHelix &, const MpdMiniHelix &);
std::ostream &operator<<(std::ostream &, const MpdMiniHelix &);

//
//     Inline functions
//
inline Int_t MpdMiniHelix::h() const
{
   return fH;
}

inline Double_t MpdMiniHelix::dipAngle() const
{
   return fDipAngle;
}

inline Double_t MpdMiniHelix::curvature() const
{
   return fCurvature;
}

inline Double_t MpdMiniHelix::phase() const
{
   return fPhase;
}

inline Double_t MpdMiniHelix::x(Double_t s) const
{
   if (fSingularity)
      return fOrigin.x() - s * fCosDipAngle * fSinPhase;
   else
      return fOrigin.x() + (cos(fPhase + s * fH * fCurvature * fCosDipAngle) - fCosPhase) / fCurvature;
}

inline Double_t MpdMiniHelix::y(Double_t s) const
{
   if (fSingularity)
      return fOrigin.y() + s * fCosDipAngle * fCosPhase;
   else
      return fOrigin.y() + (sin(fPhase + s * fH * fCurvature * fCosDipAngle) - fSinPhase) / fCurvature;
}

inline Double_t MpdMiniHelix::z(Double_t s) const
{
   return fOrigin.z() + s * fSinDipAngle;
}

inline Double_t MpdMiniHelix::cx(Double_t s) const
{
   if (fSingularity)
      return -fCosDipAngle * fSinPhase;
   else
      return -sin(fPhase + s * fH * fCurvature * fCosDipAngle) * fH * fCosDipAngle;
}

inline Double_t MpdMiniHelix::cy(Double_t s) const
{
   if (fSingularity)
      return fCosDipAngle * fCosPhase;
   else
      return cos(fPhase + s * fH * fCurvature * fCosDipAngle) * fH * fCosDipAngle;
}

inline Double_t MpdMiniHelix::cz(Double_t /* s */) const
{
   return fSinDipAngle;
}
inline const TVector3 &MpdMiniHelix::origin() const
{
   return fOrigin;
}
inline TVector3 MpdMiniHelix::at(Double_t s) const
{
   return TVector3(x(s), y(s), z(s));
}
inline TVector3 MpdMiniHelix::cat(Double_t s) const
{
   return TVector3(cx(s), cy(s), cz(s));
}
inline Double_t MpdMiniHelix::pathLength(Double_t X, Double_t Y) const
{
   return fudgePathLength(TVector3(X, Y, 0));
}

inline Int_t MpdMiniHelix::bad(Double_t WorldSize) const
{

   Int_t ierr;
   if (!::finite(fDipAngle)) {
      return 11;
   }
   if (!::finite(fCurvature)) {
      return 12;
   }

   // ierr = mOrigin.bad(WorldSize);

   // The line above is commented and the StThreeVector::bad(double)
   // is rewritten here
   for (Int_t iIter = 0; iIter < 3; iIter++) {

      Double_t tmpVal;
      // Value StThreeVector.mX1[iter] ???
      switch (iIter) {
      case 0: tmpVal = fOrigin.X(); break;
      case 1: tmpVal = fOrigin.Y(); break;
      case 2: tmpVal = fOrigin.Z(); break;
      default: tmpVal = NAN;
      };

      if (!::finite(tmpVal)) {
         ierr = 10 + iIter;
      }
      if (::fabs(tmpVal) > WorldSize) {
         ierr = 20 + iIter;
      }
   } // for(Int_t iIter=0; iIter<3; iIter+)

   if (ierr) {
      return (3 + ierr * 100);
   }
   if (::fabs(fDipAngle) > 1.58) {
      return 21;
   }

   Double_t qwe = ::fabs(::fabs(fDipAngle) - M_PI / 2);
   if (qwe < 1. / WorldSize) {
      return 31;
   }

   if (::fabs(fCurvature) > WorldSize) {
      return 22;
   }
   if (fCurvature < 0) {
      return 32;
   }
   if (abs(fH) != 1) {
      return 24;
   }

   return 0;
}

#endif // #define MpdMiniHelix_h
