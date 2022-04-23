/**
 * \class MpdFemtoHelix
 * \author Grigory Nigmatkulov, May 07 2018
 *
 * Parametrization of a helix (modification of StHelix). Can also cope
 * with straight tracks, i.e. with zero curvature. This represents only
 * the mathematical model of a helix. See the SCL user guide for more details.
 *
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * email: nigmatkulov@gmail.com
 */

#ifndef MpdFemtoHelix_h
#define MpdFemtoHelix_h

/// C++ headers
#include <math.h>
#include <utility>
#include <algorithm>

/// ROOT headers
#include "TVector3.h"

/// PicoDst headers
//#ifdef _VANILLA_ROOT_
#include "SystemOfUnits.h"
//#else
//#include "StarClassLibrary/SystemOfUnits.h"
//#endif

/// Declare C++ namespaces
#if !defined(ST_NO_NAMESPACES)
using std::max;
using std::pair;
using std::swap;
#endif

//_________________
class MpdFemtoHelix {
public:
   /// Empty constructor
   MpdFemtoHelix();

   /// Constructor that takes next arguments:
   /// curvature, dip angle, phase, origin, h
   MpdFemtoHelix(double c, double dip, double phase, const TVector3 &o, int h = -1);

   /// Copy constructor
   MpdFemtoHelix(const MpdFemtoHelix &);

   /// Assignment operator (will use the one, provided by compiler)
   /// MpdFemtoHelix& operator=(const MpdFemtoHelix&);

   /// Destructor
   virtual ~MpdFemtoHelix();

   double dipAngle() const;
   double curvature() const; /// 1/R in xy-plane
   double phase() const;     /// aziumth in xy-plane measured from ring center
   double xcenter() const;   /// x-center of circle in xy-plane
   double ycenter() const;   /// y-center of circle in xy-plane
   int    h() const;         /// -sign(q*B);

   const TVector3 &origin() const; /// starting point

   /// Set helix parameters
   void setParameters(double c, double dip, double phase, const TVector3 &o, int h);

   /// coordinates of helix at point s
   double   x(double s) const;
   double   y(double s) const;
   double   z(double s) const;
   TVector3 at(double s) const;

   /// pointing vector of helix at point s
   double   cx(double s) const;
   double   cy(double s) const;
   double   cz(double s = 0) const;
   TVector3 cat(double s) const;

   /// returns period length of helix
   double period() const;

   /// path length at given r (cylindrical r)
   pair<double, double> pathLength(double r) const;

   /// path length at given r (cylindrical r, cylinder axis at x,y)
   pair<double, double> pathLength(double r, double x, double y);

   /// path length at distance of closest approach to a given point
   double pathLength(const TVector3 &p, bool scanPeriods = true) const;

   /// path length at intersection with plane
   double pathLength(const TVector3 &r, const TVector3 &n) const;

   /// path length at distance of closest approach in the xy-plane to a given point
   double pathLength(double x, double y) const;

   /// path lengths at dca between two helices
   pair<double, double> pathLength(const MpdFemtoHelix &h, double minStepSize = 10 * micrometer,
                                   double minRange = 10 * centimeter) const;

   /// minimal distance between point and helix
   double distance(const TVector3 &p, bool scanPeriods = true) const;

   /// checks for valid parametrization

   bool valid(double world = 1.e+5) const { return !bad(world); }
   int  bad(double world = 1.e+5) const;

   /// Move the origin along the helix to s which becomes then s=0
   virtual void moveOrigin(double s);

   static const double NoSolution;

protected:
   /// Set curvature of the helix
   void setCurvature(double); /// performs also various checks
   /// Set phase of the helix
   void setPhase(double);
   /// Set dip angle of the helix
   void setDipAngle(double);
   /// Value of S where distance in x-y plane is minimal
   double fudgePathLength(const TVector3 &) const;

protected:
   bool     mSingularity; /// true for straight line case (B=0)
   TVector3 mOrigin;      /// starting point of a helix
   double   mDipAngle;
   double   mCurvature;
   double   mPhase;
   int      mH; /// -sign(q*B);

   double mCosDipAngle;
   double mSinDipAngle;
   double mCosPhase;
   double mSinPhase;

   ClassDef(MpdFemtoHelix, 1);
};

//
//     Non-member functions
//
int           operator==(const MpdFemtoHelix &, const MpdFemtoHelix &);
int           operator!=(const MpdFemtoHelix &, const MpdFemtoHelix &);
std::ostream &operator<<(std::ostream &, const MpdFemtoHelix &);

//
//     Inline functions
//

inline int MpdFemtoHelix::h() const
{
   return mH;
}

inline double MpdFemtoHelix::dipAngle() const
{
   return mDipAngle;
}

inline double MpdFemtoHelix::curvature() const
{
   return mCurvature;
}

inline double MpdFemtoHelix::phase() const
{
   return mPhase;
}

inline double MpdFemtoHelix::x(double s) const
{
   if (mSingularity)
      return mOrigin.x() - s * mCosDipAngle * mSinPhase;
   else
      return mOrigin.x() + (cos(mPhase + s * mH * mCurvature * mCosDipAngle) - mCosPhase) / mCurvature;
}

inline double MpdFemtoHelix::y(double s) const
{
   if (mSingularity)
      return mOrigin.y() + s * mCosDipAngle * mCosPhase;
   else
      return mOrigin.y() + (sin(mPhase + s * mH * mCurvature * mCosDipAngle) - mSinPhase) / mCurvature;
}

inline double MpdFemtoHelix::z(double s) const
{
   return mOrigin.z() + s * mSinDipAngle;
}

inline double MpdFemtoHelix::cx(double s) const
{
   if (mSingularity)
      return -mCosDipAngle * mSinPhase;
   else
      return -sin(mPhase + s * mH * mCurvature * mCosDipAngle) * mH * mCosDipAngle;
}

inline double MpdFemtoHelix::cy(double s) const
{
   if (mSingularity)
      return mCosDipAngle * mCosPhase;
   else
      return cos(mPhase + s * mH * mCurvature * mCosDipAngle) * mH * mCosDipAngle;
}

inline double MpdFemtoHelix::cz(double /* s */) const
{
   return mSinDipAngle;
}

inline const TVector3 &MpdFemtoHelix::origin() const
{
   return mOrigin;
}

inline TVector3 MpdFemtoHelix::at(double s) const
{
   return TVector3(x(s), y(s), z(s));
}

inline TVector3 MpdFemtoHelix::cat(double s) const
{
   return TVector3(cx(s), cy(s), cz(s));
}

inline double MpdFemtoHelix::pathLength(double X, double Y) const
{
   return fudgePathLength(TVector3(X, Y, 0));
}

inline int MpdFemtoHelix::bad(double WorldSize) const
{

   int ierr;
   if (!::finite(mDipAngle)) {
      return 11;
   }
   if (!::finite(mCurvature)) {
      return 12;
   }

   // ierr = mOrigin.bad(WorldSize);

   /// The line above is commented and the StThreeVector::bad(double)
   /// is rewritten here
   for (int iIter = 0; iIter < 3; iIter++) {

      double tmpVal;
      /// Value StThreeVector.mX1[iter] ???
      switch (iIter) {
      case 0: tmpVal = mOrigin.X(); break;
      case 1: tmpVal = mOrigin.Y(); break;
      case 2: tmpVal = mOrigin.Z(); break;
      default: tmpVal = NAN;
      };

      if (!::finite(tmpVal)) {
         ierr = 10 + iIter;
      }
      if (::fabs(tmpVal) > WorldSize) {
         ierr = 20 + iIter;
      }
   } // for(int iIter=0; iIter<3; iIter+)

   if (ierr) {
      return (3 + ierr * 100);
   }
   if (::fabs(mDipAngle) > 1.58) {
      return 21;
   }

   double qwe = ::fabs(::fabs(mDipAngle) - M_PI / 2);
   if (qwe < 1. / WorldSize) {
      return 31;
   }

   if (::fabs(mCurvature) > WorldSize) {
      return 22;
   }
   if (mCurvature < 0) {
      return 32;
   }
   if (abs(mH) != 1) {
      return 24;
   }

   return 0;
}

#endif
