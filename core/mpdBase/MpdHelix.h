/***************************************************************************
 * Armen Kechechyan, September 2009
 ***************************************************************************
 * Description:
 * Parametrization of a helix
 **************************************************************************/

#ifndef MPD_HELIX_H
#define MPD_HELIX_H

#include "TVector3.h"
#include <math.h>
#include <utility>
#include <algorithm>

#if !defined(NO_NAMESPACES)
using std::pair;
using std::swap;
using std::max;
#endif

class MpdHelix {
public:
    /// curvature, dip angle, phase, origin, h
    MpdHelix(double c, double dip, double phase,
	    const TVector3 o, int h=-1);
    MpdHelix(TVector3 mom, TVector3 o, Double_t charge, Double_t Bz=0.5);
    
    virtual ~MpdHelix();
    // MpdHelix(const MpdHelix&);			// use default
    // MpdHelix& operator=(const MpdHelix&);	// use default

    double       dipAngle()   const;           
    double       curvature()  const;	/// 1/R in xy-plane
    double       phase()      const;	/// aziumth in xy-plane measured from ring center
    double       xcenter()    const;	/// x-center of circle in xy-plane
    double       ycenter()    const;	/// y-center of circle in xy-plane
    int          h()          const;	/// -sign(q*B);
    
    const TVector3 origin() const;	/// starting point

    void setParameters(double c, double dip, double phase, const TVector3 o, int h);

    /// coordinates of helix at point s
    double       x(double s)  const;
    double       y(double s)  const;
    double       z(double s)  const;

    TVector3  at(double s) const;

    /// pointing vector of helix at point s
    double       cx(double s)  const;
    double       cy(double s)  const;
    double       cz(double s)  const;
    
    TVector3  cat(double s) const;

    /// returns period length of helix
    double       period()       const;
    
    /// path length at given r (cylindrical r)
    pair<double, double> pathLength(double r)   const;
    
    /// path length at given r (cylindrical r, cylinder axis at x,y)
    pair<double, double> pathLength(double r, double x, double y);
    
    /// path length at distance of closest approach to a given point
    double       pathLength(const TVector3 p, bool scanPeriods = true) const;
    
    /// path length at intersection with plane
    double       pathLength(const TVector3 r, const TVector3 n) const;

    /// path length at distance of closest approach in the xy-plane to a given point
    double       pathLength(double x, double y) const;

    /// path lengths at dca between two helices 
    pair<double, double> pathLengths(const MpdHelix&) const;
    
    /// minimal distance between point and helix
    double       distance(const TVector3 p, bool scanPeriods = true) const;    
    
    /// checks for valid parametrization
    bool         valid(double world = 1.e+5) const {return !bad(world);}
    int            bad(double world = 1.e+5) const;
    
    /// move the origin along the helix to s which becomes then s=0
    virtual void moveOrigin(double s);
    
    static const double NoSolution;
    
protected:
    MpdHelix();
    
    void setCurvature(double);	/// performs also various checks   
    void setPhase(double);	        
    void setDipAngle(double);
    
    /// value of S where distance in x-y plane is minimal
    double fudgePathLength(const TVector3) const;
    
protected:
    bool                   mSingularity;	// true for straight line case (B=0)
    TVector3  mOrigin;
    double                 mDipAngle;
    double                 mCurvature;
    double                 mPhase;
    int                    mH;			// -sign(q*B);

    double                 mCosDipAngle;
    double                 mSinDipAngle;
    double                 mCosPhase;
    double                 mSinPhase;
#ifdef __ROOT__
  ClassDef(MpdHelix,1)
#endif
};

//
//     Non-member functions
//
int operator== (const MpdHelix&, const MpdHelix&);
int operator!= (const MpdHelix&, const MpdHelix&);
//ostream& operator<<(ostream&, const MpdHelix&);

//
//     Inline functions
//
inline int MpdHelix::h() const {return mH;}

inline double MpdHelix::dipAngle() const {return mDipAngle;}

inline double MpdHelix::curvature() const {return mCurvature;}

inline double MpdHelix::phase() const {return mPhase;}

inline double MpdHelix::x(double s) const
{
    if (mSingularity)
	return mOrigin.X() - s*mCosDipAngle*mSinPhase;
    else
	return mOrigin.X() + (cos(mPhase + s*mH*mCurvature*mCosDipAngle)-mCosPhase)/mCurvature;
}
 
inline double MpdHelix::y(double s) const
{
    if (mSingularity)
	return mOrigin.Y() + s*mCosDipAngle*mCosPhase;
    else
	return mOrigin.Y() + (sin(mPhase + s*mH*mCurvature*mCosDipAngle)-mSinPhase)/mCurvature;
}

inline double MpdHelix::z(double s) const
{
    return mOrigin.Z() + s*mSinDipAngle;
}

inline double MpdHelix::cx(double s)  const
{
    if (mSingularity)
	return -mCosDipAngle*mSinPhase;
    else
	return -sin(mPhase + s*mH*mCurvature*mCosDipAngle)*mH*mCosDipAngle;
}

inline double MpdHelix::cy(double s)  const
{
    if (mSingularity)
	return mCosDipAngle*mCosPhase;
    else
	return cos(mPhase + s*mH*mCurvature*mCosDipAngle)*mH*mCosDipAngle;
}

inline double MpdHelix::cz(double s)  const
{
    return mSinDipAngle;
}    

inline const TVector3 MpdHelix::origin() const {return mOrigin;}

inline TVector3 MpdHelix::at(double s) const
{
    return TVector3(x(s), y(s), z(s));
}

inline TVector3 MpdHelix::cat(double s) const
{
    return TVector3(cx(s), cy(s), cz(s));
}

inline double MpdHelix::pathLength(double xx, double yy) const
{
    return fudgePathLength(TVector3(xx, yy, 0));
}

inline int MpdHelix::bad(double WorldSize) const
{

    int ierr;
    if (!::finite(mDipAngle    )) 	return   11;
    if (!::finite(mCurvature   )) 	return   12;

//    ierr = mOrigin.bad(WorldSize);
   ierr = 1;
    if (ierr)                           return    3+ierr*100;

    if (::fabs(mDipAngle)  >1.58)	return   21;
    double qwe = ::fabs(::fabs(mDipAngle)-M_PI/2);
    if (qwe < 1./WorldSize      ) 	return   31; 

    if (::fabs(mCurvature) > WorldSize)	return   22;
    if (mCurvature < 0          )	return   32;

    if (abs(mH) != 1            )       return   24; 

    return 0;
}

#endif
