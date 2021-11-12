/***************************************************************************
 * Armen Kechechyan, September 2009
 ***************************************************************************
 * Description:
 * Parametrization of a helix
 **************************************************************************/

#if !defined(NO_NUMERIC_LIMITS)
#    include <limits>
#    if !defined(NO_NAMESPACES)
         using std::numeric_limits;
#    endif
#endif

#include "MpdBase.h" // for basic SI units

#define FOR_HELIX

#include "MpdHelix.h"
#include <TMath.h>

#ifdef __ROOT__
ClassImpT(MpdHelix,double);
#endif

const double MpdHelix::NoSolution = 3.e+33;

MpdHelix::MpdHelix():
		mSingularity(false),
		mOrigin(),
		mDipAngle(0),
		mCurvature(0),
		mPhase(0),
		mH(0),
		mCosDipAngle(0),
		mSinDipAngle(0),
		mCosPhase(0),
		mSinPhase(0)
	{ /*noop*/ }

MpdHelix::MpdHelix(double c, double d, double Phase,
		 const TVector3 o, int hh)
{
    setParameters(c, d, Phase, o, hh);
}

MpdHelix::MpdHelix(TVector3 mom, TVector3 o, Double_t charge,
		Double_t Bz) {
	Double_t temp_h =  -TMath::Sign(1.,charge*Bz);
	Double_t dip_angle = TMath::ATan2(mom.Pz(),mom.Pt());
	Double_t curv = TMath::Abs(charge*Bz*0.299792458/mom.Pt())*0.01;
	Double_t Phase =  mom.Phi()-temp_h*TMath::PiOver2();
	setParameters(curv,dip_angle,Phase,o,temp_h);
}
MpdHelix::~MpdHelix() { /* noop */ };

void MpdHelix::setParameters(double c, double dip, double Phase,
			    			const TVector3 o, int hh)
{
    //
    //  The order in which the parameters are set is important
    //  since setCurvature might have to adjust the others.
    //
    mH = (hh>=0) ? 1 : -1;    // Default is: positive particle
                             //             positive field
    mOrigin   = o;
    setDipAngle(dip);
    setPhase(Phase);

    //
    // Check for singularity and correct for negative curvature.           
    // May change mH and mPhase. Must therefore be set last.
    //
    setCurvature(c);

    //
    // For the case B=0, h is ill defined. In the following we
    // always assume h = +1. Since phase = psi - h * pi/2
    // we have to correct the phase in case h = -1.
    // This assumes that the user uses the same h for phase
    // as the one he passed to the constructor.
    //
    if (mSingularity && mH == -1) {
	mH = +1;
	setPhase(mPhase-M_PI);
    }
}

void MpdHelix::setCurvature(double val)
{
    if (val < 0) {
	mCurvature = -val;
	mH = -mH;
	setPhase(mPhase+M_PI);
    }
    else
	mCurvature = val;

#ifndef NO_NUMERIC_LIMITS
    if (fabs(mCurvature) <= numeric_limits<double>::epsilon())
#else
    if (fabs(mCurvature) <= static_cast<double>(0))
#endif    
	mSingularity = true;			// straight line
    else
	mSingularity = false;            	// curved
}

void MpdHelix::setPhase(double val)
{
    mPhase       = val;
    mCosPhase    = cos(mPhase);
    mSinPhase    = sin(mPhase);
    if (fabs(mPhase) > M_PI)
	mPhase = atan2(mSinPhase, mCosPhase);  // force range [-pi,pi]
}

void MpdHelix::setDipAngle(double val)
{
    mDipAngle    = val;
    mCosDipAngle = cos(mDipAngle);
    mSinDipAngle = sin(mDipAngle);
}

double MpdHelix::xcenter() const
{
    if (mSingularity)
	return 0;
    else
	return mOrigin.X()-mCosPhase/mCurvature;
}

double MpdHelix::ycenter() const
{
    if (mSingularity)
	return 0;
    else
	return mOrigin.Y()-mSinPhase/mCurvature;
}

double MpdHelix::fudgePathLength(const TVector3 p) const
{
    double s;
    double dx = p.X()-mOrigin.X();
    double dy = p.Y()-mOrigin.Y();
    
    if (mSingularity) {
	s = (dy*mCosPhase - dx*mSinPhase)/mCosDipAngle;
    }
    else {
	s = atan2(dy*mCosPhase - dx*mSinPhase,
		  1/mCurvature + dx*mCosPhase+dy*mSinPhase)/
	    (mH*mCurvature*mCosDipAngle);
    }
    return s;
}

double MpdHelix::distance(const TVector3 p, bool scanPeriods) const
{
    return (this->at(pathLength(p,scanPeriods))-p).Mag();
}

double MpdHelix::pathLength(const TVector3 p, bool scanPeriods) const 
{
    //
    //  Returns the path length at the distance of closest 
    //  approach between the helix and point p. 
    //  For the case of B=0 (straight line) the path length
    //  can be calculated analytically. For B>0 there is
    //  unfortunately no easy solution to the problem.
    //  Here we use the Newton method to find the root of the
    //  referring equation. The 'fudgePathLength' serves
    //  as a starting value.
    //
    double s;
    double dx = p.X()-mOrigin.X();
    double dy = p.Y()-mOrigin.Y();
    double dz = p.Z()-mOrigin.Z();

    if (mSingularity) {
	s = mCosDipAngle*(mCosPhase*dy-mSinPhase*dx) +
	    mSinDipAngle*dz;
    }
    else { //
// #ifndef NO_NAMESPACES
 	{
// 	    using namespace units;
// #endif
//	    const double MaxPrecisionNeeded = 0.000001;  // micrometer;
	    const double MaxPrecisionNeeded = SI_units::micrometer;
	    const int    MaxIterations      = 100;

	    //
	    // The math is taken from Maple with C(expr,optimized) and
	    // some hand-editing. It is not very nice but efficient.
	    //
	    double t34 = mCurvature*mCosDipAngle*mCosDipAngle;
	    double t41 = mSinDipAngle*mSinDipAngle;
	    double t6, t7, t11, t12, t19;

	    //
	    // Get a first guess by using the dca in 2D. Since
	    // in some extreme cases we might be off by n periods
	    // we add (subtract) periods in case we get any closer.
	    // 
	    s = fudgePathLength(p);

	    if (scanPeriods) {
	        double ds = period();
	        int    j, jmin = 0;
	        double d, dmin = (at(s) - p).Mag();
	        for(j=1; j<MaxIterations; j++) {
		  if ((d = (at(s+j*ds) - p).Mag()) < dmin) {
		      dmin = d;
		      jmin = j;
		  }
		  else
		      break;
	        }
	        for(j=-1; -j<MaxIterations; j--) {
		  if ((d = (at(s+j*ds) - p).Mag()) < dmin) {
		      dmin = d;
		      jmin = j;
		  }
		  else
		      break;
	        }
	        if (jmin) s += jmin*ds;
	    }
	    
	    //
	    // Newtons method:
	    // Stops after MaxIterations iterations or if the required
	    // precision is obtained. Whatever comes first.
	    //
	    double sOld = s;
	    for (int i=0; i<MaxIterations; i++) {
		t6  = mPhase+s*mH*mCurvature*mCosDipAngle;
		t7  = cos(t6);
		t11 = dx-(1/mCurvature)*(t7-mCosPhase);
		t12 = sin(t6);
		t19 = dy-(1/mCurvature)*(t12-mSinPhase);
		s  -= (t11*t12*mH*mCosDipAngle-t19*t7*mH*mCosDipAngle -
		       (dz-s*mSinDipAngle)*mSinDipAngle)/
		    (t12*t12*mCosDipAngle*mCosDipAngle+t11*t7*t34 +
		     t7*t7*mCosDipAngle*mCosDipAngle +
		     t19*t12*t34+t41);
		if (fabs(sOld-s) < MaxPrecisionNeeded) break;
		sOld = s;
	    }
//#ifndef NO_NAMESPACES
	}
//#endif
    }
    return s;
}

double MpdHelix::period() const
{
    if (mSingularity)
#ifndef NO_NUMERIC_LIMITS
            return numeric_limits<double>::max();
#else
            return DBL_MAX;
#endif    
    else	
	return fabs(2*M_PI/(mH*mCurvature*mCosDipAngle)); 
}

pair<double, double> MpdHelix::pathLength(double r) const
{
    pair<double,double> value;
    pair<double,double> VALUE(999999999.,999999999.);
    //
    // The math is taken from Maple with C(expr,optimized) and
    // some hand-editing. It is not very nice but efficient.
    // 'first' is the smallest of the two solutions (may be negative)
    // 'second' is the other.
    //
    if (mSingularity) {
	double t1 = mCosDipAngle*(mOrigin.X()*mSinPhase-mOrigin.Y()*mCosPhase);
	double t12 = mOrigin.Y()*mOrigin.Y();
	double t13 = mCosPhase*mCosPhase;
	double t15 = r*r;
	double t16 = mOrigin.X()*mOrigin.X();
	double t20 = -mCosDipAngle*mCosDipAngle*(2.0*mOrigin.X()*mSinPhase*mOrigin.Y()*mCosPhase +
				 t12-t12*t13-t15+t13*t16);
	if (t20<0.) return VALUE;
	t20 = ::sqrt(t20);
	value.first  = (t1-t20)/(mCosDipAngle*mCosDipAngle);
	value.second = (t1+t20)/(mCosDipAngle*mCosDipAngle);
    }
    else {
	double t1 = mOrigin.Y()*mCurvature;
	double t2 = mSinPhase;
	double t3 = mCurvature*mCurvature;
	double t4 = mOrigin.Y()*t2;
	double t5 = mCosPhase;
	double t6 = mOrigin.X()*t5;
	double t8 = mOrigin.X()*mOrigin.X();
	double t11 = mOrigin.Y()*mOrigin.Y();
	double t14 = r*r;
	double t15 = t14*mCurvature;
	double t17 = t8*t8;
	double t19 = t11*t11;
	double t21 = t11*t3;
	double t23 = t5*t5;
	double t32 = t14*t14;
	double t35 = t14*t3;
	double t38 = 8.0*t4*t6 - 4.0*t1*t2*t8 - 4.0*t11*mCurvature*t6 +
	             4.0*t15*t6 + t17*t3 + t19*t3 + 2.0*t21*t8 + 4.0*t8*t23 -
	             4.0*t8*mOrigin.X()*mCurvature*t5 - 4.0*t11*t23 -
	             4.0*t11*mOrigin.Y()*mCurvature*t2 + 4.0*t11 - 4.0*t14 +
	             t32*t3 + 4.0*t15*t4 - 2.0*t35*t11 - 2.0*t35*t8;
	double t40 = (-t3*t38);
	if (t40<0.) return VALUE;
	t40 = ::sqrt(t40);
	
	double t43 = mOrigin.X()*mCurvature;
	double t45 = 2.0*t5 - t35 + t21 + 2.0 - 2.0*t1*t2 -2.0*t43 - 2.0*t43*t5 + t8*t3;
	double t46 = mH*mCosDipAngle*mCurvature;
	
	value.first = (-mPhase + 2.0*atan((-2.0*t1 + 2.0*t2 + t40)/t45))/t46;
	value.second = -(mPhase + 2.0*atan((2.0*t1 - 2.0*t2 + t40)/t45))/t46;

	//
	//   Solution can be off by +/- one period, select smallest
	//
	double p = period();
    if (!std::isnan(value.first)) {
	    if (fabs(value.first-p) < fabs(value.first)) value.first = value.first-p;
	    else if (fabs(value.first+p) < fabs(value.first)) value.first = value.first+p;
	}
    if (!std::isnan(value.second)) {
	    if (fabs(value.second-p) < fabs(value.second)) value.second = value.second-p;
	    else if (fabs(value.second+p) < fabs(value.second)) value.second = value.second+p;
	}
    }
    if (value.first > value.second)
	swap(value.first,value.second);
    return(value);
}

pair<double, double> MpdHelix::pathLength(double r, double xx, double yy)
{
    double x0 = mOrigin.X();
    double y0 = mOrigin.Y();
    mOrigin.SetX(x0-xx);
    mOrigin.SetY(y0-yy);
    pair<double, double> result = this->pathLength(r);
    mOrigin.SetX(x0);
    mOrigin.SetY(y0);
    return result;  
}

double MpdHelix::pathLength(const TVector3 r,
		           const TVector3 n) const
{
    //
    // Vector 'r' defines the position of the center and
    // vector 'n' the normal vector of the plane.
    // For a straight line there is a simple analytical
    // solution. For curvatures > 0 the root is determined
    // by Newton method. In case no valid s can be found
    // the max. largest value for s is returned.
    //
    double s;

    if (mSingularity) {
	double t = n.Z()*mSinDipAngle +
	           n.Y()*mCosDipAngle*mCosPhase -
	           n.X()*mCosDipAngle*mSinPhase;
	if (t == 0)
	    s = NoSolution;
	else
	    s = ((r - mOrigin)*n)/t;
    }
    else {
//        const double MaxPrecisionNeeded = 0.000001;  // micrometer;
        const double MaxPrecisionNeeded = SI_units::micrometer;
        const int    MaxIterations      = 20;
        	
	double A = mCurvature*((mOrigin - r)*n) -
	           n.X()*mCosPhase - 
	           n.Y()*mSinPhase;
	double t = mH*mCurvature*mCosDipAngle;
	double u = n.Z()*mCurvature*mSinDipAngle;
	
	double a, f, fp;
	double sOld = s = 0;  
	double shiftOld = 0;
	double shift;
//		(cos(angMax)-1)/angMax = 0.1
        const double angMax = 0.21;
        double deltas = fabs(angMax/(mCurvature*mCosDipAngle));
//              dampingFactor = exp(-0.5);
	double dampingFactor = 0.60653;
	int i;

	for (i=0; i<MaxIterations; i++) {
	    a  = t*s+mPhase;
            double sina = sin(a);
            double cosa = cos(a);
	    f  = A +
		 n.X()*cosa +
		 n.Y()*sina +
		 u*s;
	    fp = -n.X()*sina*t +
		  n.Y()*cosa*t +
		  u;
            if ( fabs(fp)*deltas <= fabs(f) ) { //too big step
               int sgn = 1;
               if (fp<0.) sgn = -sgn;
               if (f <0.) sgn = -sgn;
	       shift = sgn*deltas;
               if (shift<0) shift*=0.9;  // don't get stuck shifting +/-deltas
            } else {
               shift = f/fp;
            }
	    s -= shift;
	    shiftOld = shift;
	    if (fabs(sOld-s) < MaxPrecisionNeeded) break;
	    sOld = s;
	}
        if (i == MaxIterations) return NoSolution;
    }
    return s;
}

pair<double, double>
MpdHelix::pathLengths(const MpdHelix& hh) const
{

    //
    //	Cannot handle case where one is a helix
    //  and the other one is a straight line.
    //
    if (mSingularity != hh.mSingularity) 
	return pair<double, double>(NoSolution, NoSolution);

    double s1, s2;

    if (mSingularity) {
	//
	//  Analytic solution
	//
	TVector3 dv = hh.mOrigin - mOrigin;
	TVector3 a(-mCosDipAngle*mSinPhase,
				mCosDipAngle*mCosPhase,
				mSinDipAngle);
	TVector3 b(-hh.mCosDipAngle*hh.mSinPhase,
				hh.mCosDipAngle*hh.mCosPhase,
				hh.mSinDipAngle);	
	double ab = a*b;
	double g  = dv*a;
	double k  = dv*b;
	s2 = (k-ab*g)/(ab*ab-1.);
	s1 = g+s2*ab;
	return pair<double, double>(s1, s2);
    }
    else {	
	//
	//  First step: get dca in the xy-plane as start value
	//
	double dx = hh.xcenter() - xcenter();
	double dy = hh.ycenter() - ycenter();
	double dd = ::sqrt(dx*dx + dy*dy);
	double r1 = 1/curvature();
	double r2 = 1/hh.curvature();
	
	double cosAlpha = (r1*r1 + dd*dd - r2*r2)/(2*r1*dd);
	
	double s;
	double xx, yy;
	if (fabs(cosAlpha) < 1) {           // two solutions
	    double sinAlpha = sin(acos(cosAlpha));
	    xx = xcenter() + r1*(cosAlpha*dx - sinAlpha*dy)/dd;
	    yy = ycenter() + r1*(sinAlpha*dx + cosAlpha*dy)/dd;
	    s = pathLength(xx, yy);
	    xx = xcenter() + r1*(cosAlpha*dx + sinAlpha*dy)/dd;
	    yy = ycenter() + r1*(cosAlpha*dy - sinAlpha*dx)/dd;
	    double a = pathLength(xx, yy);
	    if (hh.distance(at(a)) < hh.distance(at(s))) s = a;
	}
	else {                              // no intersection (or exactly one)
	    int rsign = ((r2-r1) > dd ? -1 : 1); // set -1 when *this* helix is
                                                   // completely contained in the other              
	    xx = xcenter() + rsign*r1*dx/dd;
	    yy = ycenter() + rsign*r1*dy/dd;
	    s = pathLength(xx, yy);
	}
	
	//
	//   Second step: scan in decreasing intervals around seed 's'
	// 
//	const double MinStepSize = 0.001;  // 10*micrometer;
//	const double MinRange    = 10.;  // 10*centimeter;    
	const double MinStepSize = 10*SI_units::micrometer;
	const double MinRange    = 10*SI_units::centimeter;    
	double dmin              = hh.distance(at(s));
	double range             = max(2*dmin, MinRange);
	double ds                = range/10;
	double slast=-999999, ss, d;
	s1 = s - range/2.;
	s2 = s + range/2.;
	double d1 = 0, d2 = 0, d3 = 0, s11 = 0, s12 = 0, s13 = 0;
	while (ds > MinStepSize) {
		int i = 0;
	    for (ss=s1; ss<s2+ds; ss+=ds) {
			d = hh.distance(at(ss));
			i++;
			if (i > 10){
		        if (i == 11) {d1 = d; s11 = s;}
                if (i == 12) {d2 = d; s12 = s;}
                if (i == 13) {d3 = d; s13 = s;}
                if (i == 14 && d1 == d2 && d2 == d3 && s11 == s12 && s12 == s13) {
                    printf("REPEATING VALUES IN MpdHelix\n");
                    return pair<double, double>(s, hh.pathLength(at(s)));
                }
            }
			if (d < dmin) {
		    	dmin = d;
		    	s = ss;
			}	
			slast = ss;
	    }
	    //
	    //  In the rare cases where the minimum is at the
	    //  the border of the current range we shift the range
	    //  and start all over, i.e we do not decrease 'ds'.
	    //  Else we decrease the search intervall around the
	    //  current minimum and redo the scan in smaller steps.
	    //
	    if (s == s1) {
		d = 0.8*(s2-s1);
		s1 -= d;
		s2 -= d;
	    }
	    else if (s == slast) {
		d = 0.8*(s2-s1);
		s1 += d;
		s2 += d;
	    }
	    else {           
		s1 = s-ds;
		s2 = s+ds;
		ds /= 10;
	    }
	}
	return pair<double, double>(s, hh.pathLength(at(s)));
    }
}


void MpdHelix::moveOrigin(double s)
{
    if (mSingularity)
	mOrigin	= at(s);
    else {
	TVector3 newOrigin = at(s);
	double newPhase = atan2(newOrigin.Y() - ycenter(),
				newOrigin.X() - xcenter());
	mOrigin = newOrigin;
	setPhase(newPhase);	        
    }
}

int operator== (const MpdHelix& a, const MpdHelix& b)
{
    //
    // Checks for numerical identity only !
    //
    return (a.origin()    == b.origin()    &&
	    a.dipAngle()  == b.dipAngle()  &&
	    a.curvature() == b.curvature() &&
	    a.phase()     == b.phase()     &&
	    a.h()         == b.h());
}

int operator!= (const MpdHelix& a, const MpdHelix& b) {return !(a == b);}

// ostream& operator<<(ostream& os, const MpdHelix& h)
// {
//     return os << '('
// 	      << "curvature = "  << h.curvature() << ", " 
// 	      << "dip angle = "  << h.dipAngle()  << ", "
// 	      << "phase = "      << h.phase()     << ", "  
// 	      << "h = "          << h.h()         << ", "    
// 	      << "origin = "     << h.origin()    << ')';
// }
// 
