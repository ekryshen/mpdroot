//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_GEO_UTILS_H
#define __MPD_TOF_GEO_UTILS_H

#include <TVector3.h>
#include <TList.h>
#include <TGeoMatrix.h>

#include "IntervalTree.h"
//------------------------------------------------------------------------------------------------------------------------
class LRectangle // convex quadrangle
{
   bool IsInvalid;

   inline void getRPhi(const TVector3 &point, Double_t &Rmin, Double_t &Rmax, Double_t &Phimin, Double_t &Phimax)
   {
      const Double_t R = point.Perp(), Phi = point.Phi();
      if (R < Rmin) Rmin = R;
      if (R > Rmax) Rmax = R;
      if (Phi < Phimin) Phimin = Phi;
      if (Phi > Phimax) Phimax = Phi;
   }

public:
   enum Side_t { kRight = 0, kLeft = 1, kInvalid = -1 };

   Int_t    volumeUID;
   TVector3 A, B, C, D, center, perp; // [cm]

   LRectangle() : IsInvalid(true), volumeUID(kInvalid){};
   LRectangle(Int_t uid, const TVector3 &a, const TVector3 &b, const TVector3 &c, const TVector3 &d,
              bool check = false);

   TVector3 GetCenter() const { return (A + B + C + D) * 0.25; }
   bool     isInvalid() const { return IsInvalid; }
   Double_t GetIntersectionSquare(const LRectangle &v) const;
   void     GetRPhiRanges(Double_t &Rmin, Double_t &Rmax, Double_t &Phimin, Double_t &Phimax);
   Double_t DistanceFromPointToLineSegment(const TVector3 *pos, const TVector3 &P1, const TVector3 &P2) const;
   Double_t DistanceFromPointToLine(const TVector3 *pos, const TVector3 &P1, const TVector3 &P2) const;
   Double_t MinDistanceToEdge(const TVector3 *pos, Side_t &side) const;
   Double_t Angle(const LRectangle &r) const { return perp.Angle(r.perp); }; // [rad]
   void     Rotate(const TRotation &rot);
   void     Shift(const TVector3 &shift)
   {
      A += shift;
      B += shift;
      C += shift;
      D += shift;
   };

   void        Dump(const char *comment = nullptr, std::ostream &out = std::cout) const;
   static void Test(const char *comment = nullptr, std::ostream &out = std::cout);

   //-----------------------------------------------------
   inline double Square() const { return (B - A).Cross(D - A).Mag(); }
   //-----------------------------------------------------
   inline bool IsParallelPlane(const LRectangle &v) const
   {
      //		return (perp == v.perp);
      const static double epsilon = 1.e-6;
      return ((std::fabs(perp.X() - v.perp.X()) < epsilon) && (std::fabs(perp.Y() - v.perp.Y()) < epsilon) &&
              (std::fabs(perp.Z() - v.perp.Z()) < epsilon));
   }
   //-----------------------------------------------------
   inline std::pair<bool, double> IsSamePlane(const LRectangle &v) const
   {
      if (!IsParallelPlane(v)) return {false, 0.};

      const static double epsilon = 1.e-6;
      double              D1 = A * perp, D2 = v.A * v.perp, distance = fabs(D2 - D1);
      if (distance < epsilon) return {true, 0.};

      return {false, distance};
   }
   //-----------------------------------------------------
   inline void InitCenterPerp()
   {
      center = (A + B + C + D) * 0.25;
      perp   = (B - A).Cross(D - A).Unit();
   }
   //-----------------------------------------------------
   inline void CheckInValid()
   {
      IsInvalid = false;

      // Convex Polygon Definition: A polygon that has all interior angles less than 180°
      ;
      // Sum of Interior Angles, sum = 180*(n-2) degree, where n is the number of sides
      // A square has 4 sides, 	so interior angles sum = 360°
      ;
      // Rectangle check - all angles == 90 degree
      TVector3 ab = A - B, bc = B - C, cd = C - D, da = D - A;
      if (ab.Dot(bc) != 0. || bc.Dot(cd) != 0. || cd.Dot(da) != 0. || da.Dot(ab) != 0.) {
         std::cerr << "\n ---> ERROR: invalid Rectangle.";
         Dump("", std::cerr);
         IsInvalid = true;
      }
   }
   //-----------------------------------------------------
   inline TVector3 ProjectionToPlane(const TVector3 &point) const
   {
      Double_t d = perp * (A - point); // d = DotProduct ( N, A - point )
      return point + perp * d;
   }
   //-----------------------------------------------------
   inline LRectangle ProjectionToPlane(const LRectangle &v) const
   {
      LRectangle ret(v.volumeUID, ProjectionToPlane(v.A), ProjectionToPlane(v.B), ProjectionToPlane(v.C),
                     ProjectionToPlane(v.D));
      ret.InitCenterPerp();
      return ret;
   }
   //-----------------------------------------------------
   // valid ONLY if all(point and rectangle) at the same plane
   inline bool IsPointInside(const TVector3 &point) const
   {
      TVector3 v1  = A - B;
      TVector3 v2  = C - B;
      TVector3 v3  = point - B;
      double   mag = v1.Mag(), proj = v3 * v1 / mag;
      if (proj < 0. || proj > mag) return false;

      mag = v2.Mag(), proj = v3 * v2 / mag;
      if (proj < 0. || proj > mag) return false;
      return true;
   }
   //-----------------------------------------------------
};
//------------------------------------------------------------------------------------------------------------------------
class LStrip : public LRectangle {
public:
   Int_t          sectorID, detectorID, stripID;
   Int_t          neighboring[2]; // dim same as  Side_t enum
   TGeoCombiTrans fMatrix;

   LStrip();
   LStrip(Int_t uid, Int_t sector, Int_t detector, Int_t strip);

   void SetIDs(Int_t uid, Int_t sector, Int_t detector, Int_t strip)
   {
      volumeUID  = uid;
      sectorID   = sector;
      detectorID = detector;
      stripID    = strip;
   }

   inline bool IsSameDetector(const LStrip &strip) const
   {
      return (sectorID == strip.sectorID && detectorID == strip.detectorID);
   }
   inline bool operator==(const LStrip &rhs) const
   {
      return (sectorID == rhs.sectorID && detectorID == rhs.detectorID && stripID == rhs.stripID);
   }

   void     Dump(const char *comment = nullptr, std::ostream &out = std::cout) const;
   Double_t Distance(Side_t side, const LStrip &strip) const;
};
//------------------------------------------------------------------------------------------------------------------------
class MpdTofHitProducerQA;

typedef std::map<Int_t, LStrip>    TmStrips; // pair<suid, LStrip>
typedef Interval<LRectangle *>     Tinterval;
typedef IntervalTree<LRectangle *> TintervalTree;

class MpdTofGeoUtils {
   TmStrips      mStrips;                   //! mapping strips by suid
   TintervalTree mDetectorsZ, mStripsZ;     //! detectors Z[cm] location interval tree
   TintervalTree mDetectorsPhi, mStripsPhi; //! detectors Phi[rads] location interval tree

   static MpdTofGeoUtils *instance;

   void localToMaster(const TGeoMatrix *matrix, Double_t *local, TVector3 &position, std::pair<double, double> &Z,
                      std::pair<double, double> &Phi) const;

public:
   static MpdTofGeoUtils *Instance()
   {
      if (instance == nullptr) instance = new MpdTofGeoUtils;
      return instance;
   }

   // thresh, [cm] <--- thresh. distance between neighbor strips, (see QA histo)
   void FindNeighborStrips(Double_t thresh, MpdTofHitProducerQA *pQA, bool forced = false);
   void ParseTGeoManager(MpdTofHitProducerQA *, bool forced = false, const char *flnm = nullptr);
   // search by strip unique ID.
   const LStrip *FindStrip(Int_t suid) const;

   bool IsPointInsideStrips(const TVector3 &position, std::vector<Tinterval> &intersect, double Zerror = 0.,
                            double PhiError = 0.);
   bool IsPointInsideDetectors(const TVector3 &position, std::vector<Tinterval> &intersect, double Zerror = 0.,
                               double PhiError = 0.);

   const TintervalTree *GetDetZ() { return &mDetectorsZ; };
   const TintervalTree *GetDetPhi() { return &mDetectorsPhi; };

   const TintervalTree *GetStripZ() { return &mStripsZ; };
   const TintervalTree *GetStripPhi() { return &mStripsPhi; };
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_GEO_UTILS_H
