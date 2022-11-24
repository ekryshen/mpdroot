//-------------------------------------------------------------------------------------------------
// Description:
//      BaseTpcSectorGeo class is the base class for TPC Geometry,
//      with its configuration in ideal, non-misaligned state.
//      Descriptive drawing:
//      https://git.jinr.ru/nica/docs/-/blob/main/docs/mpdroot/coding/geometry/TPC_coordinates.svg
//
//      Any custom TPC geometry must inherit from this class.
//
//      Note: If you are overriding any method, add in front of it the 'virtual' keyword.
//            All class members are constants,
//            i.e. 'const' qualifier at the end of all methods is omitted.
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Alexander Bychkov, Slavomir Hnatic
//      JINR, October, 2022
//-------------------------------------------------------------------------------------------------

#ifndef BASETPCSECTORGEO_HH
#define BASETPCSECTORGEO_HH

// ROOT Class Headers ---------------
#include <TObject.h>
#include <TMath.h>
#include <TVector2.h>
#include <TVector3.h>

class BaseTpcSectorGeo : public TObject {
public:
   // Constructors/Destructors ---------
   BaseTpcSectorGeo();
   virtual ~BaseTpcSectorGeo();

   /* transformations */
   // sectors: angle and number
   double SectorAxisAngleRad(int iSector) { return iSector * SECTOR_PHI_RAD; } // angle of sector axis in radians
   int    SectorNumberFromGlobal(const TVector3 &globalXYZ);

   // padplane: padrows to local coordinates and vice versa
   TVector2                  PadRow2Local(double pad, double row);
   TVector2                  PadRowCenter2Local(int padNumber, int rowNumber);
   std::pair<double, double> Local2PadRow(const TVector3 &localXYZ);

   // time and z axis conversions
   int    Time2TimeBin(double time) { return TMath::FloorNint(time / TIMEBIN_LENGTH); }
   double TimeBin2Time(int timeBin) { return timeBin * TIMEBIN_LENGTH; }
   double TimeBin2Z(int timeBin, double velocity) { return velocity * TimeBin2Time(timeBin); }
   int    Z2TimeBin(double z, double velocity) { return Time2TimeBin(z / velocity); }

   // global to local coordinates and vice versa
   TVector3                 Global2Local(const TVector3 &globalXYZ, int &iSector);
   std::pair<TVector3, int> Global2Local(const TVector3 &globalXYZ);
   TVector3                 Local2Global(const TVector3 &localXYZ, int iSector);
   TVector3                 Local2Global(const std::pair<TVector3, int> &localPosition);

   /* constants */
   enum PadArea : int { inner, outer };
   enum PadAreaBoundary : int { lowerEdge, midBoundary, upperEdge };

   const int    SECTOR_COUNT      = 24;                     // total number of sectors
   const int    SECTOR_COUNT_HALF = 12;                     // number of sectors on one side
   const double SECTOR_PHI_RAD    = 30 * TMath::DegToRad(); // sector angle
   const double SECTOR_PHI0_RAD   = -SECTOR_PHI_RAD / 2.;   // first sector angle

   const double Z_MIN        = 0.01;                 // half of membrane thickness
   const double DRIFT_LENGTH = 163.879;              // drift length
   const double Z_MAX        = Z_MIN + DRIFT_LENGTH; // distance to pads in cm

   const int    TIMEBIN_COUNT  = 310;  // number of timebins
   const double TIMEBIN_LENGTH = 100.; // timebin length, ns
   const double TIMEDRIFT_MAX  = TIMEBIN_COUNT * TIMEBIN_LENGTH;

   const std::vector<int>    ROW_COUNT{27, 26};    // number of padrows for inner & outer regions
   const std::vector<double> PAD_HEIGHT{1.2, 1.8}; // height of pads in inner & outer regions
   const std::vector<double> PAD_WIDTH{0.5, 0.5};  // width of pads in inner & outer regions
   const std::vector<int>    PAD_COUNT{20, 21, 21, 22, 23, 23, 24, 24, 25, 26, 26, 27, 28, 28, 29, 30, 30, 31,
                                    32, 32, 33, 33, 34, 35, 35, 36, 37, 38, 39, 40, 41, 41, 42, 43, 44, 45,
                                    46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62};
   // half number of pads in rows

   const double YPADPLANE_OFFSET         = 40.3; // Y of padplane's low edge in global coordinates
   const double YPADPLANE2PADAREA_OFFSET = 0.4;  // distance between padplane's lower edge and padarea
   const double YPADAREA_LOWEREDGE       = YPADPLANE_OFFSET + YPADPLANE2PADAREA_OFFSET;
   // Y length of padareas: inner, outer
   const std::vector<double> YPADAREA_LENGTH{ROW_COUNT[inner] * PAD_HEIGHT[inner],
                                             ROW_COUNT[outer] * PAD_HEIGHT[outer]};
   // Y coordinates of padareas in local & global coordinates: inner edge, boundary inner/outer, outer edge
   const std::vector<double> YPADAREA_LOCAL{0., YPADAREA_LENGTH[inner],
                                            YPADAREA_LENGTH[inner] + YPADAREA_LENGTH[outer]};
   const std::vector<double> YPADAREA_GLOBAL{YPADAREA_LOWEREDGE, YPADAREA_LOWEREDGE + YPADAREA_LOCAL[midBoundary],
                                             YPADAREA_LOWEREDGE + YPADAREA_LOCAL[upperEdge]};

protected:
private:
   ClassDef(BaseTpcSectorGeo, 1);
};

#endif
