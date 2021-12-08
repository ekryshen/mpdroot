// -------------------------------------------------------------------------
// -----          TShieldGeometry header file                          -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

/**  TShieldGeometry.h
 *@author Dmitry Sosnov (02.2015) <dmitry.e.sosnov@gmail.com>
 *
 * The TShieldGeometry is set of functions for converting geometry
 * from Root's geometry (TGeoManager) to SHIELD's geometry
 * for using at TShield.
 **/

#ifndef TSHIELD_GEOMETRY_H
#define TSHIELD_GEOMETRY_H

#include "shield.h"
#include <cstdio>
#include <cmath>
#include <set>
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TGeoNode.h"

#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoSphere.h"
#include "TGeoPcon.h"
#include "TGeoPgon.h"
#include "TGeoEltu.h"
#include "TGeoArb8.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"

// Not implemented
#include "TGeoPara.h"
#include "TGeoHype.h"
#include "TGeoTorus.h"
#include "TGeoParaboloid.h"
#include "TGeoXtru.h"
#include "TGeoHalfSpace.h"

//-----------------------------------------------------------------------------------------------------------------------------
// Types:
// Output type (tuple) with data for adding to shield
struct Geant4ShieldElement {
   std::vector<int>      zoneVector;
   std::vector<SGeoBody> bodyVector;
   MediumData            medium;
};
typedef std::vector<Geant4ShieldElement> Geant4ShieldData;
namespace tgeanttoshield {
inline double getDefaultScale()
{
   double cm = 1;
   return 1 / cm;
}

// Element for using in boolean calculations at geometry.
typedef std::pair<int, SGeoBody> zoneElement;
// vector zoneList is:
// zone, written as
// "(A1 A2 A3 A4) OR (B1 B2 B3 B4) OR (C1 C2 C3 C4)"
// where "(A1 A2)" is "A1 AND A2"
// at shield input files
// at this structure is:
// +---+-------------+
// | 1 | A1 A2 A3 A4 |
// +---+-------------+
// | 2 | B1 B2 B3 B4 |
// +---+-------------+
// | 3 | C1 C2 C3 C4 |
// +---+-------------+
typedef std::vector<std::vector<zoneElement>> zoneList;
// Type for vector of zones for Shield
typedef std::vector<std::pair<zoneList, MediumData>> zoneData;
} // namespace tgeanttoshield

// Functions:
//-----------------------------------------------------------------------------------------------------------------------------
Geant4ShieldData getDataAtPoint(TGeoManager *geoMan, Double_t x, Double_t y, Double_t z, Double_t vx, Double_t vy,
                                Double_t vz);
void             printGeometry(TGeoManager *geoMan); // at TShieldGeometryPrint.cxx
namespace tgeanttoshield {
zoneData         getCurrentZone(TGeoNavigator *nav);
Geant4ShieldData zoneDataToShield(zoneData data);
// Add volume with "outer vacuum" media, larger then top volume
std::pair<zoneList, MediumData> addOuterVacuum(TGeoNode *world);
zoneData                        getZonesAtPoint(TGeoNavigator *nav, Double_t x, Double_t y, Double_t z);

//-----------------------------------------------------------------------------------------------------------------------------
// TShieldGeometryBool.cxx
// Boolean algebra function on zoneElement and zoneData
zoneElement notElement(zoneElement el);
zoneList    notZone(zoneList list);
zoneList    orZone(zoneList list1, zoneList list2);
zoneList    andZone(zoneList list1, zoneList list2);

//-----------------------------------------------------------------------------------------------------------------------------
// TShieldGeometryPrint.cxx
// Print used values
void printElement(Element element);
void printMediumData(MediumData medium);
void printSGeoBody(SGeoBody body);
void printSGeoZone(SGeoZone zone);
void printGeant4ShieldData(Geant4ShieldData data);
void printGeant4ShieldElement(Geant4ShieldElement element);
void printGeometryBranch(TGeoNavigator *nav);

//-----------------------------------------------------------------------------------------------------------------------------
// TShieldGeometryOperators.cxx
bool            doubleNE(const double &l, const double &r);
bool            doubleEQ(const double &l, const double &r);
bool            doubleLT(const double &l, const double &r);
bool            doubleLE(const double &l, const double &r);
bool            doubleGT(const double &l, const double &r);
bool            doubleGE(const double &l, const double &r);
TGeoHMatrix     operator*(TGeoHMatrix &matrixR, const TGeoMatrix &matrixL);
TGeoHMatrix     operator*(TGeoHMatrix &matrixR, const TGeoHMatrix &matrixL);
TGeoHMatrix     operator*(const TGeoMatrix &matrixR, const TGeoMatrix &matrixL);
TGeoHMatrix     operator*(const TGeoMatrix &matrixR, const TGeoHMatrix &matrixL);
TGeoTranslation operator+(const TGeoTranslation &matrixR, const TGeoTranslation &matrixL);
TGeoTranslation operator-(const TGeoTranslation matrixR, const TGeoTranslation &matrixL);
TGeoTranslation operator+(const TGeoTranslation matrixR, const TGeoMatrix &matrixL);
TGeoTranslation operator*(const TGeoTranslation matrixR, const Double_t &scale);
void            addVectorToElement(SGeoBody &body, unsigned int position, TGeoTranslation vector);

//-----------------------------------------------------------------------------------------------------------------------------
// TShieldGeometryConvert.cxx
std::pair<zoneList, zoneList> getZoneFromBody(TGeoNode *node, TGeoHMatrix oldTransformation,
                                              double scale = getDefaultScale());
//  WARNING: extended initializer lists only available with -std=c++11 or -std=gnu++11
std::pair<zoneList, zoneList>        getZoneFromShape(TGeoShape *shape, TGeoHMatrix currTransformation,
                                                      double scale = getDefaultScale());
inline std::pair<zoneList, zoneList> boxToZones(TGeoBBox *box, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> tubeSegToZones(TGeoTubeSeg *tube, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> tubeToZones(TGeoTube *tube, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> coneSegToZones(TGeoConeSeg *cone, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> coneToZones(TGeoCone *cone, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> trd1ToZones(TGeoTrd1 *trd, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> trd2ToZones(TGeoTrd2 *trd, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> trapToZones(TGeoTrap *trap, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> sphereToZones(TGeoSphere *sphere, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> polyconeToZones(TGeoPcon *polycone, TGeoHMatrix currTransformation, double scale);
inline std::pair<zoneList, zoneList> polyhedraToZones(TGeoPgon *polyhedra, TGeoHMatrix currTransformation,
                                                      double scale);
inline std::pair<zoneList, zoneList> ellipticalTubeToZones(TGeoEltu *tube, TGeoHMatrix currTransformation,
                                                           double scale);
inline std::pair<zoneList, zoneList> arb8ToZones(TGeoArb8 *gtrap, TGeoHMatrix currTransformation, double scale);
std::pair<zoneList, zoneList>        compositeShapeToZones(TGeoCompositeShape *shape, TGeoHMatrix currTransformation,
                                                           double scale);
MediumData                           getMedium(TGeoNode *node);

} // end of namespace tgeanttoshield

// Obsolete, full converting of geometry;
Geant4ShieldData convertTGeantToShield(TGeoManager *geoMan);
namespace tgeanttoshield {
zoneData                      convertBodyRecursive(TGeoNode *node);
std::pair<zoneData, zoneList> convertBodyRecursive(TGeoNode *node, TGeoHMatrix oldTransformation);
} // namespace tgeanttoshield

#endif