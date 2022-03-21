// -------------------------------------------------------------------------
// -----          SHIELD header file                                   -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

/**  shield_geo.h
 *@author Alexander Timofeev (02.2012) <antimofeew@gmail.com>
 *@author Dmitry Sosnov (02.2015) <dmitry.e.sosnov@gmail.com>
 *
 * The SHIELD library header.
 * The SHIELD library used to connect SHIELD generator with
 * CERN Root system (class TShield).
 * Original code for SHIELD is written at JINR (Dubna) and INR RAS (Moscow).
 **/

#ifndef SHIELD_H
#error ERROR: Don`t include shield_geo.h, include shield.h instead.
#endif

#ifndef SHIELD_GEO_H
#define SHIELD_GEO_H

// Geometry callbacks
typedef double (*GCURZLCALLBACK)(double x, double y, double z, double vx, double vy, double vz);
typedef int (*GNEXTZCALLBACK)(double x, double y, double z, double vx, double vy, double vz);

struct SGeoBody {
   int    type; // Warning! Type is wrote before parameters in bodydb
   double parameters[36];
};
struct SGeoZone {
   int countELements;
   int mediaNumber;
   int definition[10000]; // Current count of elements at fortran ZONEDB is 5000, but we hawe number of body and type of
                          // zone separately.
};

SHIELD_API void shield_set_gcurzl_callback(GCURZLCALLBACK function);
SHIELD_API void shield_set_gnextz_callback(GNEXTZCALLBACK function);

SHIELD_API int shield_add_body(int type, double *parameters);
SHIELD_API int shield_add_zone(int countBodies, int *zoneParameters, struct SGeoBody *bodies, struct MediumData medium);
SHIELD_API void shield_clean_geometry();
SHIELD_API void shield_init_geometry();

#endif
