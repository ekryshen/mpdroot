// -------------------------------------------------------------------------
// -----          SHIELD header file                                   -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

/**  shield_media.h
 *@author Alexander Timofeev (02.2012) <antimofeew@gmail.com>
 *@author Dmitry Sosnov (02.2015) <dmitry.e.sosnov@gmail.com>
 *
 * The SHIELD library header.
 * The SHIELD library used to connect SHIELD generator with
 * CERN Root system (class TShield).
 * Original code for SHIELD is written at JINR (Dubna) and INR RAS (Moscow).
 **/

#ifndef SHIELD_MEDIA_H
#define SHIELD_MEDIA_H

#ifndef SHIELD_H
#error ERROR: Don`t include shield_media.h header, include shield.h instead
#endif

struct Element {
   float Nuclid;
   float Conc;        // Concentration [10^27/cm^3]
   float Density;     // Partial density [g/cm^3]
   float Z;           // Atomic number Z
   float A;           // Atomic weight A
   float PureDensity; // Density of pure material [g/cm^3]
   float ionEv;       // Ionization potential [eV]
   float reserved;    // reserved
};

struct MediumData {
   int            nType;        // Medium type, should be in range 1-4
   int            nChemEl;      // up to 24 elements
   float          Rho;          // Density of medium
   struct Element Elements[24]; // elements info
};

void            shield_send_medium(struct MediumData medium, int MediaNum);
int             _isElementsEqual(struct Element element1, struct Element element2);
int             _isMediumEqual(struct MediumData medium1, struct MediumData medium2);
int             shield_add_medium(struct MediumData medium);
SHIELD_API void shield_init_medium(void);
void            shield_clean_medium(void);

#endif
