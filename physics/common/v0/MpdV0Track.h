/*
 * MpdV0Track.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0TRACK_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0TRACK_H_

#include "MpdV0Particle.h"

class MpdV0Track : public MpdV0Particle {
   Double_t fFirstDaughterS;
   Double_t fSecondDaughterS;
   Double_t fChi2;

public:
   MpdV0Track();
   MpdV0Track(const MpdV0Track &other) = default;
   MpdV0Track &operator=(const MpdV0Track &other) = default;

   Double_t GetFirstDaughterS() const { return fFirstDaughterS; }
   Double_t GetSecondDaughterS() const { return fSecondDaughterS; }
   Double_t GetChi2() const { return fChi2; };

   void SetFirstDaughterS(Double_t firstDaughterS) { this->fFirstDaughterS = firstDaughterS; }
   void SetSecondDaughterS(Double_t secondDaughterS) { this->fSecondDaughterS = secondDaughterS; }
   void SetChi2(Double_t chi) { fChi2 = chi; };
   virtual ~MpdV0Track();
   ClassDef(MpdV0Track, 2)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0TRACK_H_ */
