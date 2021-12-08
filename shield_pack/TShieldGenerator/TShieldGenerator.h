// -------------------------------------------------------------------------
// -----          TShieldGenerator header file                         -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

/**  TShieldGenerator.h
 *@author Dmitry Sosnov (05.2015) <dmitry.e.sosnov@gmail.com>
 *
 * The TShieldGenerator is a generator for using TShield
 * after THadgen collisions.
 * Derived from TGenerator.
 **/

#ifndef TSHIELDGENERATOR_H
#define TSHIELDGENERATOR_H

#include "THadgen.h"
#include "TShield.h"
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TGenerator.h"
#include "TParticle.h"
#include "TGeoManager.h"

class TShieldGenerator : public TGenerator {
   ClassDef(TShieldGenerator, 1) // TShieldGenerator
      protected :
      // auto-cleanup mechanism
      static TShieldGenerator *fInstance;
   class TShieldGeneratorClean {
   public:
      TShieldGeneratorClean() {}
      ~TShieldGeneratorClean();
   };
   friend class TShieldGeneratorClean;

   Int_t fParticlesNum;
   bool  autoseed;

private:
   // TShieldGenerator class should be never copied
   TShieldGenerator(const TShieldGenerator &s) : TGenerator(s) {}
   TShieldGenerator &operator=(const TShieldGenerator &) { return *this; }

public:
   //         TShieldGenerator(const TGenerator &tg);
   TShieldGenerator();
   ~TShieldGenerator();

   void SetGeometry(TGeoManager *geom);
   void SetEnergy(Float_t energy);
   void SetStartPoint(Float_t x, Float_t y, Float_t z);
   void SetParticleFromPdgCode(Int_t pdg_code, Float_t A, Float_t Z);
   void SetStatistics(Int_t nstat);

   void       GenerateEvent();
   TObjArray *ImportParticles(Option_t *option = "");
   Int_t      ImportParticles(TClonesArray *particles, Option_t *option = "");

private:
   TShield     *fShield;
   TGeoManager *fGeoManager;
   Double_t     fEnergy;
   Double_t     fStartX, fStartY, fStartZ;
   Int_t        fPdgCode, fA, fZ, fNumStat;
};

#endif
