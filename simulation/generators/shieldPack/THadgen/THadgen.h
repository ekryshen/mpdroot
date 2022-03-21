// -------------------------------------------------------------------------
// -----          THadgen header file                                  -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

/**  THadgen.h
 *@author Alexander Timofeev (12.2011) <antimofeew@gmail.com>
 *@author Dmitry Sosnov (05.2015) <dmitry.e.sosnov@gmail.com>
 *
 * The THadgen is an interface (wrapper) class
 * for HADGEN generator. Original code for HADGEN (SHIELD)
 * is written at JINR (Dubna) and INR RAS (Moscow).
 * Derived from TGenerator.
 **/

#ifndef THadgen_H
#define THadgen_H

#include "hadgen.h"

#include <cmath>
#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TGenerator.h"
#include "TParticle.h"

class THadgen : public TGenerator {
   ClassDef(THadgen, 1) // THadgen
      protected :
      // auto-cleanup mechanism
      static THadgen *fInstance;
   class THadgenClean {
   public:
      THadgenClean() {}
      ~THadgenClean();
   };
   friend class THadgenClean;

   Int_t fParticlesNum;
   bool  fAutoSeed;

private:
   // THadgen class should be never copied
   THadgen(const THadgen &s) : TGenerator(s) {}
   THadgen &operator=(const THadgen &) { return *this; }

public:
   THadgen(Int_t autoseed = true); // creates a unique instance of THadgen
   virtual ~THadgen();

   static THadgen *Instance();

   // Random generator state
   void  SetRandomSeed(Int_t seed);
   void  SetRandomSeed();
   Int_t GetRandomSeed();

   void SetAutoSeed(bool flag = true); // initially true

protected:
   // the following functions should not be used
   void  SetLuxCount(Int_t luxcnt);
   Int_t GetLuxCount();

public:
   // Inreac properties
   void    SetIncidentParticle(Int_t jpart); // these functions may be more smart
   Int_t   GetIncidentParticle();
   void    SetNuclid(Int_t nuclid);
   Int_t   GetNuclid();
   void    SetEnergy(Float_t energy);
   Float_t GetEnergy();
   // Flags
   void  SetSystem(Int_t lanti = 0);
   Int_t GetSystem();
   void  SetLAntil(Int_t lantil = 0);
   Int_t GetLAntil();
   // Debug options (prInt_ting options)
   // 0 means no output, 1 - force output
   void  SetLSTAR(Int_t lstar = 0);
   Int_t GetLSTAR();
   void  SetLCASC(Int_t lcasc = 0);
   Int_t GetLCASC();

   // Monte-Carlo statistics parameter
   void  SetNumStat(Int_t nstat = 1);
   Int_t GetNumStat();

   // A and Z of projectile ion
   // notice that these params are set automaticaly for
   // Deuteron, Tritium, He3 and alpha particles
   void    SetAProj(Float_t aproj);
   Float_t GetAProj();
   void    SetZProj(Float_t zproj);
   Float_t GetZProj();
   void    SetParticleFromPdgCode(Int_t pdg_code);
   void    SetParticleFromPdgCode(Int_t pdg_code, Float_t A, Float_t Z);

   // Parameter functions derived from TGenerator
   // virtual void SetParameter(const char *name, Double_t value);
   // virtual Double_t GetParameter(const char *name);

   // Generation Functions derived from TGenerator
   void       GenerateEvent();
   void       CopyToInternalArray();
   TObjArray *ImportParticles(Option_t *option = "");
   Int_t      ImportParticles(TClonesArray *particles, Option_t *option = "");

   void FileOut(char *filename); // test function
};

#endif
