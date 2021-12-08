// -------------------------------------------------------------------------
// -----          TShield header file                                  -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

/**  TShield.h
 *@author Alexander Timofeev (02.2012) <antimofeew@gmail.com>
 *@author Dmitry Sosnov (05.2015) <dmitry.e.sosnov@gmail.com>
 *
 * The TShield is an interface (wrapper) class
 * for SHIELD transport code. Original code for SHIELD
 * is written at JINR (Dubna) and INR RAS (Moscow).
 * Derived from TGenerator.
 **/

#ifndef TSHIELD_H
#define TSHIELD_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "TObjArray.h"
#include "TClonesArray.h"
#include "TGenerator.h"
#include "TParticle.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"

#include "shield.h"
#include "hadgen.h"
#include "TShieldGeometry.h"

class TShield {
   ClassDef(TShield, 1) // TShield

      protected :
      // auto-cleanup mechanism
      static TShield *fInstance;
   class TShieldClean {
   public:
      TShieldClean() {}
      ~TShieldClean();
   };
   friend class TShieldClean;
   TGeoManager *fGeometry;

   TClonesArray *fParticlesFlyOut, *fParticlesAbsorbed;
   Int_t         ParticlesNumFlyOut, ParticlesNumAbsorbed;
   int           ParticlesNumFlyOutMax, ParticlesNumAbsorbedMax;
   bool          clearAtStart;
   bool          autoseed;
   TTree        *fTree;
   TParticle    *fCurrentParticle;

public:
   const TClonesArray *GetFlyOutArray() { return fParticlesFlyOut; }
   const TClonesArray *GetAbsorbedArray() { return fParticlesAbsorbed; }
   const TTree        *GetTree() { return fTree; }
   Int_t               GetFlyOutNum() { return ParticlesNumFlyOut; }
   Int_t               GetAbsorbedNum() { return ParticlesNumAbsorbed; }
   void                ClearArrays();

private:
   // TShield class should be never copied
   TShield(const TShield &s) {}
   TShield &operator=(const TShield &) { return *this; }

   void InitCallbacks();

   // static callbacks
   static int    GeoNextCallback(double, double, double, double, double, double);
   static double GeoDistCallback(double, double, double, double, double, double);
   static void   TreeCallback(shield_tree_node *);

   void   AddTreeParticle(shield_tree_node *);     // non-static callback to store the tree
   void   AddFliedOutParticle(shield_tree_node *); // non-static callback to store flied-out particles
   void   AddAbsorbedParticle(shield_tree_node *); // non-static callback to store absorbed particles
   double FindBoundary(double x, double y, double z, double vx, double vy,
                       double vz); // non-static callbacks to track geometry
   int    GetNextVolume(double x, double y, double z, double vx, double vy, double vz);

protected:
   // the following functions should not be used
   void  SetLuxCount(Int_t luxcnt);
   Int_t GetLuxCount();
   void  SetLCASC(Int_t lcasc = 0);
   Int_t GetLCASC();

public:
   TShield(bool clean = kTRUE, int FlyOutMax = 10000, int AbsorbedMax = 10000); // creates a unique instance of TShield
   virtual ~TShield();

   static TShield *Instance();

   // Random generator state
   void  SetRandomSeed(Int_t seed);
   void  SetRandomSeed();
   Int_t GetRandomSeed();

   void SetAutoSeed(bool flag = false); // initially false

   void         SetGeometry(TGeoManager *geom);
   TGeoManager *GetGeometry();

   // Inreac properties
   void  SetIncidentParticle(Int_t jpart); // these functions may be more smart
   Int_t GetIncidentParticle();

   // void SetNuclid(Int_t nuclid);
   // Int_t GetNuclid();

   void    SetEnergy(Float_t energy);
   Float_t GetEnergy();
   void    SetStatistics(Int_t nstat);
   Int_t   GetStatistics();

   // A and Z of projectile ion
   // notice that these params are set automaticaly for
   // Deuteron, Tritium, He3 and alpha particles
   void    SetAProj(Float_t aproj);
   Float_t GetAProj();
   void    SetZProj(Float_t zproj);
   Float_t GetZProj();
   void    SetParticleFromPdgCode(Int_t pdg_code);
   void    SetParticleFromPdgCode(Int_t pdg_code, Float_t A, Float_t Z);

   void SetStartPoint(Float_t x, Float_t y, Float_t z);
   void SetDirection(Float_t cos_theta, Float_t sin_phi, Float_t cos_phi);
   void SetDirection(Float_t phi, Float_t theta);
   void SetDirectionVector(Float_t fx, Float_t cy, Float_t cz);
   void GenerateEvent();

   void PrintGeometry();
};

#endif
