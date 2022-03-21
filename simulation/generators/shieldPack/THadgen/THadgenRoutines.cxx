// -------------------------------------------------------------------------
// -----          THadgen source file                                  -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

#include "THadgen.h"

// setters in THadgen class
//______________________________________________________________________________
void THadgen::SetRandomSeed(Int_t seed)
{
   if (seed < 0) seed = -seed; // should be positive!
   hadgen_set_randomseed(seed);
}
//______________________________________________________________________________
void THadgen::SetRandomSeed()
{
   SetRandomSeed(clock() + time(0));
}
//______________________________________________________________________________
void THadgen::SetAutoSeed(bool flag)
{
   fAutoSeed = flag;
}
//______________________________________________________________________________
void THadgen::SetLuxCount(Int_t luxcnt)
{
   hadgen_set_luxcount(luxcnt);
}
//______________________________________________________________________________
void THadgen::SetIncidentParticle(Int_t jpart)
{
   hadgen_set_incidentparticle(jpart);
}
//______________________________________________________________________________
void THadgen::SetNuclid(Int_t nuclid)
{
   hadgen_set_nuclid(nuclid);
}
//______________________________________________________________________________
void THadgen::SetEnergy(Float_t energy)
{
   hadgen_set_energy(energy);
}
//______________________________________________________________________________
void THadgen::SetSystem(Int_t lanti)
{
   hadgen_set_system(lanti);
}
//______________________________________________________________________________
void THadgenSetLAntil(Int_t lantil)
{
   hadgen_set_lantil(lantil);
}
//______________________________________________________________________________
void THadgen::SetLSTAR(Int_t lstar)
{
   hadgen_set_lstar(lstar);
}
//______________________________________________________________________________
void THadgen::SetLCASC(Int_t lcasc)
{
   hadgen_set_lcasc(lcasc);
}
//______________________________________________________________________________
void THadgen::SetNumStat(Int_t nstat)
{
   hadgen_set_statisticsnum(nstat);
}
//______________________________________________________________________________
void THadgen::SetAProj(Float_t aproj)
{
   hadgen_set_aprojectile(aproj);
}
//______________________________________________________________________________
void THadgen::SetZProj(Float_t zproj)
{
   hadgen_set_zprojectile(zproj);
}

void THadgen::SetParticleFromPdgCode(Int_t pdg_code)
{
   float a, z;
   int   type = hadgen_get_nuclei_parameters_by_pdg_code(pdg_code, &a, &z);
   if (type == 0) {
      printf("THadgen: Unknown pdg code for beam");
      return;
      //  } else if(type != 25){ //Disabled, because a and z should be null
      //         SetIncidentParticle(type); SetAProj(0); SetZProj(0);
      //     }
   } else {
      SetIncidentParticle(type);
      SetAProj(a);
      SetZProj(z);
   }
   return;
}
void THadgen::SetParticleFromPdgCode(Int_t pdg_code, Float_t A, Float_t Z)
{
   float _a, _z;
   int   type = hadgen_get_nuclei_parameters_by_pdg_code(pdg_code, &_a, &_z);
   if (type == 25) {
      SetIncidentParticle(type);
      SetAProj(_a);
      SetZProj(_z);
   } else if (type != 0) {
      SetIncidentParticle(type);
      SetAProj(_a);
      SetZProj(_z);
   } else { // if(type == 0)
      SetIncidentParticle(25);
      SetAProj(A);
      SetZProj(Z);
   }
   return;
}

// getters in THadgen
//______________________________________________________________________________
Int_t THadgen::GetRandomSeed()
{
   return hadgen_get_randomseed();
}
//______________________________________________________________________________
Int_t THadgen::GetLuxCount()
{
   return hadgen_get_luxcount();
}
//______________________________________________________________________________
Int_t THadgen::GetIncidentParticle()
{
   return hadgen_get_incidentparticle();
}
//______________________________________________________________________________
Int_t THadgen::GetNuclid()
{
   return hadgen_get_nuclid();
}
//______________________________________________________________________________
Float_t THadgen::GetEnergy()
{
   return hadgen_get_energy();
}
//______________________________________________________________________________
Int_t THadgen::GetSystem()
{
   return hadgen_get_system();
}
//______________________________________________________________________________
Int_t THadgen::GetLAntil()
{
   return hadgen_get_lantil();
}
//______________________________________________________________________________
Int_t THadgen::GetLSTAR()
{
   return hadgen_get_lstar();
}
//______________________________________________________________________________
Int_t THadgen::GetLCASC()
{
   return hadgen_get_lcasc();
}
//______________________________________________________________________________
Int_t THadgen::GetNumStat()
{
   return hadgen_get_statisticsnum();
}
//______________________________________________________________________________
Float_t THadgen::GetAProj()
{
   return hadgen_get_aprojectile();
}
//______________________________________________________________________________
Float_t THadgen::GetZProj()
{
   return hadgen_get_zprojectile();
}
