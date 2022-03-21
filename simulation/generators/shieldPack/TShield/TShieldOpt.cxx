// -------------------------------------------------------------------------
// -----          TShield source file                                  -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

#include "TShield.h"

void TShield::SetAutoSeed(bool flag)
{
   autoseed = flag;
}

void TShield::SetRandomSeed()
{
   SetRandomSeed(clock() + time(0));
}

void TShield::SetRandomSeed(Int_t seed)
{
   if (seed < 0) seed = -seed;
   shield_set_randomseed(seed);
}

Int_t TShield::GetRandomSeed()
{
   return shield_get_randomseed();
}

void TShield::SetLuxCount(Int_t luxcnt)
{
   shield_set_luxcount(luxcnt);
}

Int_t TShield::GetLuxCount()
{
   return shield_get_luxcount();
}

void TShield::SetLCASC(Int_t lcasc)
{
   shield_set_lcasc(lcasc);
}

Int_t TShield::GetLCASC()
{
   return shield_get_lcasc();
}

void TShield::SetIncidentParticle(Int_t jpart)
{
   shield_set_incidentparticle(jpart);
}

Int_t TShield::GetIncidentParticle()
{
   return shield_get_incidentparticle();
}

void TShield::SetEnergy(Float_t energy)
{
   shield_set_energy(energy);
}

Float_t TShield::GetEnergy()
{
   return shield_get_energy();
}

void SetStatistics(Int_t nstat)
{
   shield_set_statistics(nstat);
}
Int_t GetStatistics()
{
   return shield_get_statistics();
}

void TShield::SetAProj(Float_t aproj)
{
   shield_set_aproj(aproj);
}

void TShield::SetZProj(Float_t zproj)
{
   shield_set_zproj(zproj);
}

Float_t TShield::GetAProj()
{
   return shield_get_aproj();
}

Float_t TShield::GetZProj()
{
   return shield_get_zproj();
}

void TShield::SetStartPoint(Float_t x, Float_t y, Float_t z)
{
   shield_set_runpoint(x, y, z);
}

void TShield::SetDirectionVector(Float_t cx, Float_t cy, Float_t cz)
{
   Float_t norm = sqrt(cx * cx + cy * cy + cz * cz);
   cx /= norm;
   cy /= norm;
   cz /= norm;
   Float_t sin_theta = sqrt(1 - cz * cz);
   if (sin_theta != 0.0) {
      shield_set_rundirection(cz, cy / sin_theta, cx / sin_theta);
   } else {
      shield_set_rundirection(1, 0, 1);
   }
}

void TShield::SetDirection(Float_t cos_theta, Float_t sin_phi, Float_t cos_phi)
{
   shield_set_rundirection(cos_theta, sin_phi, cos_phi);
}

void TShield::SetDirection(Float_t phi, Float_t theta)
{
   shield_set_rundirection(cos(theta), sin(phi), cos(phi));
}

void TShield::SetParticleFromPdgCode(Int_t pdg_code)
{
   float a, z;
   int   type = hadgen_get_nuclei_parameters_by_pdg_code(pdg_code, &a, &z);
   if (type == 0) {
      printf("TShield: Unknown pdg code for beam");
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
void TShield::SetParticleFromPdgCode(Int_t pdg_code, Float_t A, Float_t Z)
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