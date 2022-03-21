// -------------------------------------------------------------------------
// -----          HADGEN source file                                   -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

#define HADGEN_LIB_INTERNAL
#include "hadgen.h"
#include "hadgen_common_blocks.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

const char *hadgen_get_particle_name(int type)
{
   static char _name[20];
   switch (type) {
   case PARTICLE_Neutron: strcpy(_name, "Neutron"); break;
   case PARTICLE_AntiNeutron: strcpy(_name, "Anti-neutron"); break;
   case PARTICLE_AntiProton: strcpy(_name, "Anti-proton"); break;
   case PARTICLE_Proton: strcpy(_name, "Proton"); break;
   case PARTICLE_PIminus: strcpy(_name, "Pi-"); break;
   case PARTICLE_PIplus: strcpy(_name, "Pi+"); break;
   case PARTICLE_PI0: strcpy(_name, "Pi0"); break;
   case PARTICLE_Kminus: strcpy(_name, "K-"); break;
   case PARTICLE_Kplus: strcpy(_name, "K+"); break;
   case PARTICLE_K0Long: strcpy(_name, "K0 Long-lived"); break;
   case PARTICLE_K0Short: strcpy(_name, "K0 Short-lived"); break;
   case PARTICLE_Gamma: strcpy(_name, "Gamma"); break;
   case PARTICLE_Deuterium: strcpy(_name, "Deuterium"); break;
   case PARTICLE_Tritium: strcpy(_name, "Tritium"); break;
   case PARTICLE_He3: strcpy(_name, "He-3"); break;
   case PARTICLE_Alpha: strcpy(_name, "Alpha"); break;
   case PARTICLE_NONE:
   default: strcpy(_name, "Unknown");
   }
   return _name;
}

int hadgen_get_pdg_code(int particle_type)
{
   //   printf("pdg: jpart = %i", particle_type);
   switch (particle_type) {
   case PARTICLE_Neutron: return 2112;
   case PARTICLE_AntiNeutron: return -2112;
   case PARTICLE_AntiProton: return -2212;
   case PARTICLE_Proton: return 2212;
   case PARTICLE_PIminus: return -211;
   case PARTICLE_PIplus: return 211;
   case PARTICLE_PI0: return 111;
   case PARTICLE_Kminus: return -321;
   case PARTICLE_Kplus: return 321;
   case PARTICLE_K0Long: return 130; // Long-lived only
   case PARTICLE_K0Short: return 311;
   case PARTICLE_Gamma: return 22;
   case PARTICLE_Deuterium: return hadgen_get_pdg_code_nuclei(2, 1);
   case PARTICLE_Tritium: return hadgen_get_pdg_code_nuclei(3, 1);
   case PARTICLE_He3: return hadgen_get_pdg_code_nuclei(3, 2);
   case PARTICLE_Alpha: return hadgen_get_pdg_code_nuclei(4, 2);
   case PARTICLE_Electron: return 11;
   case PARTICLE_Positron: return -11;
   case PARTICLE_MUminus: return 13;
   case PARTICLE_MUplus: return -13;
   case PARTICLE_ElNeutrino: return 12;
   case PARTICLE_ElAntiNeutrino: return -12;
   case PARTICLE_MuNeutrino: return 14;
   case PARTICLE_MuAntiNeutrino: return -14;
   case PARTICLE_NONE:
   default: return 0;
   }
}
int hadgen_get_particle_type(int pdg_code)
{
   const int particlePDGs[26] = {0,
                                 2112,
                                 2212,
                                 -211,
                                 211,
                                 111,
                                 -2112,
                                 -2212,
                                 -321,
                                 321,
                                 130,
                                 -130,
                                 22,
                                 11,
                                 -11,
                                 13,
                                 -13,
                                 12,
                                 -12,
                                 14,
                                 -14,
                                 hadgen_get_pdg_code_nuclei(2, 1),
                                 hadgen_get_pdg_code_nuclei(3, 1),
                                 hadgen_get_pdg_code_nuclei(3, 2),
                                 hadgen_get_pdg_code_nuclei(4, 2),
                                 0};
   int       i                = 0;
   for (i = 0; i < 25; i++) {
      if (particlePDGs[i] == pdg_code) {
         return i;
      }
   }
   return 25;
}

const int Z_CODE      = 10000;
const int A_CODE      = 10;
const int NUCLEI_CODE = 1000000000;

int hadgen_get_pdg_code_nuclei(float _A, float _Z)
{
   //    printf("pdg for %3.0f %3.0f", _A, _Z);
   int A = _A, Z = _Z;
   // check for ambiguities
   if ((A == 1) && (Z == 1)) return hadgen_get_pdg_code(PARTICLE_Proton);
   if (A > 250) return -99901; // Why?
   if (Z > 125) return -99902;
   if (Z < -125) return -99903;
   if (A <= 0) return -99904;
   if (A < abs(Z)) return -99905; // Weight number should be greater than charge number
   int result = NUCLEI_CODE;
   if (Z > 0)
      result += A_CODE * A + Z_CODE * Z;
   else
      result += A_CODE * A - Z_CODE * Z;
   if (Z < 0) result = -result; // for antinuclei
   return result;
}

int hadgen_get_nuclei_parameters_by_pdg_code(int pdg_code, float *_A, float *_Z)
{
   int type = hadgen_get_particle_type(pdg_code);
   if (type != 25) {
      *_A = 0;
      *_Z = 0;
      return type;
   }
   int result           = NUCLEI_CODE;
   int antimatterFactor = (pdg_code > 0) ? 1 : -1;
   switch (pdg_code) {
   case -99901:
   case -99902:
   case -99903:
   case -99904:
   case -99905:
      type = 0;
      *_A  = 0;
      *_Z  = 0;
      break;
   default:
      result = abs(pdg_code) - result;
      if (Z_CODE > A_CODE) {
         *_Z = antimatterFactor * ceil(result / Z_CODE);
         *_A = ceil((result % Z_CODE) / A_CODE);
      } else {
         *_A = ceil(result / A_CODE);
         *_Z = antimatterFactor * ceil((result % A_CODE) / Z_CODE);
      }
      break;
   }
   return type;
}

void hadgen_set_randomseed(int seed)
{
   random.IXFIRS = seed;
}
void hadgen_set_luxcount(int luxcnt)
{
   other.LUXCNT = luxcnt;
}
void hadgen_set_incidentparticle(int jpart)
{
   // check for validity TODO
   inreac.JPART = jpart;
}
void hadgen_set_nuclid(int nuclid)
{
   // check for validity TODO
   inreac.NUCLID = nuclid;
}
void hadgen_set_energy(float energy)
{ // input energy in GeV
   if (energy > 0)
      inreac.TINT = energy * 1000.;
   else
      inreac.TINT = -energy * 1000.;
}
void hadgen_set_system(int lanti)
{
   if (lanti != 1) lanti = 0;
   antlab.lanti = lanti;
}
void hadgen_set_lantil(int lantil)
{
   antil.lantil = lantil;
}
void hadgen_set_lstar(int lstar)
{
   if (lstar != 1) lstar = 0;
   debug.LSTAR = lstar;
}
void hadgen_set_lcasc(int lcasc)
{
   if (lcasc != 1) lcasc = 0;
   debug.LCASC = lcasc;
}
void hadgen_set_statisticsnum(int nstat)
{
   if (nstat < 0) nstat = -nstat;
   if (nstat == 0) nstat = 1;
   other.NSTAT = nstat;
}
void hadgen_set_aprojectile(float aproj)
{
   // Check for validity will be performed later
   hiproj.APROJ = aproj;
}
void hadgen_set_zprojectile(float zproj)
{
   hiproj.ZPROJ = zproj;
}

int hadgen_get_randomseed()
{
   return random.IXFIRS;
}
int hadgen_get_luxcount()
{
   return other.LUXCNT;
}
int hadgen_get_incidentparticle()
{
   return inreac.JPART;
}
int hadgen_get_nuclid()
{
   return inreac.NUCLID;
}
float hadgen_get_energy()
{
   return inreac.TINT * 0.001; // output energy in GeV
}
int hadgen_get_system()
{
   return antlab.lanti;
}
int hadgen_get_lantil()
{
   return antil.lantil;
}
int hadgen_get_LSTAR()
{
   return debug.LSTAR;
}
int hadgen_get_LCASC()
{
   return debug.LCASC;
}
int hadgen_get_statisticsnum()
{
   return other.NSTAT;
}
float hadgen_get_aprojectile()
{
   return hiproj.APROJ;
}
float hadgen_get_zprojectile()
{
   return hiproj.ZPROJ;
}
