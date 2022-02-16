/*
 * MpdV0StandardCandidateMonitor.cxx
 *
 *  Created on: 11 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0StandardCandidateMonitor.h"

#include <RtypesCore.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <vector>

#include "MpdV0Particle.h"

MpdV0StandardCandidateMonitor::MpdV0StandardCandidateMonitor(MpdV0::EParticleType type)
   : MpdV0CandidateMonitor(4, 1), fType(type)
{
   SetCosAxis(100);
   SetArmenterosAlphaAxis(100, -1, 1);
   SetArmenterosPtAxis(100, 0, 1);
   SetMinvAxis(100, 0.5, 1.5);

   SetDecayLenghtMonitor(100, 0, 10);
   SetDau1to2(100, 0, 10);
}

void MpdV0StandardCandidateMonitor::Fill(const MpdV0Particle &particle, Bool_t status)
{
   Double_t mass = 0;
   switch (fType) {
   case MpdV0::EParticleType::k0Short: {
      mass = particle.GetK0Mass();
   } break;
   case MpdV0::EParticleType::kAntiLambda: {
      mass = particle.GetAntiLambdaMass();
   } break;
   default: {
      mass = particle.GetLambdaMass();
   } break;
   }

   Fill1D(0, mass, status);
   Fill1D(1, particle.GetCosAngle(), status);
   Fill1D(2, particle.GetDecayLenght(), status);
   Fill1D(3, particle.GetDau1to2(), status);
   Fill2D(0, particle.GetAplhaArm(), particle.GetPtArm(), status);
}

void MpdV0StandardCandidateMonitor::Init()
{
   MakeHistogram1d("Minv", "m_{V0} [GeV/c^2]", 0);
   MakeHistogram1d("Cos", "cos [rad]", 1);
   MakeHistogram1d("DL", "DecayLenght [cm]", 2);
   MakeHistogram1d("Dau12", "Dau1to2 [cm]", 3);
   MakeHistogram2d("Armenteros", "#alpha", "p_T", 0);
}

MpdV0StandardCandidateMonitor::~MpdV0StandardCandidateMonitor() {}
