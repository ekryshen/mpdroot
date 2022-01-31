/*
 * MpdV0Particle.cxx
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0Particle.h"

#include <TMath.h>
#include <iostream>

MpdV0Particle::MpdV0Particle()
   : fDau1to2(0), fAplhaArm(0), fPtArm(0), fCosAngle(0), fDecLenght(0), fPosDaugherIndex(-1),
     fNegDaughterIndex(-1), fPdg(0)
{
}

MpdV0Particle::~MpdV0Particle() {}

void MpdV0Particle::Recalculate(const TVector3 &vertex)
{
   fMomentum     = fMomPosDaughter + fMomNegDaughter;
   Double_t pTot = fMomentum.Mag();

   Double_t pPosTot = fMomPosDaughter.Mag();

   Double_t MomPosAlongV0 = fMomPosDaughter * fMomentum / pTot;
   Double_t MomNegALongV0 = fMomNegDaughter * fMomentum / pTot;
   /* std::cout << "****" << std::endl;
    std::cout << MomPosAlongV0 << " " << MomNegALongV0 << " " << fMomFirstDaughter.Mag() << " "
              << fMomSecondDaughter.Mag() << std::endl;
    TVector3 a      = fMomFirstDaughter;
    TVector3 b      = fMomSecondDaughter;
    TVector3 c      = fMomentum;
    auto     print  = [](TVector3 &v) { std::cout << Form("%4.4f %4.4f %4.4f", v.X(), v.Y(), v.Z()) << std::endl; };
    auto     print2 = [](TVector3 &v, TVector3 &v2) {
       std::cout << Form("%4.4f %4.4f %4.4f\t%4.4f", v.X() * v2.X(), v.Y() * v2.Y(), v.Z() * v2.Z(), v * v2)
                 << std::endl;
    };
    print(a);
    print(b);
    print(c);
    print2(a, c);
    print2(b, c);*/
   SetAlphaArm((MomPosAlongV0 - MomNegALongV0) / (MomPosAlongV0 + MomNegALongV0));
   SetPtArm(TMath::Sqrt(pPosTot * pPosTot - MomPosAlongV0 * MomPosAlongV0));

   TVector3 pozV0 = fDecayPoint - vertex;

   Double_t t = -(pozV0 * fMomentum) / (pTot * pTot);
   fDca.SetXYZ(pozV0.X() + t * fMomentum.X(), pozV0.Y() + t * fMomentum.Y(), pozV0.Z() + t * fMomentum.Z());
   TVector3 dca_rel = fDca - pozV0;
   fCosAngle        = pozV0 * fMomentum / (pTot * pozV0.Mag());
   fDecLenght       = dca_rel.Mag();
}

const Double_t MpdV0Particle::GetMinv(MpdV0::EParticleType type) const
{
   Double_t ps   = fMomentum.Mag2();
   Double_t p1_2 = fMomPosDaughter.Mag2();
   Double_t p2_2 = fMomNegDaughter.Mag2();
   Double_t e1, e2;
   switch (type) {
   case MpdV0::EParticleType::k0Short: {
      e1 = TMath::Sqrt(p1_2 + 0.019479785);
      e2 = TMath::Sqrt(p2_2 + 0.019479785);
      return 0.49761;
   } break;
   case MpdV0::EParticleType::kLambda: {
      e1 = TMath::Sqrt(p1_2 + 0.880354511);
      e2 = TMath::Sqrt(p2_2 + 0.019479785);
   } break;
   case MpdV0::EParticleType::kAntiLambda: {
      e1 = TMath::Sqrt(p1_2 + 0.019479785);
      e2 = TMath::Sqrt(p2_2 + 0.880354511);
   } break;
   case MpdV0::EParticleType::kPdgHypo: {
      switch (fPdg) {
      case 310: {
         e1 = TMath::Sqrt(p1_2 + 0.019479785);
         e2 = TMath::Sqrt(p2_2 + 0.019479785);
      } break;
      case 3122: {
         e1 = TMath::Sqrt(p1_2 + 0.880354511);
         e2 = TMath::Sqrt(p2_2 + 0.019479785);
      } break;
      case -3122: {
         e1 = TMath::Sqrt(p1_2 + 0.019479785);
         e2 = TMath::Sqrt(p2_2 + 0.880354511);
      } break;
      }
   } break;
   }
   return TMath::Sqrt((e1 + e2) * (e1 + e2) - ps);
}
