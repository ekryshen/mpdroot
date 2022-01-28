/*
 * MpdV0Particle.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0PARTICLE_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0PARTICLE_H_

#include <TObject.h>
#include <TVector3.h>

#include "MpdV0Namespace.h"

class MpdV0Particle : public TObject {
   TVector3 fMomentum;
   TVector3 fDca;
   TVector3 fDecayPoint;
   TVector3 fMomPosDaughter;
   TVector3 fMomNegDaughter;
   Double_t fDau1to2;
   Double_t fAplhaArm;
   Double_t fPtArm;
   Double_t fCosAngle;
   Double_t fDecLenght;
   Int_t    fPosDaugherIndex;
   Int_t    fNegDaughterIndex;
   Int_t    fPdg;

public:
   MpdV0Particle();
   MpdV0Particle(const MpdV0Particle &other) = default;
   MpdV0Particle &operator=(const MpdV0Particle &other) = default;
   virtual ~MpdV0Particle();

   Int_t GetPositiveDaughterIndex() const { return fPosDaugherIndex; }
   Int_t GetNegativeDaugterIndex() const { return fNegDaughterIndex; }
   /**
    * returns assumed pdg
    */
   Int_t GetPdg() const { return fPdg; }

   Double_t GetDau1to2() const { return fDau1to2; }
   Double_t GetAplhaArm() const { return fAplhaArm; }
   Double_t GetPtArm() const { return fPtArm; }
   Double_t GetCosAngle() const { return fCosAngle; }
   Double_t GetDecayLenght() const { return fDecLenght; }

   const TVector3 &GetDca() const { return fDca; }
   const TVector3 &GetDecayPoint() const { return fDecayPoint; }
   const TVector3 &GetMomPositiveDaughter() const { return fMomPosDaughter; }
   const TVector3 &GetMomNegativeDautgher() const { return fMomNegDaughter; }
   const TVector3 &GetMomentum() const { return fMomentum; }
   const Double_t  GetMinv(MpdCommonV0::EParticleType type) const;
   Double_t        GetLambdaMass() const { return GetMinv(MpdCommonV0::EParticleType::kLambda); };
   Double_t        GetAntiLambdaMass() const { return GetMinv(MpdCommonV0::EParticleType::kAntiLambda); };
   Double_t        GetK0Mass() const { return GetMinv(MpdCommonV0::EParticleType::k0Short); };

   void SetPdg(Int_t pid) { fPdg = pid; }
   void SetAlphaArm(Double_t alpha) { fAplhaArm = alpha; }
   void SetDau1to2(Double_t dau1to2) { fDau1to2 = dau1to2; }
   void SetPositiveDaughterIndex(Int_t dauId) { fPosDaugherIndex = dauId; }
   void SetNegativeDaughterIndex(Int_t dauId) { fNegDaughterIndex = dauId; }
   void SetDecayPoint(const TVector3 &decayPoint) { fDecayPoint = decayPoint; }
   void SetMomPositiveDaughter(const TVector3 &mom) { fMomPosDaughter = mom; }
   void SetMomNegativeDaughter(const TVector3 &mom) { fMomNegDaughter = mom; }
   void SetMomentum(const TVector3 &momentum) { fMomentum = momentum; }
   void SetPtArm(Double_t ptArm) { this->fPtArm = ptArm; }
   void Recalculate(const TVector3 &vertex);
   ClassDef(MpdV0Particle, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0PARTICLE_H_ */
