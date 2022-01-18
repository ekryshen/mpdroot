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
   TVector3 fMomFirstDaughter;
   TVector3 fMomSecondDaughter;
   Double_t fDau1to2;
   Double_t fAplhaArm;
   Double_t fPtArm;
   Double_t fCosAngle;
   Double_t fDecLenght;
   Int_t    fFirstDaugherIndex;
   Int_t    fSecondDaughterIndex;

public:
   MpdV0Particle();
   MpdV0Particle(const MpdV0Particle &other) = default;
   MpdV0Particle &operator=(const MpdV0Particle &other) = default;
   virtual ~MpdV0Particle();

   Int_t GetFirstDaughterndex() const { return fFirstDaugherIndex; }
   Int_t GetSecondDaugterIndex() const { return fSecondDaughterIndex; }

   Double_t GetDau1to2() const { return fDau1to2; }
   Double_t GetAplhaArm() const { return fAplhaArm; }
   Double_t GetPtArm() const { return fPtArm; }
   Double_t GetCosAngle() const { return fCosAngle; }
   Double_t GetDecayLenght() const { return fDecLenght; }

   const TVector3 &GetDca() const { return fDca; }
   const TVector3 &GetDecayPoint() const { return fDecayPoint; }
   const TVector3 &GetMomFirstDaughter() const { return fMomFirstDaughter; }
   const TVector3 &GetMomSecondDautgher() const { return fMomSecondDaughter; }
   const TVector3 &GetMomentum() const { return fMomentum; }
   const Double_t  GetMinv(MpdCommonV0::EParticleType type) const;
   Double_t        GetLambdaMass() const { return GetMinv(MpdCommonV0::EParticleType::kLambda); };
   Double_t        GetAntiLambdaMass() const { return GetMinv(MpdCommonV0::EParticleType::kAntiLambda); };
   Double_t        GetK0Mass() const { return GetMinv(MpdCommonV0::EParticleType::k0Short); };

   void SetAlphaArm(Double_t alpha) { fAplhaArm = alpha; }
   void SetDau1to2(Double_t dau1to2) { fDau1to2 = dau1to2; }
   void SetFirstDaughterIndex(Int_t dauId) { fFirstDaugherIndex = dauId; }
   void SetSecondDaughterIndex(Int_t dauId) { fSecondDaughterIndex = dauId; }
   void SetDecayPoint(const TVector3 &decayPoint) { fDecayPoint = decayPoint; }
   void SetMomFirstDaughter(const TVector3 &mom) { fMomFirstDaughter = mom; }
   void SetMomSecondDaughter(const TVector3 &mom) { fMomSecondDaughter = mom; }
   void SetMomentum(const TVector3 &momentum) { fMomentum = momentum; }
   void SetPtArm(Double_t ptArm) { this->fPtArm = ptArm; }
   void Recalculate(const TVector3 &vertex);
   ClassDef(MpdV0Particle, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0PARTICLE_H_ */
