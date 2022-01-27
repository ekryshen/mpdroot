/*
 * MpdV0CandidateCutBasicBasic.h
 *
 *  Created on: 27 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUTBASIC_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUTBASIC_H_

#include <Rtypes.h>
#include <RtypesCore.h>

#include "MpdV0CandidateCut.h"

class KFParticleTopoReconstructor;
class MpdV0Track;

/**
 * basic class for V0 particle selection
 */
class MpdV0CandidateCutBasic : public MpdV0CandidateCut {
   Double_t fDca1to2[2];
   Double_t fDecayLenght[2];
   Double_t fDcaPrimVtx[2];

public:
   MpdV0CandidateCutBasic() : fDca1to2{0, 2.0}, fDecayLenght{0, 1E+10}, fDcaPrimVtx{0, 1E+10} {};
   void SetDcaBetweenDaughters(Double_t min, Double_t max)
   {
      fDca1to2[0] = min;
      fDca1to2[1] = max;
   }
   void SetDecayLenght(Double_t min, Double_t max)
   {
      fDecayLenght[0] = min;
      fDecayLenght[1] = max;
   }
   void SetMaxDcaPrimVtx(Double_t min, Double_t max)
   {
      fDcaPrimVtx[0] = min;
      fDcaPrimVtx[1] = max;
   };
   virtual void SetupKF(KFParticleTopoReconstructor *kf) const;
   virtual ~MpdV0CandidateCutBasic(){};
   virtual Bool_t Pass(MpdV0Track &track);
   MpdV0CandidateCutBasic(const MpdV0CandidateCutBasic &other) = default;
   MpdV0CandidateCutBasic &operator=(const MpdV0CandidateCutBasic &other) = default;
   ClassDef(MpdV0CandidateCutBasic, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUTBASIC_H_ */
