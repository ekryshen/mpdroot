/*
 * MpdV0CandicateCut.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUT_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUT_H_

#include <Rtypes.h>
#include <RtypesCore.h>

#include "MpdV0FinderCut.h"

class KFParticleTopoReconstructor;

class MpdV0Track;

/**
 * abstract class for rejecting bad V0 candidates
 */
class MpdV0CandidateCut : public MpdV0FinderCut {

public:
   MpdV0CandidateCut(){};
   virtual void SetupKF(KFParticleTopoReconstructor *kf) const {};
   virtual ~MpdV0CandidateCut(){};
   virtual Bool_t Pass(MpdV0Track &track)            = 0;
   MpdV0CandidateCut(const MpdV0CandidateCut &other) = default;
   MpdV0CandidateCut &operator=(const MpdV0CandidateCut &other) = default;
   ClassDef(MpdV0CandidateCut, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUT_H_ */
