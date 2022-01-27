/*
 * MpdV0CandidateCutKalman.h
 *
 *  Created on: 27 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_CANDIDATECUTS_MPDV0CANDIDATECUTKALMAN_H_
#define MPDROOT_PHYSICS_COMMON_V0_CANDIDATECUTS_MPDV0CANDIDATECUTKALMAN_H_

#include <RtypesCore.h>

#include "MpdV0CandidateCutBasic.h"

class MpdV0CandidateCutKalman : public MpdV0CandidateCutBasic {
   Double_t fChiPrim;
   Double_t fChiTopo;
   Double_t fLDL;

public:
   MpdV0CandidateCutKalman() : fChiPrim(1E+6), fChiTopo(1E+6), fLDL(0){};
   void         SetChi2Prim(Double_t val) { fChiPrim = val; };
   void         SetChi2Topo(Double_t val) { fChiTopo = val; };
   void         SetLDL(Double_t val) { fLDL = val; };
   virtual void SetupKF(KFParticleTopoReconstructor *kf) const;
   virtual ~MpdV0CandidateCutKalman(){};
   ClassDef(MpdV0CandidateCutKalman, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_CANDIDATECUTS_MPDV0CANDIDATECUTKALMAN_H_ */
