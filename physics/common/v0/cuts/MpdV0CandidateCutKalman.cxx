/*
 * MpdV0CandidateCutKalman.cxx
 *
 *  Created on: 27 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0CandidateCutKalman.h"

#include "KFParticleFinder.h"
#include "KFParticleTopoReconstructor.h"
#include "MpdV0CandidateCut.h"

void MpdV0CandidateCutKalman::SetupKF(KFParticleTopoReconstructor *kf) const
{
   MpdV0CandidateCut::SetupKF(kf);
   kf->GetKFParticleFinder()->SetChi2Cut2D(fChiTopo);
   kf->GetKFParticleFinder()->SetChiPrimaryCut2D(fChiPrim);
   kf->GetKFParticleFinder()->SetLdLCut2D(fLDL);
}
