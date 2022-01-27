/*
 * MpdV0CandidateCutBasicBasic.cxx
 *
 *  Created on: 27 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0CandidateCutBasic.h"

#include <TVector3.h>

#include "KFParticleTopoReconstructor.h"
#include "KFParticleFinder.h"

#include "MpdV0Particle.h"
#include "MpdV0Track.h"

void MpdV0CandidateCutBasic::SetupKF(KFParticleTopoReconstructor *kf) const
{
   kf->GetKFParticleFinder()->SetLCut(fDecayLenght[0]);
}

Bool_t MpdV0CandidateCutBasic::Pass(MpdV0Track &track)
{

   if (track.GetDau1to2() < fDca1to2[0]) return kFALSE;
   if (track.GetDau1to2() > fDca1to2[1]) return kFALSE;
   if (track.GetDecayLenght() < fDecayLenght[0]) return kFALSE;
   if (track.GetDecayLenght() > fDecayLenght[1]) return kFALSE;
   Double_t dca = track.GetDca().Mag();
   if (dca < fDcaPrimVtx[0]) return kFALSE;
   if (dca > fDcaPrimVtx[1]) return kFALSE;
   return kTRUE;
}
