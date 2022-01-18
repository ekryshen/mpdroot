/*
 * MpdV0CandicateCut.cxx
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0CandidateCut.h"

#include "MpdV0Particle.h"
#include "MpdV0Track.h"

Bool_t MpdV0CandidateCut::Pass(MpdV0Track &track)
{

   if (track.GetDau1to2() > fDca1to2) return kFALSE;
   if (track.GetDecayLenght() < fDecayLenght) return kFALSE;
   if (track.GetDca().Mag() > fDcaPrimVtx) return kFALSE;
   return kTRUE;
}
