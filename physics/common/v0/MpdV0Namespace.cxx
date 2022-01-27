/*
 * MpdV0Namespace.cxx
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0Namespace.h"

namespace MpdCommonV0 {
Double_t GetMass(EParticleType type)
{
   switch (type) {
   case EParticleType::k0Short: {
      return 0.49761;
   } break;
   case EParticleType::kLambda: {
      return 1.115683;
   } break;
   case EParticleType::kAntiLambda: {
      return 1.115683;
   } break;
   }
   return 0;
}
} // namespace MpdCommonV0
