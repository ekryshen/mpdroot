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

#include <vector>

namespace MpdV0 {
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
std::vector<Int_t> GetPdgs(EParticleType type)
{
   std::vector<Int_t> res;
   switch (type) {
   case EParticleType::k0Short: {
      res.push_back(310); // K0s
      res.push_back(311); // K0
   } break;
   case EParticleType::kLambda: {
      res.push_back(3122);
   } break;
   case EParticleType::kAntiLambda: {
      res.push_back(-3122);
   } break;
   case EParticleType::kPdgHypo: {
      res.push_back(310);   // K0s
      res.push_back(311);   // K0
      res.push_back(3122);  // Lambda
      res.push_back(-3122); // Anti-lambda
   } break;
   }
   return res;
}
} // namespace MpdCommonV0
