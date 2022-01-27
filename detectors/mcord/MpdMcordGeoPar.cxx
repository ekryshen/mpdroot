/*
 * MpdMcordPar.cxx
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdMcordGeoPar.h"

#include "FairParamList.h" // for FairParamList

#include "TObjArray.h"

ClassImp(MpdMcordGeoPar);

MpdMcordGeoPar::~MpdMcordGeoPar()
{
   // TODO Auto-generated destructor stub
}

MpdMcordGeoPar::MpdMcordGeoPar(const char *name, const char *title, const char *context)
   : FairParGenericSet(name, title, context), fGeoSensNodes(new TObjArray()), fGeoPassNodes(new TObjArray())
{
}

void MpdMcordGeoPar::clear(void)
{
   delete fGeoSensNodes;
   delete fGeoPassNodes;
}

void MpdMcordGeoPar::putParams(FairParamList *l)
{
   if (!l) {
      return;
   }
   l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
   l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}

Bool_t MpdMcordGeoPar::getParams(FairParamList *l)
{
   if (!l) {
      return kFALSE;
   }
   if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) {
      return kFALSE;
   }
   if (!l->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) {
      return kFALSE;
   }
   return kTRUE;
}
