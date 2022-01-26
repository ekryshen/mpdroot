using namespace std;
#include "MpdFfdGeoPar.h"
#include "FairParamList.h"
#include <iostream>
#include <iomanip>
#include "TObjArray.h"

ClassImp(MpdFfdGeoPar)

MpdFfdGeoPar::MpdFfdGeoPar(const char *name, const char *title, const char *context)
   : FairParGenericSet(name, title, context)
{

   fGeoSensNodes = new TObjArray();
   fGeoPassNodes = new TObjArray();
}

MpdFfdGeoPar::~MpdFfdGeoPar(void)
{
   if (fGeoSensNodes) delete fGeoSensNodes;
   if (fGeoPassNodes) delete fGeoPassNodes;
}

void MpdFfdGeoPar::clear(void)
{
   if (fGeoSensNodes) delete fGeoSensNodes;
   if (fGeoPassNodes) delete fGeoPassNodes;
   fGeoSensNodes = fGeoPassNodes = 0x0;
}

void MpdFfdGeoPar::putParams(FairParamList *l)
{
   if (!l) return;
   l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
   l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}

Bool_t MpdFfdGeoPar::getParams(FairParamList *l)
{
   if (!l) return kFALSE;
   if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   if (!l->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;

   return kTRUE;
}
