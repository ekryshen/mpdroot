using namespace std;
#include "MpdStsGeoPar.h"
#include "FairParamList.h"
#include <iostream>
#include <iomanip>
#include "TObjArray.h"

ClassImp(MpdStsGeoPar);

MpdStsGeoPar::MpdStsGeoPar(const char *name, const char *title, const char *context)
   : FairParGenericSet(name, title, context)
{

   fGeoSensNodes = new TObjArray();
   fGeoPassNodes = new TObjArray();
}

MpdStsGeoPar::~MpdStsGeoPar(void) {}

void MpdStsGeoPar::clear(void)
{
   if (fGeoSensNodes) delete fGeoSensNodes;
   if (fGeoPassNodes) delete fGeoPassNodes;
   fGeoSensNodes = fGeoPassNodes = 0x0; // AZ
}

void MpdStsGeoPar::putParams(FairParamList *l)
{
   if (!l) return;
   // l->addBinary("FairGeoNodes Sensitive List", fGeoSensNodes);
   // l->addBinary("FairGeoNodes Passive List", fGeoPassNodes);
   l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
   l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}

Bool_t MpdStsGeoPar::getParams(FairParamList *l)
{
   if (!l) return kFALSE;
   // if (!l->fillBinary("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   // if (!l->fillBinary("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;
   if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   if (!l->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;
   return kTRUE;
}
