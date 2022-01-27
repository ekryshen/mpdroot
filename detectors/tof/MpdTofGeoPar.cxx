//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <iomanip>

#include "TObjArray.h"

#include "FairParamList.h"
#include "MpdTofGeoPar.h"

using namespace std;
//------------------------------------------------------------------------------------------------------------------------
MpdTofGeoPar::MpdTofGeoPar(const char *name, const char *title, const char *context)
   : FairParGenericSet(name, title, context)
{
   fGeoSensNodes = new TObjArray();
   fGeoPassNodes = new TObjArray();
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofGeoPar::~MpdTofGeoPar(void) {}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofGeoPar::clear(void)
{
   if (fGeoSensNodes) delete fGeoSensNodes;
   if (fGeoPassNodes) delete fGeoPassNodes;
   fGeoSensNodes = fGeoPassNodes = 0x0; // AZ
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofGeoPar::putParams(FairParamList *list)
{
   if (!list) return;
   // list->addBinary("FairGeoNodes Sensitive List", fGeoSensNodes);
   // list->addBinary("FairGeoNodes Passive List", fGeoPassNodes);
   list->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
   list->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdTofGeoPar::getParams(FairParamList *list)
{
   if (!list) return kFALSE;
   // if(!list->fillBinary("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   // if(!list->fillBinary("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;
   if (!list->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   if (!list->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;

   return kTRUE;
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdTofGeoPar);
