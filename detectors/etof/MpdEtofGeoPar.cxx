//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <iomanip>

#include "MpdEtofGeoPar.h"
#include "FairParamList.h"

using namespace std;
//------------------------------------------------------------------------------------------------------------------------
MpdEtofGeoPar::MpdEtofGeoPar(const char *name, const char *title, const char *context)
   : FairParGenericSet(name, title, context)
{
   fGeoSensNodes = new TObjArray();
   fGeoPassNodes = new TObjArray();
}
//------------------------------------------------------------------------------------------------------------------------
MpdEtofGeoPar::~MpdEtofGeoPar(void)
{
   delete fGeoSensNodes;
   delete fGeoPassNodes;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofGeoPar::clear(void)
{
   if (fGeoSensNodes) fGeoSensNodes->Clear();
   if (fGeoPassNodes) fGeoPassNodes->Clear();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofGeoPar::putParams(FairParamList *l)
{
   if (!l) return;
   // l->addBinary("FairGeoNodes Sensitive List", fGeoSensNodes);
   // l->addBinary("FairGeoNodes Passive List", fGeoPassNodes);
   l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
   l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdEtofGeoPar::getParams(FairParamList *l)
{
   if (!l) return kFALSE;
   // if (!l->fillBinary("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   // if (!l->fillBinary("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;
   if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   if (!l->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;

   return kTRUE;
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdEtofGeoPar)
