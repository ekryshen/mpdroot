/*************************************************************************************
 *
 *         Class MpdZdcGeoPar
 *
 *  Adopted for MPD by:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

using namespace std;
#include "MpdZdcGeoPar.h"
#include "FairParamList.h"
#include <iostream>
#include <iomanip>

ClassImp(MpdZdcGeoPar);

MpdZdcGeoPar::MpdZdcGeoPar(const char *name, const char *title, const char *context)
   : FairParGenericSet(name, title, context)
{
   fGeoSensNodes = new TObjArray();
   fGeoPassNodes = new TObjArray();
}

MpdZdcGeoPar::~MpdZdcGeoPar(void) {}

// probably the next funtions can be deleted

void MpdZdcGeoPar::clear(void)
{
   if (fGeoSensNodes) delete fGeoSensNodes;
   if (fGeoPassNodes) delete fGeoPassNodes;
}

void MpdZdcGeoPar::putParams(FairParamList *l)
{
   if (!l) return;
   l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
   l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}

Bool_t MpdZdcGeoPar::getParams(FairParamList *l)
{
   if (!l) return kFALSE;
   if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
   if (!l->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;

   return kTRUE;
}
