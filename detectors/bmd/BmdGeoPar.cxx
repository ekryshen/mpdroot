using namespace std;
#include "TObjArray.h"
#include "BmdGeoPar.h"
#include "FairParamList.h"
#include <iostream>
#include <iomanip>

ClassImp(BmdGeoPar)

BmdGeoPar::BmdGeoPar(const char* name,const char* title,const char* context)
           : FairParGenericSet(name,title,context)
{
  fGeoSensNodes = new TObjArray();
  fGeoPassNodes = new TObjArray();
}

BmdGeoPar::~BmdGeoPar(void) 
{
  if(fGeoSensNodes) delete fGeoSensNodes;
  if(fGeoPassNodes) delete fGeoPassNodes;
}

void BmdGeoPar::clear(void)
{
  if(fGeoSensNodes) 
    delete fGeoSensNodes;
  if(fGeoPassNodes) 
    delete fGeoPassNodes;
  fGeoSensNodes = fGeoPassNodes = 0x0; //AZ
}

void BmdGeoPar::putParams(FairParamList* l)
{
  if (!l) 
    return;
  l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
  l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}

Bool_t BmdGeoPar::getParams(FairParamList* l)
{
  if (!l) 
    return kFALSE;
   if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes))
    return kFALSE;
   if (!l->fillObject("FairGeoNodes Passive List",   fGeoPassNodes))
    return kFALSE;

  return kTRUE;
}
