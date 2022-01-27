/*************************************************************************************
 *
 *         MpdZdcDigiPar
 *    Container class for MpdZdc digitisation parameters
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  04-Apr-2008
 *
 ************************************************************************************/

#include "MpdZdcDigiPar.h"
#include "FairParamList.h"

ClassImp(MpdZdcDigiPar);

MpdZdcDigiPar::MpdZdcDigiPar(const char *name, const char *title, const char *context)
   : FairParGenericSet(name, title, context)
{
}

void MpdZdcDigiPar::putParams(FairParamList *list)
{
   if (!list) return;
   list->add("ADCBits", fADCBits);
   list->add("ADCResolution", fADCResolution);
   list->add("EnergyDigiThreshold", fEnergyDigiThreshold);
}
Bool_t MpdZdcDigiPar::getParams(FairParamList *list)
{
   if (!list) return kFALSE;
   if (!list->fill("ADCBits", &fADCBits)) return kFALSE;
   if (!list->fill("ADCResolution", &fADCResolution)) return kFALSE;
   if (!list->fill("EnergyDigiThreshold", &fEnergyDigiThreshold)) return kFALSE;
   return kTRUE;
}
