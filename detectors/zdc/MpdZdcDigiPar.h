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

#ifndef MPDZDCDIGIPAR_H
#define MPDZDCDIGIPAR_H

#include "FairParGenericSet.h"
#include "TObjArray.h"
class MpdZdcDigiPar : public FairParGenericSet {

public:
   MpdZdcDigiPar(const char *name = "MpdZdcDigiPar", const char *title = "ZDC Digi Parameters",
                 const char *context = "TestDefaultContext");
   virtual ~MpdZdcDigiPar(){};

   void   putParams(FairParamList *list);
   Bool_t getParams(FairParamList *list);

   inline Int_t    GetADCBits() { return fADCBits; };
   inline Double_t GetADCResolution() { return fADCResolution; };
   inline Double_t GetEnergyDigiThreshold() { return fEnergyDigiThreshold; };

private:
   Int_t    fADCBits;       // Hardware parameter
   Double_t fADCResolution; // Hardware parameter
   Double_t fEnergyDigiThreshold;

   ClassDef(MpdZdcDigiPar, 1);
};
#endif // MPDZDCDIGIPAR_H
