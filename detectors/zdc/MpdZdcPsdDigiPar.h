/*************************************************************************************
 *
 *         MpdZdcPsdDigiPar
 *    Container class for MpdZdc digitisation parameters
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  04-Apr-2008
 *
 ************************************************************************************/

#ifndef MPDZDCPSDDIGIPAR_H
#define MPDZDCPSDDIGIPAR_H

#include "FairParGenericSet.h"
#include "TObjArray.h"

class MpdZdcPsdDigiPar : public MpdZdcDigiPar {

public:
   MpdZdcPsdDigiPar(const char *name = "MpdZdcPsdDigiPar", const char *title = "ZDCPSD Digi Parameters",
                    const char *context = "TestDefaultContext")
      : MpdZdcDigiPar(name, title, context){};
   virtual ~MpdZdcPsdDigiPar(){};

   ClassDef(MpdZdcPsdDigiPar, 1);
};
#endif // MPDZDCPSDDIGIPAR_H
