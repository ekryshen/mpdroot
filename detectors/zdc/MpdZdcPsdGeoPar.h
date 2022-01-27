/*************************************************************************************
 *
 *         Class MpdZdcPsdGeoPar
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

#ifndef MPDGEOZDCPSDPAR_H
#define MPDGEOZDCPSDPAR_H

#include "MpdZdcGeoPar.h"

class MpdZdcPsdGeoPar : public MpdZdcGeoPar {
public:
   MpdZdcPsdGeoPar(const char *name = "MpdZdcPsdGeoPar", const char *title = "ZdcPsd Geometry Parameters",
                   const char *context = "ZdcPsdDefaultContext")
      : MpdZdcGeoPar(name, title, context){};
   ~MpdZdcPsdGeoPar(void){};

   ClassDef(MpdZdcPsdGeoPar, 1);
};

#endif /* !MPDGEOZDCPSDPAR_H */
