/*************************************************************************************
 *
 *         Class MpdZdcGeoPar
 *
 *  Adopted for MPD by:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

#ifndef MPDGEOZDCPAR_H
#define MPDGEOZDCPAR_H

#include "FairParGenericSet.h"
#include "TObjArray.h"

class MpdZdcGeoPar : public FairParGenericSet {
public:
   TObjArray *fGeoSensNodes; /** List of FairGeoNodes for sensitive  volumes */
   TObjArray *fGeoPassNodes; /** List of FairGeoNodes for passive  volumes */

   MpdZdcGeoPar(const char *name = "MpdZdcGeoPar", const char *title = "Zdc Geometry Parameters",
                const char *context = "ZdcDefaultContext");
   ~MpdZdcGeoPar(void);
   void       clear(void);
   void       putParams(FairParamList *);
   Bool_t     getParams(FairParamList *);
   TObjArray *GetGeoSensitiveNodes() { return fGeoSensNodes; }
   TObjArray *GetGeoPassiveNodes() { return fGeoPassNodes; }

   ClassDef(MpdZdcGeoPar, 1);
};

#endif /* !MPDGEOZDCPAR_H */
